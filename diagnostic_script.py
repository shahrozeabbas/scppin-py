"""
Diagnostic script for scPPIN module detection.
Analyzes why modules are small and identifies potential issues.
"""

import pandas as pd
import numpy as np
import networkx as nx
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from scppin import scPPIN
from scppin.core import fit_bum, compute_node_scores, filter_network_by_pvalues
from scppin.core.scoring import shift_scores_for_pcst, get_minimum_score
from scppin.core.pcst_solver import prepare_edge_costs

def print_section(title):
    """Print a formatted section header."""
    print("\n" + "=" * 80)
    print(f"  {title}")
    print("=" * 80)

def print_subsection(title):
    """Print a formatted subsection header."""
    print(f"\n--- {title} ---")

def load_data():
    """Load network and p-values from test files."""
    print_section("LOADING DATA")
    
    # Load network
    print("\nLoading network...")
    network_file = 'tests/denoised_mg_ad_fava_network.tsv'
    network_df = pd.read_csv(network_file, sep='\t')
    print(f"  Network file: {network_file}")
    print(f"  Total edges in file: {len(network_df)}")
    print(f"  Columns: {list(network_df.columns)}")
    
    # Load p-values
    print("\nLoading p-values...")
    pvalues_file = 'tests/magma_gene_symbol_results.tsv'
    pvalues_df = pd.read_csv(pvalues_file, sep='\t')
    print(f"  P-values file: {pvalues_file}")
    print(f"  Total genes in file: {len(pvalues_df)}")
    print(f"  Columns: {list(pvalues_df.columns)}")
    
    # Extract p-values dictionary
    pvalues = dict(zip(pvalues_df['GENE'], pvalues_df['P']))
    print(f"  Genes with p-values: {len(pvalues)}")
    
    # Check for invalid p-values
    invalid_pvals = {g: p for g, p in pvalues.items() if p <= 0 or p > 1 or np.isnan(p) or np.isinf(p)}
    if invalid_pvals:
        print(f"  WARNING: {len(invalid_pvals)} invalid p-values found!")
        print(f"  Sample invalid p-values: {dict(list(invalid_pvals.items())[:5])}")
        # Filter out invalid p-values
        pvalues = {g: p for g, p in pvalues.items() if g not in invalid_pvals}
    
    return network_df, pvalues

def analyze_network_connectivity(network):
    """Analyze network connectivity metrics."""
    print_section("NETWORK CONNECTIVITY ANALYSIS")
    
    print(f"\nBasic network metrics:")
    print(f"  Nodes: {network.number_of_nodes()}")
    print(f"  Edges: {network.number_of_edges()}")
    print(f"  Density: {nx.density(network):.6f}")
    print(f"  Average degree: {2 * network.number_of_edges() / network.number_of_nodes():.4f}")
    
    # Check if graph is connected
    is_connected = nx.is_connected(network)
    print(f"\nConnectivity:")
    print(f"  Is connected: {is_connected}")
    
    # Analyze connected components
    components = list(nx.connected_components(network))
    component_sizes = [len(c) for c in components]
    component_sizes.sort(reverse=True)
    
    print(f"  Number of connected components: {len(components)}")
    print(f"  Largest component size: {component_sizes[0] if component_sizes else 0}")
    print(f"  Second largest component size: {component_sizes[1] if len(component_sizes) > 1 else 0}")
    print(f"  Smallest component size: {component_sizes[-1] if component_sizes else 0}")
    
    # Component size distribution
    print(f"\nComponent size distribution:")
    print(f"  Components with 1 node (isolated): {sum(1 for s in component_sizes if s == 1)}")
    print(f"  Components with 2 nodes: {sum(1 for s in component_sizes if s == 2)}")
    print(f"  Components with 3-10 nodes: {sum(1 for s in component_sizes if 3 <= s <= 10)}")
    print(f"  Components with >10 nodes: {sum(1 for s in component_sizes if s > 10)}")
    
    if len(component_sizes) <= 20:
        print(f"\n  All component sizes: {component_sizes}")
    else:
        print(f"\n  Top 10 component sizes: {component_sizes[:10]}")
        print(f"  Bottom 10 component sizes: {component_sizes[-10:]}")
    
    # Analyze largest component
    if component_sizes:
        largest_component = max(components, key=len)
        largest_subgraph = network.subgraph(largest_component)
        print(f"\nLargest component details:")
        print(f"  Nodes: {largest_subgraph.number_of_nodes()}")
        print(f"  Edges: {largest_subgraph.number_of_edges()}")
        print(f"  Density: {nx.density(largest_subgraph):.6f}")
        print(f"  Average degree: {2 * largest_subgraph.number_of_edges() / largest_subgraph.number_of_nodes():.4f}")
        
        # Check if largest component is a tree
        is_tree = largest_subgraph.number_of_edges() == largest_subgraph.number_of_nodes() - 1
        print(f"  Is tree: {is_tree}")
        
        # Average path length (if not too large)
        if largest_subgraph.number_of_nodes() <= 1000:
            try:
                avg_path_length = nx.average_shortest_path_length(largest_subgraph)
                print(f"  Average path length: {avg_path_length:.4f}")
            except:
                print(f"  Average path length: N/A (disconnected or too large)")
        
        # Diameter (if not too large)
        if largest_subgraph.number_of_nodes() <= 1000:
            try:
                diameter = nx.diameter(largest_subgraph)
                print(f"  Diameter: {diameter}")
            except:
                print(f"  Diameter: N/A (disconnected or too large)")
    
    return components, component_sizes

def analyze_network(network_df, pvalues):
    """Analyze network structure and gene overlap."""
    print_section("NETWORK ANALYSIS")
    
    # Get unique genes in network
    network_genes = set(network_df['Protein_1'].unique()) | set(network_df['Protein_2'].unique())
    print(f"\nNetwork genes:")
    print(f"  Unique genes in network: {len(network_genes)}")
    print(f"  Sample network genes: {list(network_genes)[:10]}")
    
    # Get genes with p-values
    pvalue_genes = set(pvalues.keys())
    print(f"\nP-value genes:")
    print(f"  Genes with p-values: {len(pvalue_genes)}")
    print(f"  Sample p-value genes: {list(pvalue_genes)[:10]}")
    
    # Check overlap
    overlap = network_genes & pvalue_genes
    print(f"\nOverlap:")
    print(f"  Genes in both network and p-values: {len(overlap)}")
    print(f"  Network-only genes: {len(network_genes - pvalue_genes)}")
    print(f"  P-value-only genes: {len(pvalue_genes - network_genes)}")
    
    if len(overlap) == 0:
        print("\n  ERROR: No overlap between network and p-value genes!")
        print("  This could be due to:")
        print("    - Case sensitivity (e.g., 'ST6GAL1' vs 'st6gal1')")
        print("    - Different naming conventions")
        print("    - Missing genes in one dataset")
        
        # Check case sensitivity
        network_genes_upper = {g.upper() for g in network_genes}
        pvalue_genes_upper = {g.upper() for g in pvalue_genes}
        overlap_upper = network_genes_upper & pvalue_genes_upper
        print(f"\n  Case-insensitive overlap: {len(overlap_upper)}")
        if len(overlap_upper) > 0:
            print("  SUGGESTION: Gene names may need case normalization")
    
    return overlap

def analyze_scoring(network, pvalues, fdr=0.01):
    """Analyze node scoring and prize calculation."""
    print_section("SCORING ANALYSIS")
    
    # Filter network to genes with p-values
    print("\nFiltering network to genes with p-values...")
    network_filtered = filter_network_by_pvalues(network, pvalues, missing_data_score=False)
    print(f"  Nodes before filtering: {network.number_of_nodes()}")
    print(f"  Nodes after filtering: {network_filtered.number_of_nodes()}")
    print(f"  Edges before filtering: {network.number_of_edges()}")
    print(f"  Edges after filtering: {network_filtered.number_of_edges()}")
    
    # Analyze network connectivity AFTER filtering
    print("\n" + "="*80)
    print("  NETWORK CONNECTIVITY AFTER FILTERING")
    print("="*80)
    components_after, component_sizes_after = analyze_network_connectivity(network_filtered)
    
    if network_filtered.number_of_nodes() == 0:
        print("\n  ERROR: No nodes remain after filtering!")
        return None, None, None, None
    
    # Fit BUM model
    print("\nFitting BUM model...")
    pvalues_array = np.array(list(pvalues.values()))
    lambda_param, alpha, success = fit_bum(pvalues_array)
    print(f"  Lambda: {lambda_param:.6f}")
    print(f"  Alpha: {alpha:.6f}")
    print(f"  Convergence: {success}")
    
    if not success:
        print("  WARNING: BUM model did not converge!")
    
    # Compute node scores
    print("\nComputing node scores...")
    node_scores = compute_node_scores(
        network_filtered,
        pvalues,
        lambda_param,
        alpha,
        fdr,
        missing_data_score=False
    )
    
    print(f"  Nodes with scores: {len(node_scores)}")
    if node_scores:
        scores_array = np.array(list(node_scores.values()))
        print(f"  Min score: {np.min(scores_array):.6f}")
        print(f"  Max score: {np.max(scores_array):.6f}")
        print(f"  Mean score: {np.mean(scores_array):.6f}")
        print(f"  Median score: {np.median(scores_array):.6f}")
        print(f"  Std score: {np.std(scores_array):.6f}")
        
        # Count positive/negative scores
        positive_scores = sum(1 for s in scores_array if s > 0)
        negative_scores = sum(1 for s in scores_array if s < 0)
        zero_scores = sum(1 for s in scores_array if s == 0)
        print(f"  Positive scores: {positive_scores}")
        print(f"  Negative scores: {negative_scores}")
        print(f"  Zero scores: {zero_scores}")
        
        # Show top/bottom scores
        sorted_scores = sorted(node_scores.items(), key=lambda x: x[1], reverse=True)
        print(f"\n  Top 10 scores:")
        for gene, score in sorted_scores[:10]:
            pval = pvalues.get(gene, 'N/A')
            print(f"    {gene}: score={score:.6f}, p-value={pval}")
        
        print(f"\n  Bottom 10 scores:")
        for gene, score in sorted_scores[-10:]:
            pval = pvalues.get(gene, 'N/A')
            print(f"    {gene}: score={score:.6f}, p-value={pval}")
    
    # Shift scores to prizes
    print("\nShifting scores to prizes...")
    node_prizes = shift_scores_for_pcst(node_scores)
    min_score = get_minimum_score(node_scores)
    base_cost = -min_score
    
    print(f"  Min score (before shift): {min_score:.6f}")
    print(f"  Base cost: {base_cost:.6f}")
    
    if node_prizes:
        prizes_array = np.array(list(node_prizes.values()))
        print(f"  Min prize: {np.min(prizes_array):.6f}")
        print(f"  Max prize: {np.max(prizes_array):.6f}")
        print(f"  Mean prize: {np.mean(prizes_array):.6f}")
        print(f"  Median prize: {np.median(prizes_array):.6f}")
        
        # Count zero prizes
        zero_prizes = sum(1 for p in prizes_array if p == 0)
        non_zero_prizes = sum(1 for p in prizes_array if p > 0)
        print(f"  Zero prizes: {zero_prizes}")
        print(f"  Non-zero prizes: {non_zero_prizes}")
        
        # Show top prizes
        sorted_prizes = sorted(node_prizes.items(), key=lambda x: x[1], reverse=True)
        print(f"\n  Top 10 prizes:")
        for gene, prize in sorted_prizes[:10]:
            score = node_scores.get(gene, 'N/A')
            print(f"    {gene}: prize={prize:.6f}, score={score}")
    
    return network_filtered, node_scores, node_prizes, base_cost

def analyze_edge_costs(network, base_cost, edge_weight_attr=None, node_prizes=None):
    """Analyze edge costs."""
    print_section("EDGE COST ANALYSIS")
    
    # If node_prizes provided, use mean prize as base_cost (matching new implementation)
    if node_prizes is not None:
        prizes_array = np.array(list(node_prizes.values()))
        actual_base_cost = float(np.mean(prizes_array))
        print(f"\nNote: Using mean prize ({actual_base_cost:.6f}) as base_cost (new implementation)")
        print(f"  Old base_cost would have been: {base_cost:.6f}")
        base_cost = actual_base_cost
    
    edge_costs = prepare_edge_costs(
        network,
        base_cost,
        edge_weight_attr=edge_weight_attr,
        edge_weight_scale=1.0,
        c0=None
    )
    
    if edge_costs:
        costs_array = np.array(list(edge_costs.values()))
        print(f"\nEdge costs:")
        print(f"  Total edges: {len(edge_costs)}")
        print(f"  Min cost: {np.min(costs_array):.6f}")
        print(f"  Max cost: {np.max(costs_array):.6f}")
        print(f"  Mean cost: {np.mean(costs_array):.6f}")
        print(f"  Median cost: {np.median(costs_array):.6f}")
        print(f"  Base cost: {base_cost:.6f}")
    
    return edge_costs

def analyze_prize_cost_balance(node_prizes, edge_costs):
    """Analyze the balance between prizes and costs."""
    print_section("PRIZE-COST BALANCE ANALYSIS")
    
    if not node_prizes or not edge_costs:
        print("  Cannot analyze: missing prizes or costs")
        return
    
    avg_prize = np.mean(list(node_prizes.values()))
    avg_cost = np.mean(list(edge_costs.values()))
    max_prize = max(node_prizes.values())
    
    print(f"\nBalance metrics:")
    print(f"  Average prize: {avg_prize:.6f}")
    print(f"  Average edge cost: {avg_cost:.6f}")
    print(f"  Max prize: {max_prize:.6f}")
    print(f"  Prize/Cost ratio (avg): {avg_prize/avg_cost:.6f}")
    print(f"  Max prize / Avg cost: {max_prize/avg_cost:.6f}")
    
    # Analyze how many nodes can "pay" for edges
    print(f"\nNode-edge economics:")
    nodes_above_avg_cost = sum(1 for p in node_prizes.values() if p > avg_cost)
    nodes_above_2x_cost = sum(1 for p in node_prizes.values() if p > 2 * avg_cost)
    print(f"  Nodes with prize > avg cost: {nodes_above_avg_cost}")
    print(f"  Nodes with prize > 2x avg cost: {nodes_above_2x_cost}")
    
    if avg_prize < avg_cost:
        print("\n  WARNING: Average prize is less than average edge cost!")
        print("  This means most nodes cannot 'pay' for their edges.")
        print("  The PCST solver will likely only include high-prize nodes.")
    elif avg_prize == avg_cost:
        print("\n  ✓ Average prize equals average edge cost (good balance)")
    else:
        print("\n  ✓ Average prize is greater than average edge cost (good balance)")
    
    if max_prize < avg_cost:
        print("\n  ERROR: Even the highest prize is less than average edge cost!")
        print("  The PCST solver will struggle to form any module.")
        print("  SUGGESTIONS:")
        print("    - Lower FDR threshold (e.g., 0.05 or 0.1)")
        print("    - Check if p-values are too large (not significant)")
        print("    - Verify BUM model parameters")

def analyze_high_prize_connectivity(network, node_prizes, avg_cost):
    """Analyze connectivity of high-prize nodes."""
    print_section("HIGH-PRIZE NODE CONNECTIVITY ANALYSIS")
    
    # Get high-prize nodes
    high_prize_nodes = {node: prize for node, prize in node_prizes.items() if prize > avg_cost}
    print(f"\nHigh-prize nodes (prize > avg cost {avg_cost:.6f}): {len(high_prize_nodes)}")
    
    if len(high_prize_nodes) == 0:
        print("  No high-prize nodes found!")
        return
    
    # Check connectivity
    print(f"\nConnectivity analysis:")
    isolated_nodes = []
    connected_nodes = []
    
    for node, prize in sorted(high_prize_nodes.items(), key=lambda x: x[1], reverse=True)[:20]:
        if node in network:
            neighbors = list(network.neighbors(node))
            degree = len(neighbors)
            
            if degree == 0:
                isolated_nodes.append((node, prize))
            else:
                # Check if neighbors also have high prizes
                high_prize_neighbors = [n for n in neighbors if n in high_prize_nodes]
                connected_nodes.append((node, prize, degree, len(high_prize_neighbors)))
    
    if isolated_nodes:
        print(f"\n  Isolated high-prize nodes (no edges): {len(isolated_nodes)}")
        for node, prize in isolated_nodes[:10]:
            print(f"    {node}: prize={prize:.6f}")
    
    if connected_nodes:
        print(f"\n  Connected high-prize nodes (top 10):")
        for node, prize, degree, hp_neighbors in connected_nodes[:10]:
            print(f"    {node}: prize={prize:.6f}, degree={degree}, high-prize neighbors={hp_neighbors}")
    
    # Check if high-prize nodes form a connected component
    high_prize_subgraph = network.subgraph(high_prize_nodes.keys())
    if high_prize_subgraph.number_of_nodes() > 0:
        components = list(nx.connected_components(high_prize_subgraph))
        print(f"\n  High-prize nodes form {len(components)} connected component(s)")
        print(f"  Largest component: {len(max(components, key=len))} nodes")
        if len(components) > 1:
            print(f"  Component sizes: {[len(c) for c in sorted(components, key=len, reverse=True)[:5]]}")

def run_module_detection(network_df, pvalues, fdr=0.01, use_edge_weights=False):
    """Run full module detection and analyze results."""
    if use_edge_weights:
        print_section("MODULE DETECTION (WITH EDGE WEIGHTS)")
    else:
        print_section("MODULE DETECTION")
    
    # Create analyzer
    analyzer = scPPIN()
    
    # Load network
    print("\nLoading network into analyzer...")
    if use_edge_weights:
        # Load with edge weights from Score column
        # Create a DataFrame with edges and weights
        edges_df = network_df[['Protein_1', 'Protein_2', 'Score']].copy()
        edges_df.columns = ['source', 'target', 'weight']
        analyzer.load_network(edges_df, weight_column='weight')
        print(f"  Network loaded with edge weights: {analyzer.network.number_of_nodes()} nodes, {analyzer.network.number_of_edges()} edges")
        print(f"  Edge weights set: {len(analyzer.edge_weights)} edges")
    else:
        edges = [(row['Protein_1'], row['Protein_2']) for _, row in network_df.iterrows()]
        analyzer.load_network(edges)
        print(f"  Network loaded: {analyzer.network.number_of_nodes()} nodes, {analyzer.network.number_of_edges()} edges")
    
    # Set node weights
    print("\nSetting node weights...")
    analyzer.set_node_weights(pvalues)
    print(f"  Node weights set: {len(analyzer.node_weights)} genes")
    print(f"  Network after filtering: {analyzer.network.number_of_nodes()} nodes, {analyzer.network.number_of_edges()} edges")
    
    # Detect module
    print(f"\nDetecting module with FDR={fdr}...")
    try:
        if use_edge_weights:
            module = analyzer.detect_module(fdr=fdr, edge_weight_attr='weight', edge_weight_scale=0.5)
        else:
            module = analyzer.detect_module(fdr=fdr)
        
        print(f"\nModule detected:")
        print(f"  Nodes: {module.number_of_nodes()}")
        print(f"  Edges: {module.number_of_edges()}")
        
        if module.number_of_nodes() > 0:
            print(f"  Genes in module: {sorted(list(module.nodes()))}")
            
            # Show node attributes
            print(f"\n  Node details:")
            for node in module.nodes():
                score = module.nodes[node].get('score', 'N/A')
                prize = module.nodes[node].get('prize', 'N/A')
                pval = pvalues.get(node, 'N/A')
                print(f"    {node}: score={score}, prize={prize}, p-value={pval}")
        else:
            print("  WARNING: Empty module detected!")
            
    except Exception as e:
        print(f"\n  ERROR during module detection: {e}")
        import traceback
        traceback.print_exc()
        module = None
    
    return module

def main():
    """Run full diagnostic analysis."""
    print("\n" + "=" * 80)
    print("  scPPIN DIAGNOSTIC SCRIPT")
    print("=" * 80)
    
    # Load data
    network_df, pvalues = load_data()
    
    # Analyze network
    overlap = analyze_network(network_df, pvalues)
    
    if len(overlap) == 0:
        print("\n\nSTOPPING: No overlap between network and p-value genes.")
        print("Please check gene name matching.")
        return
    
    # Build network
    print_section("BUILDING NETWORK")
    edges = [(row['Protein_1'], row['Protein_2']) for _, row in network_df.iterrows()]
    network = nx.Graph(edges)
    print(f"  Network built: {network.number_of_nodes()} nodes, {network.number_of_edges()} edges")
    
    # Analyze network connectivity BEFORE filtering
    components_before, component_sizes_before = analyze_network_connectivity(network)
    
    # Analyze scoring
    network_filtered, node_scores, node_prizes, base_cost = analyze_scoring(
        network, pvalues, fdr=0.01
    )
    
    if network_filtered is None:
        print("\n\nSTOPPING: Network filtering resulted in empty network.")
        return
    
    # Analyze edge costs (pass node_prizes to use mean prize as base_cost)
    edge_costs = analyze_edge_costs(network_filtered, base_cost, edge_weight_attr=None, node_prizes=node_prizes)
    
    # Analyze prize-cost balance
    analyze_prize_cost_balance(node_prizes, edge_costs)
    
    # Analyze high-prize node connectivity
    if edge_costs:
        avg_cost = np.mean(list(edge_costs.values()))
        analyze_high_prize_connectivity(network_filtered, node_prizes, avg_cost)
    
    # Run module detection without edge weights
    module = run_module_detection(network_df, pvalues, fdr=0.01, use_edge_weights=False)
    
    # Run module detection with edge weights
    print("\n" + "="*80)
    print("  TESTING WITH EDGE WEIGHTS")
    print("="*80)
    module_weighted = run_module_detection(network_df, pvalues, fdr=0.01, use_edge_weights=True)
    
    # Summary
    print_section("SUMMARY")
    print("\nKey findings:")
    if module and module.number_of_nodes() > 0:
        print(f"  ✓ Module detected (no edge weights): {module.number_of_nodes()} nodes")
    else:
        print(f"  ✗ Module detection failed or returned empty module (no edge weights)")
    
    if module_weighted and module_weighted.number_of_nodes() > 0:
        print(f"  ✓ Module detected (with edge weights): {module_weighted.number_of_nodes()} nodes")
    elif module_weighted:
        print(f"  ✗ Module detection with edge weights returned empty module")
    
    print("\nRecommendations:")
    print("  1. High-prize nodes are poorly connected (63 components)")
    print("     → This limits module size regardless of parameters")
    print("  2. Try using edge weights to reduce costs for high-confidence edges")
    print("  3. Try lowering FDR threshold (e.g., 0.05 or 0.1) to get more high-prize nodes")
    print("  4. Consider if the network structure matches your biological question")
    
    if node_prizes:
        zero_prizes = sum(1 for p in node_prizes.values() if p == 0)
        total_prizes = len(node_prizes)
        if zero_prizes > total_prizes * 0.5:
            print(f"  ⚠ {zero_prizes}/{total_prizes} nodes have zero prizes (may indicate scoring issue)")
    
    if edge_costs:
        avg_cost = np.mean(list(edge_costs.values()))
        if node_prizes:
            avg_prize = np.mean(list(node_prizes.values()))
            if abs(avg_prize - avg_cost) < 1e-6:
                print(f"  ✓ Average prize ({avg_prize:.6f}) = average cost ({avg_cost:.6f})")
                print(f"     Balance is good. If modules are still small, check network connectivity.")
            elif avg_prize < avg_cost:
                print(f"  ⚠ Average prize ({avg_prize:.6f}) < average cost ({avg_cost:.6f})")
                print(f"     This may cause small modules. Try lowering FDR threshold.")
            else:
                print(f"  ✓ Average prize ({avg_prize:.6f}) > average cost ({avg_cost:.6f})")
                print(f"     Balance is good. If modules are still small, check network connectivity.")
    
    print("\n" + "=" * 80)

if __name__ == '__main__':
    main()

