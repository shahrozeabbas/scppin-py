#!/usr/bin/env python
"""
Diagnostic test to understand module structure and connectivity.

Questions:
1. What 3 nodes does PCST select? Are they top 3 by prize?
2. Are those nodes directly connected?
3. Can we reach the next high-prize node from them?
4. What happens with num_clusters > 1?
"""

import pandas as pd
import networkx as nx
import numpy as np
from scppin.core.scoring import compute_node_scores
from scppin.core.pcst_solver import detect_functional_module_core, solve_pcst, prepare_edge_costs
from scppin.core.bum_model import fit_bum

def load_pvalues(filepath):
    df = pd.read_csv(filepath, sep='\t', index_col=0)
    pvalues = dict(zip(df['GENE'], df['P']))
    return pvalues

def load_network_with_weights(filepath):
    df = pd.read_csv(filepath, sep='\t', index_col=0)
    G = nx.Graph()
    for _, row in df.iterrows():
        G.add_edge(row['Protein_1'], row['Protein_2'], weight=row['Score'])
    return G

def main():
    print("="*80)
    print("MODULE STRUCTURE AND CONNECTIVITY DIAGNOSTICS")
    print("="*80)
    
    # Load data
    pvalues = load_pvalues('tests/magma_gene_symbol_results.tsv')
    network = load_network_with_weights('tests/denoised_mg_ad_fava_network.tsv')
    
    genes_in_pvalues = set(pvalues.keys())
    genes_in_network = set(network.nodes())
    overlap = genes_in_pvalues & genes_in_network
    network_filtered = network.subgraph(overlap).copy()
    pvalues_filtered = {k: pvalues[k] for k in overlap}
    
    pvalues_array = np.array(list(pvalues_filtered.values()))
    lambda_param, alpha, success = fit_bum(pvalues_array)
    
    fdr = 0.01
    node_scores = compute_node_scores(
        network_filtered,
        pvalues_filtered,
        lambda_param,
        alpha,
        fdr,
        missing_data_score=False
    )
    
    min_score = min(node_scores.values())
    prizes = {node: score - min_score for node, score in node_scores.items()}
    base_cost = np.mean(list(prizes.values()))
    
    # Get top prizes for reference
    sorted_by_prize = sorted(prizes.items(), key=lambda x: x[1], reverse=True)
    
    print(f"\n[1] TOP 10 PRIZES FOR REFERENCE")
    print("-"*80)
    for i, (node, prize) in enumerate(sorted_by_prize[:10], 1):
        print(f"  {i:2d}. {node:15s}: prize={prize:8.4f}")
    
    # [2] Test with num_clusters = 1 (default)
    print(f"\n[2] SOLVING WITH num_clusters=1 (CURRENT DEFAULT)")
    print("-"*80)
    
    module_1 = detect_functional_module_core(
        network_filtered,
        node_scores,
        edge_weight_attr=None,
        c0=0.01
    )
    
    module_1_nodes = list(module_1.nodes())
    print(f"Module nodes: {module_1_nodes}")
    print(f"Module size: {len(module_1_nodes)} nodes")
    print(f"Module edges: {module_1.number_of_edges()}")
    
    print(f"\nNode details:")
    for node in module_1_nodes:
        prize = prizes.get(node, 0)
        rank = next(i for i, (n, p) in enumerate(sorted_by_prize) if n == node) + 1
        print(f"  {node:15s}: prize={prize:8.4f}, rank={rank}")
    
    # Check connectivity within module
    print(f"\nConnectivity within module:")
    if len(module_1_nodes) > 1:
        module_subgraph = network_filtered.subgraph(module_1_nodes)
        print(f"  Is connected: {nx.is_connected(module_subgraph)}")
        
        # Print edges
        for u, v in module_subgraph.edges():
            weight = network_filtered[u][v]['weight']
            print(f"    {u} -- {v} (weight={weight:.4f})")
    
    # Check nearest high-prize nodes
    print(f"\nNearest high-prize nodes to this module:")
    other_high_prize = [n for n, p in sorted_by_prize[:20] if n not in module_1_nodes]
    
    for node in other_high_prize[:5]:
        prize = prizes[node]
        
        # Find shortest path from any module node to this node
        min_distance = float('inf')
        closest_module_node = None
        
        for module_node in module_1_nodes:
            try:
                distance = nx.shortest_path_length(network_filtered, module_node, node)
                if distance < min_distance:
                    min_distance = distance
                    closest_module_node = module_node
            except nx.NetworkXNoPath:
                pass
        
        if min_distance != float('inf'):
            print(f"  {node:15s} (prize={prize:8.4f}): distance={min_distance} from {closest_module_node}")
        else:
            print(f"  {node:15s} (prize={prize:8.4f}): NO PATH")
    
    # [3] Test with multiple clusters
    print(f"\n[3] SOLVING WITH MULTIPLE CLUSTERS")
    print("-"*80)
    
    for num_clusters in [1, 2, 3, 4, 5]:
        print(f"\nnum_clusters = {num_clusters}:")
        
        # Manually create edge costs and call solve_pcst
        edge_costs = prepare_edge_costs(
            network_filtered,
            base_cost,
            edge_weight_attr=None,
            c0=0.01
        )
        
        solution_nodes = solve_pcst(
            network_filtered,
            prizes,
            edge_costs,
            num_clusters=num_clusters
        )
        
        print(f"  Total nodes in solution: {len(solution_nodes)}")
        
        # Group into connected components
        solution_subgraph = network_filtered.subgraph(solution_nodes)
        components = list(nx.connected_components(solution_subgraph))
        components_sorted = sorted(components, key=len, reverse=True)
        
        for i, component in enumerate(components_sorted, 1):
            comp_nodes = list(component)
            comp_prize = sum(prizes[n] for n in comp_nodes)
            comp_edges = network_filtered.subgraph(comp_nodes).number_of_edges()
            
            print(f"    Component {i}: {len(comp_nodes)} nodes, {comp_edges} edges, prize={comp_prize:.4f}")
            
            # Show nodes
            comp_nodes_ranked = sorted(comp_nodes, key=lambda n: prizes[n], reverse=True)
            for node in comp_nodes_ranked[:5]:
                prize = prizes[node]
                rank = next(j for j, (n, p) in enumerate(sorted_by_prize) if n == node) + 1
                print(f"      - {node:15s}: prize={prize:8.4f} (rank {rank})")
            
            if len(comp_nodes_ranked) > 5:
                print(f"      ... and {len(comp_nodes_ranked) - 5} more")
    
    # [4] Analysis: why stop at 3?
    print(f"\n[4] ANALYSIS: WHY DOES PCST STOP AT 3 NODES?")
    print("-"*80)
    
    print(f"\nManual exploration of module growth:")
    print(f"\nStarting with top node: {sorted_by_prize[0][0]} (prize={sorted_by_prize[0][1]:.4f})")
    
    # Can we connect top 1 to top 2?
    top_1 = sorted_by_prize[0][0]
    top_2 = sorted_by_prize[1][0]
    
    try:
        path_1_2 = nx.shortest_path(network_filtered, top_1, top_2)
        print(f"\nTop 1 to Top 2: shortest path length = {len(path_1_2) - 1}")
        print(f"  Path: {' -> '.join(path_1_2)}")
    except nx.NetworkXNoPath:
        print(f"\nTop 1 to Top 2: NO PATH")
    
    # Can we connect top 1 to top 3?
    top_3 = sorted_by_prize[2][0]
    try:
        path_1_3 = nx.shortest_path(network_filtered, top_1, top_3)
        print(f"\nTop 1 to Top 3: shortest path length = {len(path_1_3) - 1}")
        print(f"  Path: {' -> '.join(path_1_3)}")
    except nx.NetworkXNoPath:
        print(f"\nTop 1 to Top 3: NO PATH")
    
    # Can we connect top 2 to top 3?
    try:
        path_2_3 = nx.shortest_path(network_filtered, top_2, top_3)
        print(f"\nTop 2 to Top 3: shortest path length = {len(path_2_3) - 1}")
        print(f"  Path: {' -> '.join(path_2_3)}")
    except nx.NetworkXNoPath:
        print(f"\nTop 2 to Top 3: NO PATH")
    
    # Can we connect any of top 3 to top 4?
    top_4 = sorted_by_prize[3][0]
    
    print(f"\nConnecting top 3 to Top 4 ({top_4}):")
    for i, (node, _) in enumerate(sorted_by_prize[:3], 1):
        try:
            path = nx.shortest_path(network_filtered, node, top_4)
            print(f"  Top {i} ({node}) to Top 4: path length = {len(path) - 1}")
        except nx.NetworkXNoPath:
            print(f"  Top {i} ({node}) to Top 4: NO PATH")
    
    print("\n" + "="*80)

if __name__ == '__main__':
    main()

