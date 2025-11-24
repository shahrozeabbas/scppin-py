#!/usr/bin/env python
"""
Test script for min-max normalized module detection.

This script tests the new min-max normalization for node prizes
using the real data files from the tests directory.
"""

import pandas as pd
import networkx as nx
import numpy as np
from scppin import scPPIN

def load_network_with_weights(network_file):
    """Load network from TSV file with edge weights."""
    print(f"Loading network from {network_file}...")
    df = pd.read_csv(network_file, sep='\t')
    print(f"  Loaded {len(df)} edges")
    print(f"  Columns: {list(df.columns)}")
    
    # Create network
    G = nx.Graph()
    for _, row in df.iterrows():
        G.add_edge(row['Protein_1'], row['Protein_2'], weight=row['Score'])
    
    print(f"  Network: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")
    
    # Check edge weight range
    weights = [G[u][v]['weight'] for u, v in G.edges()]
    print(f"  Edge weights: min={min(weights):.4f}, max={max(weights):.4f}, mean={np.mean(weights):.4f}")
    
    return G

def load_pvalues(pvalue_file):
    """Load p-values from TSV file."""
    print(f"\nLoading p-values from {pvalue_file}...")
    df = pd.read_csv(pvalue_file, sep='\t')
    print(f"  Loaded {len(df)} genes")
    print(f"  Columns: {list(df.columns)}")
    
    # Create p-value dictionary
    pvalues = dict(zip(df['GENE'], df['P']))
    
    # Check p-value range
    pvals_array = np.array(list(pvalues.values()))
    print(f"  P-values: min={np.min(pvals_array):.6f}, max={np.max(pvals_array):.6f}")
    print(f"  Significant (p<0.05): {np.sum(pvals_array < 0.05)}")
    print(f"  Significant (p<0.01): {np.sum(pvals_array < 0.01)}")
    print(f"  Significant (p<0.001): {np.sum(pvals_array < 0.001)}")
    
    return pvalues

def test_module_detection(network, pvalues, fdr=0.01, edge_weight_scale=1.0):
    """Test module detection with min-max normalized prizes."""
    print(f"\n{'='*70}")
    print(f"Testing Module Detection (FDR={fdr}, edge_weight_scale={edge_weight_scale})")
    print(f"{'='*70}")
    
    # Create analyzer
    analyzer = scPPIN()
    
    # Set network and weights
    analyzer.network = network.copy()
    analyzer.set_node_weights(pvalues)
    
    # Detect module
    print("\nDetecting module...")
    try:
        module = analyzer.detect_module(
            fdr=fdr,
            edge_weight_attr='weight',
            edge_weight_scale=edge_weight_scale
        )
        
        print(f"\n✓ Module detected successfully!")
        print(f"  Nodes: {module.number_of_nodes()}")
        print(f"  Edges: {module.number_of_edges()}")
        
        if module.number_of_nodes() > 0:
            # Get node prizes and scores
            prizes = [module.nodes[node].get('prize', 0) for node in module.nodes()]
            scores = [module.nodes[node].get('score', 0) for node in module.nodes()]
            
            print(f"\n  Node prizes (min-max normalized):")
            print(f"    Min: {min(prizes):.4f}")
            print(f"    Max: {max(prizes):.4f}")
            print(f"    Mean: {np.mean(prizes):.4f}")
            print(f"    Median: {np.median(prizes):.4f}")
            
            print(f"\n  Original node scores:")
            print(f"    Min: {min(scores):.4f}")
            print(f"    Max: {max(scores):.4f}")
            print(f"    Mean: {np.mean(scores):.4f}")
            
            # Show top nodes
            print(f"\n  Top 10 nodes by prize:")
            node_prize_pairs = [(node, module.nodes[node].get('prize', 0)) 
                               for node in module.nodes()]
            node_prize_pairs.sort(key=lambda x: x[1], reverse=True)
            for i, (node, prize) in enumerate(node_prize_pairs[:10], 1):
                score = module.nodes[node].get('score', 0)
                print(f"    {i}. {node}: prize={prize:.4f}, score={score:.4f}")
            
            # Check connectivity
            if module.number_of_nodes() > 1:
                is_connected = nx.is_connected(module)
                print(f"\n  Module is connected: {is_connected}")
                if not is_connected:
                    components = list(nx.connected_components(module))
                    print(f"  Number of components: {len(components)}")
                    for i, comp in enumerate(components, 1):
                        print(f"    Component {i}: {len(comp)} nodes")
        
        return module
        
    except Exception as e:
        print(f"\n✗ Error during module detection: {e}")
        import traceback
        traceback.print_exc()
        return None

def main():
    """Main test function."""
    print("="*70)
    print("Testing Min-Max Normalization for Node Prizes")
    print("="*70)
    
    # Load data
    network_file = 'tests/denoised_mg_ad_fava_network.tsv'
    pvalue_file = 'tests/magma_gene_symbol_results.tsv'
    
    network = load_network_with_weights(network_file)
    pvalues = load_pvalues(pvalue_file)
    
    # Check overlap
    network_genes = set(network.nodes())
    pvalue_genes = set(pvalues.keys())
    overlap = network_genes & pvalue_genes
    print(f"\nGene overlap:")
    print(f"  Network genes: {len(network_genes)}")
    print(f"  P-value genes: {len(pvalue_genes)}")
    print(f"  Overlap: {len(overlap)} ({100*len(overlap)/len(network_genes):.1f}% of network)")
    
    # Test with different FDR thresholds
    print("\n" + "="*70)
    print("Testing with different FDR thresholds")
    print("="*70)
    
    for fdr in [0.05, 0.01, 0.001]:
        module = test_module_detection(network, pvalues, fdr=fdr, edge_weight_scale=1.0)
        if module is not None and module.number_of_nodes() > 0:
            print(f"\n✓ FDR={fdr}: Found module with {module.number_of_nodes()} nodes")
        else:
            print(f"\n✗ FDR={fdr}: No module found or empty module")
    
    # Test with different edge weight scales
    print("\n" + "="*70)
    print("Testing with different edge weight scales (FDR=0.01)")
    print("="*70)
    
    for scale in [0.0, 0.5, 1.0]:
        module = test_module_detection(network, pvalues, fdr=0.01, edge_weight_scale=scale)
        if module is not None and module.number_of_nodes() > 0:
            print(f"\n✓ edge_weight_scale={scale}: Found module with {module.number_of_nodes()} nodes")
        else:
            print(f"\n✗ edge_weight_scale={scale}: No module found or empty module")
    
    print("\n" + "="*70)
    print("Testing Complete!")
    print("="*70)

if __name__ == '__main__':
    main()

