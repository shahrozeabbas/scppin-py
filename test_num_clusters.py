#!/usr/bin/env python
"""
Test the num_clusters parameter with the exposed API.
"""

import pandas as pd
import networkx as nx
import numpy as np
from scppin.core.scoring import compute_node_scores
from scppin.core.pcst_solver import detect_functional_module_core
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
    print("TESTING num_clusters PARAMETER")
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
    
    # Get top prizes for reference
    min_score = min(node_scores.values())
    prizes = {node: score - min_score for node, score in node_scores.items()}
    sorted_by_prize = sorted(prizes.items(), key=lambda x: x[1], reverse=True)
    
    print(f"\nTop 10 prizes for reference:")
    for i, (node, prize) in enumerate(sorted_by_prize[:10], 1):
        print(f"  {i:2d}. {node:15s}: prize={prize:8.4f}")
    
    # Test with different num_clusters
    print(f"\n" + "="*80)
    print("RESULTS WITH DIFFERENT num_clusters")
    print("="*80)
    
    for num_clusters in [1, 2, 3, 4, 5]:
        print(f"\nnum_clusters = {num_clusters}:")
        print("-"*80)
        
        module = detect_functional_module_core(
            network_filtered,
            node_scores,
            edge_weight_attr=None,
            c0=0.01,
            num_clusters=num_clusters
        )
        
        module_nodes = list(module.nodes())
        print(f"Total nodes in module: {len(module_nodes)}")
        print(f"Total edges in module: {module.number_of_edges()}")
        
        # Analyze components
        if len(module_nodes) > 1:
            components = list(nx.connected_components(module))
            print(f"Number of connected components: {len(components)}")
            
            for i, component in enumerate(sorted(components, key=len, reverse=True), 1):
                comp_nodes = list(component)
                comp_prize = sum(prizes[n] for n in comp_nodes)
                comp_edges = module.subgraph(comp_nodes).number_of_edges()
                
                # Find ranks
                comp_nodes_ranked = sorted(comp_nodes, key=lambda n: prizes[n], reverse=True)
                top_rank = next(j for j, (n, p) in enumerate(sorted_by_prize) if n == comp_nodes_ranked[0]) + 1
                
                print(f"\n  Component {i}:")
                print(f"    Nodes: {len(comp_nodes)}, Edges: {comp_edges}")
                print(f"    Total prize: {comp_prize:.4f}")
                print(f"    Top node: {comp_nodes_ranked[0]} (prize={prizes[comp_nodes_ranked[0]]:.4f}, rank {top_rank})")
                
                if len(comp_nodes) <= 10:
                    print(f"    All nodes:")
                    for node in comp_nodes_ranked:
                        prize = prizes[node]
                        rank = next(j for j, (n, p) in enumerate(sorted_by_prize) if n == node) + 1
                        print(f"      - {node:15s}: prize={prize:8.4f} (rank {rank})")
                else:
                    print(f"    Top 5 nodes:")
                    for node in comp_nodes_ranked[:5]:
                        prize = prizes[node]
                        rank = next(j for j, (n, p) in enumerate(sorted_by_prize) if n == node) + 1
                        print(f"      - {node:15s}: prize={prize:8.4f} (rank {rank})")
                    print(f"    ... and {len(comp_nodes) - 5} more")
        else:
            print(f"  Single node: {module_nodes[0]}")
    
    print(f"\n" + "="*80)

if __name__ == '__main__':
    main()

