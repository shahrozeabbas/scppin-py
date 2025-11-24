#!/usr/bin/env python
"""
Comprehensive test comparing num_clusters values 1-5 with FDR=0.1.

Tests across different num_clusters to understand module selection patterns.
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
    print("="*100)
    print("COMPREHENSIVE num_clusters COMPARISON WITH FDR=0.1")
    print("="*100)
    
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
    
    # Use FDR=0.1
    fdr = 0.1
    print(f"\nFDR threshold: {fdr}")
    
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
    
    print(f"Total nodes with scores: {len(node_scores)}")
    print(f"Base cost (mean prize): {np.mean(list(prizes.values())):.6f}")
    
    print(f"\nTop 20 prizes for reference:")
    print(f"{'Rank':<6} {'Node':<15} {'Prize':<12} {'Score':<12} {'P-value':<12}")
    print("-"*70)
    for i, (node, prize) in enumerate(sorted_by_prize[:20], 1):
        score = node_scores[node]
        pval = pvalues_filtered[node]
        print(f"{i:<6} {node:<15} {prize:<12.6f} {score:<12.6f} {pval:<12.6e}")
    
    # Test with different num_clusters
    print(f"\n" + "="*100)
    print("RESULTS BY num_clusters")
    print("="*100)
    
    results_summary = []
    
    for num_clusters in range(1, 6):
        print(f"\n{'='*100}")
        print(f"num_clusters = {num_clusters}")
        print(f"{'='*100}")
        
        module = detect_functional_module_core(
            network_filtered,
            node_scores,
            edge_weight_attr=None,
            c0=0.01,
            num_clusters=num_clusters
        )
        
        module_nodes = list(module.nodes())
        print(f"\nOverall Statistics:")
        print(f"  Total nodes: {len(module_nodes)}")
        print(f"  Total edges: {module.number_of_edges()}")
        total_prize = sum(prizes[n] for n in module_nodes)
        print(f"  Total prize: {total_prize:.6f}")
        
        # Analyze components
        if len(module_nodes) > 0:
            if nx.is_connected(module):
                num_components = 1
            else:
                num_components = nx.number_connected_components(module)
            
            print(f"  Connected components: {num_components}")
            
            components = list(nx.connected_components(module))
            components_sorted = sorted(components, key=lambda c: sum(prizes[n] for n in c), reverse=True)
            
            print(f"\nComponent Details:")
            for i, component in enumerate(components_sorted, 1):
                comp_nodes = list(component)
                comp_prize = sum(prizes[n] for n in comp_nodes)
                comp_edges = module.subgraph(comp_nodes).number_of_edges()
                
                # Find ranks of nodes
                comp_nodes_ranked = sorted(comp_nodes, key=lambda n: prizes[n], reverse=True)
                ranks = [next(j for j, (n, p) in enumerate(sorted_by_prize) if n == nn) + 1 
                        for nn in comp_nodes_ranked]
                
                avg_rank = np.mean(ranks)
                best_rank = min(ranks)
                
                print(f"\n  Component {i}:")
                print(f"    Nodes: {len(comp_nodes)}, Edges: {comp_edges}")
                print(f"    Total prize: {comp_prize:.6f}")
                print(f"    Best node rank: {best_rank}, Avg rank: {avg_rank:.1f}")
                
                print(f"    Top 5 nodes:")
                for j, node in enumerate(comp_nodes_ranked[:5], 1):
                    prize = prizes[node]
                    rank = next(k for k, (n, p) in enumerate(sorted_by_prize) if n == node) + 1
                    score = node_scores[node]
                    pval = pvalues_filtered[node]
                    print(f"      {j}. {node:<15} rank={rank:4d} prize={prize:10.6f} score={score:10.6f} p={pval:.2e}")
                
                if len(comp_nodes) > 5:
                    print(f"    ... and {len(comp_nodes) - 5} more nodes")
                
                results_summary.append({
                    'num_clusters': num_clusters,
                    'component': i,
                    'comp_nodes': len(comp_nodes),
                    'comp_prize': comp_prize,
                    'best_rank': best_rank,
                    'avg_rank': avg_rank,
                    'top_node': comp_nodes_ranked[0]
                })
    
    # Summary comparison
    print(f"\n" + "="*100)
    print("SUMMARY COMPARISON")
    print("="*100)
    
    print(f"\n{'num_clusters':<15} {'Component':<12} {'Nodes':<10} {'Prize':<15} {'Best Rank':<12} {'Avg Rank':<12} {'Top Node':<15}")
    print("-"*95)
    
    for result in results_summary:
        print(f"{result['num_clusters']:<15} {result['component']:<12} {result['comp_nodes']:<10} "
              f"{result['comp_prize']:<15.6f} {result['best_rank']:<12} {result['avg_rank']:<12.1f} "
              f"{result['top_node']:<15}")
    
    # Find best configurations
    print(f"\n" + "="*100)
    print("KEY FINDINGS")
    print("="*100)
    
    # Best single component by prize
    best_by_prize = max(results_summary, key=lambda r: r['comp_prize'])
    print(f"\nBest component by total prize:")
    print(f"  Config: num_clusters={best_by_prize['num_clusters']}, component {best_by_prize['component']}")
    print(f"  Prize: {best_by_prize['comp_prize']:.6f}, Nodes: {best_by_prize['comp_nodes']}")
    print(f"  Best rank: {best_by_prize['best_rank']}, Avg rank: {best_by_prize['avg_rank']:.1f}")
    
    # Best by rank
    best_by_rank = min(results_summary, key=lambda r: r['best_rank'])
    print(f"\nBest component by best node rank:")
    print(f"  Config: num_clusters={best_by_rank['num_clusters']}, component {best_by_rank['component']}")
    print(f"  Best rank: {best_by_rank['best_rank']}, Prize: {best_by_rank['comp_prize']:.6f}")
    print(f"  Nodes: {best_by_rank['comp_nodes']}")
    
    # Trend analysis
    print(f"\nTrend across num_clusters:")
    for nc in range(1, 6):
        nc_results = [r for r in results_summary if r['num_clusters'] == nc]
        if nc_results:
            best_component = max(nc_results, key=lambda r: r['comp_prize'])
            print(f"  num_clusters={nc}: Best rank={best_component['best_rank']:4d}, "
                  f"Prize={best_component['comp_prize']:8.4f}, Nodes={best_component['comp_nodes']:3d}")
    
    print(f"\n" + "="*100)

if __name__ == '__main__':
    main()

