#!/usr/bin/env python
"""
Test different combinations of pruning methods and num_clusters.
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

# Load data
pvalues = load_pvalues('tests/magma_gene_symbol_results.tsv')
network = load_network_with_weights('tests/denoised_mg_ad_fava_network.tsv')

genes_in_pvalues = set(pvalues.keys())
genes_in_network = set(network.nodes())
overlap = genes_in_pvalues & genes_in_network
network_filtered = network.subgraph(overlap).copy()
pvalues_filtered = {k: pvalues[k] for k in overlap}

# Fit BUM
pvalues_array = np.array(list(pvalues_filtered.values()))
lambda_param, alpha, success = fit_bum(pvalues_array)

# Compute scores with FDR=0.1
fdr = 0.1
node_scores = compute_node_scores(
    network_filtered,
    pvalues_filtered,
    lambda_param,
    alpha,
    fdr,
    missing_data_score=False
)

# Calculate prizes for reference
min_score = min(node_scores.values())
prizes = {node: score - min_score for node, score in node_scores.items()}
base_cost = np.median(list(prizes.values()))

print("="*120)
print("TESTING COMBINATIONS OF PRUNING METHODS AND num_clusters")
print("="*120)

# Test matrix
pruning_methods = ['gw', 'strong']
num_clusters_values = [1, 2, 3, 4, 5]

results = []

for pruning_method in pruning_methods:
    for num_clusters in num_clusters_values:
        module = detect_functional_module_core(
            network_filtered,
            node_scores,
            edge_weight_attr=None,
            c0=0.1,
            num_clusters=num_clusters,
            pruning=pruning_method
        )
        
        module_nodes = list(module.nodes())
        module_edges = module.number_of_edges()
        module_prize = sum(module.nodes[n].get('prize', 0) for n in module_nodes)
        module_cost = module_edges * base_cost
        module_benefit = module_prize - module_cost
        
        results.append({
            'pruning': pruning_method,
            'num_clusters': num_clusters,
            'nodes': len(module_nodes),
            'edges': module_edges,
            'prize': module_prize,
            'cost': module_cost,
            'benefit': module_benefit
        })
        
        print(f"Pruning={pruning_method:7s}, num_clusters={num_clusters} → Nodes={len(module_nodes):2d}, Edges={module_edges:2d}, Prize={module_prize:>8.3f}, Cost={module_cost:>8.3f}, Benefit={module_benefit:>8.3f}")

# Summary table
print(f"\n" + "="*120)
print("SUMMARY TABLE")
print("="*120)

print(f"\n{'Pruning':<10} {'Clusters':<12} {'Nodes':<8} {'Edges':<8} {'Prize':<12} {'Cost':<12} {'Benefit':<12}")
print("-" * 100)

for r in results:
    print(f"{r['pruning']:<10} {r['num_clusters']:<12} {r['nodes']:<8} {r['edges']:<8} {r['prize']:<12.6f} {r['cost']:<12.6f} {r['benefit']:<12.6f}")

# Find best configuration
best = max(results, key=lambda x: x['benefit'])
print(f"\n" + "="*120)
print(f"BEST CONFIGURATION: Pruning='{best['pruning']}', num_clusters={best['num_clusters']}")
print(f"  Nodes: {best['nodes']}, Benefit: {best['benefit']:.6f}")
print("="*120)

