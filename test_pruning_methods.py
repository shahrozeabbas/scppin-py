#!/usr/bin/env python
"""
Test different PCST pruning methods to see if any find larger/better modules.

Available pruning methods:
- 'gw': Goemans-Williamson (faster, approximate)
- 'strong': Strong pruning (potentially slower, possibly better quality)
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

print("="*100)
print("TESTING DIFFERENT PCST PRUNING METHODS")
print("="*100)

print(f"\nSetup:")
print(f"  FDR: {fdr}")
print(f"  Network: {network_filtered.number_of_nodes()} nodes, {network_filtered.number_of_edges()} edges")
print(f"  Nodes with scores: {len(node_scores)}")

# Test different pruning methods
pruning_methods = ['gw', 'strong']
results = []

for pruning_method in pruning_methods:
    print(f"\n" + "="*100)
    print(f"PRUNING METHOD: '{pruning_method}'")
    print("="*100)
    
    module = detect_functional_module_core(
        network_filtered,
        node_scores,
        edge_weight_attr=None,
        c0=0.1,
        num_clusters=1,
        pruning=pruning_method
    )
    
    module_nodes = list(module.nodes())
    module_edges = module.number_of_edges()
    module_prize = sum(module.nodes[n].get('prize', 0) for n in module_nodes)
    
    print(f"\nResults:")
    print(f"  Nodes selected: {len(module_nodes)}")
    print(f"  Edges in module: {module_edges}")
    print(f"  Total prize: {module_prize:.6f}")
    
    # Base cost
    base_cost = np.median(list(prizes.values()))
    module_cost = module_edges * base_cost
    module_benefit = module_prize - module_cost
    
    print(f"  Base cost: {base_cost:.6f}")
    print(f"  Total cost: {module_cost:.6f}")
    print(f"  Net benefit: {module_benefit:.6f}")
    
    # Show nodes
    if module_nodes:
        sorted_nodes = sorted(module_nodes, key=lambda n: module.nodes[n].get('prize', 0), reverse=True)
        print(f"\n  Top 5 nodes by prize:")
        for i, node in enumerate(sorted_nodes[:5], 1):
            prize = module.nodes[node].get('prize', 0)
            score = module.nodes[node].get('score', 0)
            pval = pvalues_filtered.get(node, 'N/A')
            print(f"    {i}. {node:15s}: prize={prize:>10.6f}, score={score:>10.6f}, p={pval}")
    
    results.append({
        'pruning': pruning_method,
        'nodes': len(module_nodes),
        'edges': module_edges,
        'prize': module_prize,
        'cost': module_cost,
        'benefit': module_benefit,
        'node_list': module_nodes
    })

# Summary comparison
print(f"\n" + "="*100)
print("SUMMARY COMPARISON")
print("="*100)

print(f"\n{'Pruning':<10} {'Nodes':<8} {'Edges':<8} {'Prize':<12} {'Cost':<12} {'Benefit':<12}")
print("-" * 70)

for r in results:
    print(f"{r['pruning']:<10} {r['nodes']:<8} {r['edges']:<8} {r['prize']:<12.6f} {r['cost']:<12.6f} {r['benefit']:<12.6f}")

# Find best
best = max(results, key=lambda x: x['benefit'])
worst = min(results, key=lambda x: x['benefit'])

print(f"\n" + "="*100)
print("FINDINGS")
print("="*100)

if best['benefit'] > worst['benefit']:
    diff = best['benefit'] - worst['benefit']
    pct_better = (diff / abs(worst['benefit'])) * 100 if worst['benefit'] != 0 else float('inf')
    print(f"\nBest pruning method: '{best['pruning']}'")
    print(f"  Benefit: {best['benefit']:.6f}")
    print(f"  vs '{worst['pruning']}': {diff:+.6f} ({pct_better:+.1f}%)")
else:
    print(f"\nAll pruning methods produced the same or similar results")
    print(f"  Benefit range: {worst['benefit']:.6f} to {best['benefit']:.6f}")

print("\n" + "="*100)

