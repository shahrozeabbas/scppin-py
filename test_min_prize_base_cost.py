#!/usr/bin/env python
"""
Test PCST with minimum prize as base_cost instead of median.

This makes edges very cheap, which should encourage PCST to explore larger modules.
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

print("="*100)
print("TEST: MINIMUM PRIZE AS BASE_COST")
print("="*100)

# Calculate prizes for reference
min_score = min(node_scores.values())
prizes = {node: score - min_score for node, score in node_scores.items()}
prizes_list = list(prizes.values())

print(f"\nPrize statistics:")
print(f"  Min prize: {np.min(prizes_list):.6f}")
print(f"  Median prize: {np.median(prizes_list):.6f}")
print(f"  Mean prize: {np.mean(prizes_list):.6f}")

positive_prizes = [p for p in prizes_list if p > 0]
if positive_prizes:
    min_positive = np.min(positive_prizes)
    print(f"  Min positive prize: {min_positive:.6f}")
else:
    print(f"  No positive prizes")

print(f"\nBase cost comparison:")
print(f"  Median-based (old): {np.median(prizes_list):.6f}")
if positive_prizes:
    min_positive = np.min(positive_prizes)
    print(f"  Min-based (new): {min_positive:.6f}")
    reduction = (1 - min_positive / np.median(prizes_list)) * 100
    print(f"  Reduction: {reduction:.1f}%")
else:
    print(f"  Min-based (new): N/A")

# Test with different num_clusters
print(f"\n" + "="*100)
print("RESULTS WITH MINIMUM PRIZE BASE_COST")
print("="*100)

for num_clusters in [1, 2, 3, 4, 5]:
    print(f"\nnum_clusters = {num_clusters}:")
    
    module = detect_functional_module_core(
        network_filtered,
        node_scores,
        edge_weight_attr=None,
        c0=0.1,
        num_clusters=num_clusters,
        pruning='gw'
    )
    
    module_nodes = list(module.nodes())
    module_edges = module.number_of_edges()
    module_prize = sum(module.nodes[n].get('prize', 0) for n in module_nodes)
    
    # Calculate cost with the new base_cost
    if positive_prizes:
        base_cost = np.min(positive_prizes)
        module_cost = module_edges * base_cost
        module_benefit = module_prize - module_cost
    else:
        base_cost = np.median(prizes_list)
        module_cost = module_edges * base_cost
        module_benefit = module_prize - module_cost
    
    print(f"  Nodes: {len(module_nodes)}")
    print(f"  Edges: {module_edges}")
    print(f"  Prize: {module_prize:.6f}")
    print(f"  Cost (base={base_cost:.6f}): {module_cost:.6f}")
    print(f"  Benefit: {module_benefit:.6f}")
    
    if module_nodes:
        top_nodes = sorted(module_nodes, key=lambda n: module.nodes[n].get('prize', 0), reverse=True)[:3]
        top_nodes_str = ', '.join([f'{n} ({module.nodes[n].get("prize", 0):.3f})' for n in top_nodes])
        print(f"  Top nodes: {top_nodes_str}")

print("\n" + "="*100)
print("COMPARISON WITH OLD MEDIAN-BASED APPROACH")
print("="*100)

print("\nOld (median base_cost=0.729):")
print("  num_clusters=3: nodes=3, prize=4.377, benefit=3.648")

min_bc = min(positive_prizes) if positive_prizes else 0
print(f"\nNew (min prize base_cost={min_bc:.6f}):")
print("  num_clusters=1-5: tested above")

print("\n" + "="*100)

