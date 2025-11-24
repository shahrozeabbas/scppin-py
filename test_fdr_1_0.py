#!/usr/bin/env python
"""
Test with FDR=1.0, num_clusters=1, and default c0.

FDR=1.0 means essentially no filtering - all p-values should contribute.
This will show us what happens when ALL nodes are considered significant.
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

# FDR = 1.0 (essentially no filtering)
fdr = 1.0
print("="*80)
print(f"TEST: FDR=1.0, num_clusters=1, c0=0.1 (default)")
print("="*80)
print(f"\nFDR threshold: {fdr}")
print(f"BUM parameters: lambda={lambda_param:.6f}, alpha={alpha:.6f}")

# Compute node scores with FDR=1.0
node_scores = compute_node_scores(
    network_filtered,
    pvalues_filtered,
    lambda_param,
    alpha,
    fdr,
    missing_data_score=False
)

print(f"\nNode scores with FDR=1.0:")
print(f"  Total nodes with scores: {len(node_scores)}")

if node_scores:
    scores_array = np.array(list(node_scores.values()))
    print(f"  Min score: {np.min(scores_array):.6f}")
    print(f"  Max score: {np.max(scores_array):.6f}")
    print(f"  Mean score: {np.mean(scores_array):.6f}")
    print(f"  Median score: {np.median(scores_array):.6f}")
    print(f"  Std dev: {np.std(scores_array):.6f}")

# Run PCST with default c0=0.1
print(f"\nRunning PCST with num_clusters=1, c0=0.1 (default)...")
module = detect_functional_module_core(
    network_filtered,
    node_scores,
    edge_weight_attr=None,
    c0=0.1,  # default
    num_clusters=1
)

print(f"\nResults:")
print(f"  Total nodes in module: {len(module.nodes())}")
print(f"  Total edges in module: {module.number_of_edges()}")

if len(module.nodes()) > 0:
    module_prize = sum(module.nodes[n].get('prize', 0) for n in module.nodes())
    print(f"  Total prize: {module_prize:.6f}")
    
    # Calculate prizes for reference
    min_score = min(node_scores.values())
    prizes = {node: score - min_score for node, score in node_scores.items()}
    base_cost = np.median(list(prizes.values()))
    print(f"\nBase cost (median of prizes): {base_cost:.6f}")
    
    # Show all nodes
    nodes_list = list(module.nodes())
    sorted_nodes = sorted(nodes_list, key=lambda n: module.nodes[n].get('prize', 0), reverse=True)
    
    print(f"\nNodes in module:")
    for i, node in enumerate(sorted_nodes, 1):
        prize = module.nodes[node].get('prize', 0)
        score = module.nodes[node].get('score', 0)
        pval = pvalues_filtered.get(node, 'N/A')
        print(f"  {i}. {node}: prize={prize:.6f}, score={score:.6f}, p={pval}")
else:
    print("  Empty module returned!")

# Also show top 20 nodes overall for context
print(f"\n" + "="*80)
print("Top 20 nodes by prize (for context):")
print("="*80)
min_score = min(node_scores.values())
prizes = {node: score - min_score for node, score in node_scores.items()}
sorted_all_nodes = sorted(prizes.items(), key=lambda x: x[1], reverse=True)

for i, (node, prize) in enumerate(sorted_all_nodes[:20], 1):
    score = node_scores[node]
    pval = pvalues_filtered.get(node, 'N/A')
    print(f"  {i:2d}. {node:15s}: prize={prize:.6f}, score={score:.6f}, p={pval}")

