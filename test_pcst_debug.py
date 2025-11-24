#!/usr/bin/env python
"""
Debug PCST solver to understand why it returns small modules.
"""

import pandas as pd
import networkx as nx
import numpy as np
from scppin.core.scoring import compute_node_scores
from scppin.core.bum_model import fit_bum
import pcst_fast

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

# Calculate prizes and base_cost
min_score = min(node_scores.values())
prizes = {node: score - min_score for node, score in node_scores.items()}
base_cost = np.median(list(prizes.values()))

print("="*100)
print("PCST SOLVER DEBUG")
print("="*100)

print(f"\nBase cost: {base_cost:.6f}")
print(f"Network: {network_filtered.number_of_nodes()} nodes, {network_filtered.number_of_edges()} edges")
print(f"Nodes with prizes: {len(prizes)}")

# Examine prize distribution
prizes_array = np.array(list(prizes.values()))
print(f"\nPrize distribution:")
print(f"  Min: {np.min(prizes_array):.6f}")
print(f"  Max: {np.max(prizes_array):.6f}")
print(f"  Mean: {np.mean(prizes_array):.6f}")
print(f"  Median: {np.median(prizes_array):.6f}")
print(f"  > 0: {np.sum(prizes_array > 0)} nodes")
print(f"  > 1: {np.sum(prizes_array > 1)} nodes")
print(f"  > 10: {np.sum(prizes_array > 10)} nodes")

# Create index mapping
nodes = list(network_filtered.nodes())
node_to_idx = {node: idx for idx, node in enumerate(nodes)}
idx_to_node = {idx: node for node, idx in node_to_idx.items()}

# Prepare prizes array
prizes_vec = np.zeros(len(nodes))
for node, prize in prizes.items():
    if node in node_to_idx:
        prizes_vec[node_to_idx[node]] = prize

# Prepare edges and costs
num_edges = network_filtered.number_of_edges()
edges = np.zeros((num_edges, 2), dtype=int)
costs = np.zeros(num_edges, dtype=float)

for i, (u, v) in enumerate(network_filtered.edges()):
    edges[i] = [node_to_idx[u], node_to_idx[v]]
    costs[i] = base_cost

print(f"\nEdge costs:")
print(f"  All edges: {base_cost:.6f}")
print(f"  Total edges: {len(edges)}")

# Test with different parameters
print(f"\n" + "="*100)
print("TESTING DIFFERENT PCST PARAMETERS")
print("="*100)

test_configs = [
    {"num_clusters": 1, "pruning": "gw", "verbosity": 1},
    {"num_clusters": 1, "pruning": "strong", "verbosity": 1},
    {"num_clusters": 2, "pruning": "gw", "verbosity": 1},
    {"num_clusters": 5, "pruning": "gw", "verbosity": 1},
]

for config in test_configs:
    print(f"\nConfig: {config}")
    
    vertices, edges_in_solution = pcst_fast.pcst_fast(
        edges,
        prizes_vec,
        costs,
        -1,  # root_idx
        config["num_clusters"],
        config["pruning"],
        config["verbosity"]
    )
    
    solution_nodes = set()
    for edge_idx in edges_in_solution:
        if 0 <= edge_idx < len(edges):
            u, v = edges[edge_idx]
            solution_nodes.add(u)
            solution_nodes.add(v)
    
    for v_idx in vertices:
        if 0 <= v_idx < len(nodes):
            solution_nodes.add(v_idx)
    
    solution_node_names = [idx_to_node[idx] for idx in solution_nodes]
    solution_prize = sum(prizes[node] for node in solution_node_names)
    solution_edges = len(edges_in_solution)
    solution_cost = solution_edges * base_cost
    solution_benefit = solution_prize - solution_cost
    
    print(f"  Nodes selected: {len(solution_nodes)}")
    print(f"  Edges selected: {solution_edges}")
    print(f"  Total prize: {solution_prize:.6f}")
    print(f"  Total cost: {solution_cost:.6f}")
    print(f"  Net benefit: {solution_benefit:.6f}")
    print(f"  Top node: {sorted(solution_node_names, key=lambda n: prizes[n], reverse=True)[0] if solution_node_names else 'N/A'}")

print("\n" + "="*100)

