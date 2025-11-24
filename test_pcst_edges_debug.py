#!/usr/bin/env python
"""
Deep dive into pcst_fast return values.
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
print("PCST RETURN VALUE ANALYSIS")
print("="*100)

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

print(f"\nNetwork setup:")
print(f"  Nodes: {len(nodes)}")
print(f"  Edges: {len(edges)}")
print(f"  All edges have cost: {base_cost:.6f}")

# Call pcst_fast
vertices, edges_in_solution = pcst_fast.pcst_fast(
    edges,
    prizes_vec,
    costs,
    -1,  # root_idx
    1,   # num_clusters
    "gw",
    0    # verbosity
)

print(f"\npcst_fast return values:")
print(f"  vertices: type={type(vertices)}, len={len(vertices)}")
print(f"  edges_in_solution: type={type(edges_in_solution)}, len={len(edges_in_solution)}")
print(f"  vertices[:10] = {vertices[:10]}")
print(f"  edges_in_solution[:10] = {edges_in_solution[:10]}")

# Reconstruct solution correctly
solution_indices = set()

# Add from edges_in_solution
for edge_idx in edges_in_solution:
    if 0 <= edge_idx < len(edges):
        u, v = edges[edge_idx]
        solution_indices.add(u)
        solution_indices.add(v)

# Add from vertices
for v_idx in vertices:
    if 0 <= v_idx < len(nodes):
        solution_indices.add(v_idx)

print(f"\nReconstructed solution:")
print(f"  Total node indices from edges: {len(set(u for edge_idx in edges_in_solution if 0 <= edge_idx < len(edges) for u, v in [edges[edge_idx]]) | set(v for edge_idx in edges_in_solution if 0 <= edge_idx < len(edges) for u, v in [edges[edge_idx]]))}")
print(f"  Total node indices from vertices: {len([v for v in vertices if 0 <= v < len(nodes)])}")
print(f"  Total unique nodes: {len(solution_indices)}")
print(f"  Unique nodes: {[idx_to_node[idx] for idx in sorted(list(solution_indices))]}")

# Count actual edges in solution
actual_edges_in_solution = set()
for edge_idx in edges_in_solution:
    if 0 <= edge_idx < len(edges):
        actual_edges_in_solution.add(edge_idx)

print(f"\nEdge statistics:")
print(f"  edges_in_solution indices: {len(actual_edges_in_solution)} (should be # edges connecting the nodes)")
print(f"  With {len(solution_indices)} nodes, max possible edges: {len(solution_indices) - 1} (tree)")

# Get solution node names
solution_node_names = [idx_to_node[idx] for idx in solution_indices]
solution_prize = sum(prizes[node] for node in solution_node_names)
solution_cost = len(actual_edges_in_solution) * base_cost
solution_benefit = solution_prize - solution_cost

print(f"\nSolution metrics:")
print(f"  Nodes: {solution_node_names}")
print(f"  Prize: {solution_prize:.6f}")
print(f"  Edges: {len(actual_edges_in_solution)}")
print(f"  Cost (edges * base_cost): {solution_cost:.6f}")
print(f"  Net benefit: {solution_benefit:.6f}")

print("\n" + "="*100)

