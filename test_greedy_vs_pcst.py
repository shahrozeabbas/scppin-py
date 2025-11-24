#!/usr/bin/env python
"""
Compare what a greedy algorithm would do vs what PCST does.

This helps diagnose if PCST is truly optimal or if it's stuck in local optima.
"""

import pandas as pd
import networkx as nx
import numpy as np
from scppin.core.scoring import compute_node_scores
from scppin.core.pcst_solver import detect_functional_module_core, prepare_edge_costs
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

# Calculate prizes and base_cost
min_score = min(node_scores.values())
prizes = {node: score - min_score for node, score in node_scores.items()}
base_cost = np.median(list(prizes.values()))

print("="*100)
print("GREEDY vs PCST COMPARISON")
print("="*100)

print(f"\nBase cost (median prize): {base_cost:.6f}")
print(f"Total nodes with scores: {len(node_scores)}")

# Get PCST solution
module = detect_functional_module_core(
    network_filtered,
    node_scores,
    edge_weight_attr=None,
    c0=0.1,
    num_clusters=1
)

pcst_nodes = list(module.nodes())
pcst_prize = sum(prizes.get(n, 0) for n in pcst_nodes)
pcst_edges = module.number_of_edges()
pcst_cost = pcst_edges * base_cost

print(f"\nPCST SOLUTION:")
print(f"  Nodes: {pcst_nodes}")
print(f"  Total prize: {pcst_prize:.6f}")
print(f"  Total edges: {pcst_edges}")
print(f"  Total cost: {pcst_cost:.6f}")
print(f"  Net benefit: {pcst_prize - pcst_cost:.6f}")

# Greedy algorithm: pick highest prize node, then greedily expand
print(f"\n" + "="*100)
print("GREEDY ALGORITHM: Start with highest prize, greedily add best profit nodes")
print("="*100)

sorted_nodes = sorted(prizes.items(), key=lambda x: x[1], reverse=True)

# Start with top node
greedy_solution = [sorted_nodes[0][0]]
greedy_prize = sorted_nodes[0][1]
greedy_edges = 0
greedy_cost = 0

print(f"\nStart with top node: {greedy_solution[0]} (prize={sorted_nodes[0][1]:.6f})")

# Greedily add profitable neighbors
added = True
iteration = 1
while added and iteration <= 20:
    added = False
    
    # Find all neighbors of current solution
    neighbors = set()
    for node in greedy_solution:
        neighbors.update(network_filtered.neighbors(node))
    
    # Remove already included nodes
    neighbors = neighbors - set(greedy_solution)
    
    if not neighbors:
        break
    
    # Find best neighbor to add
    best_neighbor = None
    best_profit = 0
    best_path_cost = 0
    
    for neighbor in neighbors:
        # Cost to add this neighbor = shortest path distance * base_cost
        min_dist = float('inf')
        for node in greedy_solution:
            try:
                dist = nx.shortest_path_length(network_filtered, node, neighbor)
                min_dist = min(min_dist, dist)
            except nx.NetworkXNoPath:
                continue
        
        if min_dist == float('inf'):
            continue
        
        path_cost = min_dist * base_cost
        net_profit = prizes[neighbor] - path_cost
        
        if net_profit > best_profit:
            best_profit = net_profit
            best_neighbor = neighbor
            best_path_cost = path_cost
    
    if best_neighbor and best_profit > 0:
        greedy_solution.append(best_neighbor)
        greedy_prize += prizes[best_neighbor]
        greedy_cost += best_path_cost
        print(f"  Iteration {iteration}: Added {best_neighbor} (prize={prizes[best_neighbor]:.6f}, path_cost={best_path_cost:.6f}, net_profit={best_profit:.6f})")
        added = True
        iteration += 1
    else:
        break

# Count edges in greedy solution
greedy_subgraph = network_filtered.subgraph(greedy_solution)
greedy_edges = greedy_subgraph.number_of_edges()

print(f"\nGREEDY SOLUTION:")
print(f"  Nodes: {greedy_solution}")
print(f"  Total prize: {greedy_prize:.6f}")
print(f"  Total edges in subgraph: {greedy_edges}")
print(f"  Total cost (estimated): {greedy_cost:.6f}")
print(f"  Net benefit: {greedy_prize - greedy_cost:.6f}")

# Compare
print(f"\n" + "="*100)
print("COMPARISON")
print("="*100)

print(f"\nPCST NET BENEFIT:   {pcst_prize - pcst_cost:.6f}")
print(f"GREEDY NET BENEFIT: {greedy_prize - greedy_cost:.6f}")

if greedy_prize - greedy_cost > pcst_prize - pcst_cost:
    diff = (greedy_prize - greedy_cost) - (pcst_prize - pcst_cost)
    pct_better = (diff / (pcst_prize - pcst_cost)) * 100 if pcst_prize - pcst_cost > 0 else float('inf')
    print(f"\n⚠️  GREEDY is {diff:.6f} better ({pct_better:.1f}% improvement!)")
    print(f"   This suggests PCST is not finding the optimal solution.")
else:
    print(f"\n✓ PCST is better (or equal)")

print("\n" + "="*100)

