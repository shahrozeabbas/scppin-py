#!/usr/bin/env python
"""
Analyze the cost of connecting high-prize nodes to the current solution.

This helps understand why PCST isn't including top-prize nodes even when
they're theoretically available in the network.
"""

import pandas as pd
import networkx as nx
import numpy as np
from scppin.core.scoring import compute_node_scores
from scppin.core.bum_model import fit_bum
from collections import defaultdict

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

# Calculate prizes
min_score = min(node_scores.values())
prizes = {node: score - min_score for node, score in node_scores.items()}
base_cost = np.median(list(prizes.values()))

print("="*100)
print("COST ANALYSIS: Adding High-Prize Nodes to Current Solution")
print("="*100)

print(f"\nBase cost (median prize): {base_cost:.6f}")

# Current solution from FDR=0.1, num_clusters=1
# Need to actually compute it
from scppin.core.pcst_solver import detect_functional_module_core

module = detect_functional_module_core(
    network_filtered,
    node_scores,
    edge_weight_attr=None,
    c0=0.1,
    num_clusters=1
)

current_solution = list(module.nodes())

print(f"\nCurrent solution (num_clusters=1, FDR=0.1): {current_solution}")
print(f"Current total prize: {sum(prizes.get(n, 0) for n in current_solution):.6f}")
print(f"Current total edges: {module.number_of_edges()}")

# Get top 20 nodes by prize
top_nodes_by_prize = sorted(prizes.items(), key=lambda x: x[1], reverse=True)[:20]

print("\n" + "="*100)
print("COST OF CONNECTING EACH TOP-20 NODE TO CURRENT SOLUTION")
print("="*100)

# For each top node, find shortest path to current solution
results = []

for rank, (node, prize) in enumerate(top_nodes_by_prize, 1):
    if node in current_solution:
        results.append({
            'rank': rank,
            'node': node,
            'prize': prize,
            'in_solution': True,
            'min_distance': 0,
            'cost_path_edges': 0,
            'total_cost': 0,
            'net_benefit': prize
        })
        continue
    
    # Find shortest path to any node in current solution
    min_distance = float('inf')
    closest_node = None
    
    for solution_node in current_solution:
        try:
            distance = nx.shortest_path_length(network_filtered, node, solution_node)
            if distance < min_distance:
                min_distance = distance
                closest_node = solution_node
        except nx.NetworkXNoPath:
            continue
    
    if min_distance == float('inf'):
        # Not connected to current solution
        results.append({
            'rank': rank,
            'node': node,
            'prize': prize,
            'in_solution': False,
            'min_distance': float('inf'),
            'cost_path_edges': float('inf'),
            'total_cost': float('inf'),
            'net_benefit': float('-inf')
        })
    else:
        # Cost = number of edges on path * base_cost
        # (assuming all edges have cost = base_cost)
        cost_path_edges = min_distance
        total_cost = cost_path_edges * base_cost
        net_benefit = prize - total_cost
        
        results.append({
            'rank': rank,
            'node': node,
            'prize': prize,
            'in_solution': False,
            'min_distance': min_distance,
            'cost_path_edges': cost_path_edges,
            'total_cost': total_cost,
            'net_benefit': net_benefit,
            'closest_to': closest_node
        })

# Display results
print(f"\n{'Rank':<5} {'Node':<15} {'Prize':<12} {'Distance':<10} {'Cost':<12} {'Net Benefit':<12} {'Connected To':<15}")
print("-" * 100)

for r in results:
    rank = r['rank']
    node = r['node']
    prize = r['prize']
    
    if r['in_solution']:
        print(f"{rank:<5} {node:<15} {prize:>11.6f} {'IN SOL':<10} {0:>11.6f} {prize:>11.6f} {'':<15}")
    elif r['min_distance'] == float('inf'):
        print(f"{rank:<5} {node:<15} {prize:>11.6f} {'NOT CONN':<10} {'∞':>11} {'∞':>11} {'':<15}")
    else:
        net = r['net_benefit']
        status = "✓ PROFITABLE" if net > 0 else "✗ LOSS"
        cost_str = f"{r['cost_path_edges']:.0f} edges"
        print(f"{rank:<5} {node:<15} {prize:>11.6f} {cost_str:<10} {r['total_cost']:>11.6f} {net:>11.6f} {r['closest_to']:<15} {status}")

# Summary statistics
print("\n" + "="*100)
print("SUMMARY")
print("="*100)

in_solution = [r for r in results if r.get('in_solution', False)]
profitable = [r for r in results if not r.get('in_solution', False) and r.get('net_benefit', float('-inf')) > 0 and r.get('net_benefit', float('-inf')) != float('-inf')]
not_profitable = [r for r in results if not r.get('in_solution', False) and r.get('net_benefit', 0) <= 0 and r.get('net_benefit', float('-inf')) != float('-inf')]
not_connected = [r for r in results if r.get('min_distance', float('inf')) == float('inf')]

print(f"\nCurrent solution nodes: {current_solution}")
print(f"In top-20 results: {len(in_solution)} nodes")
print(f"Profitable to add (from top-20): {len(profitable)} nodes")
if profitable:
    best_profitable = max(profitable, key=lambda x: x['net_benefit'])
    print(f"  Best: {best_profitable['node']} (prize={best_profitable['prize']:.6f}, net benefit={best_profitable['net_benefit']:.6f})")

print(f"Not profitable to add: {len(not_profitable)} nodes")
if not_profitable:
    least_bad = max(not_profitable, key=lambda x: x['net_benefit'])
    print(f"  Best of unprofitable: {least_bad['node']} (prize={least_bad['prize']:.6f}, net benefit={least_bad['net_benefit']:.6f})")

print(f"Not connected to solution: {len(not_connected)} nodes")
if not_connected:
    top_not_connected = sorted(not_connected, key=lambda x: x['prize'], reverse=True)[0]
    print(f"  Top prize not connected: {top_not_connected['node']} (prize={top_not_connected['prize']:.6f})")

# Additional insight: what if we lowered base_cost?
print("\n" + "="*100)
print("SENSITIVITY ANALYSIS: Impact of Base Cost on Profitability")
print("="*100)

print(f"\nCurrent base_cost: {base_cost:.6f}")
print(f"\nFor each top-5 non-solution nodes, at what base_cost become profitable?")
print()

top_non_solution = [r for r in results if not r['in_solution'] and r['min_distance'] != float('inf')][:5]

for r in top_non_solution:
    # profit = prize - (num_edges * base_cost)
    # profit = 0 when: base_cost = prize / num_edges
    breakeven_cost = r['prize'] / r['cost_path_edges'] if r['cost_path_edges'] > 0 else 0
    
    if breakeven_cost < base_cost:
        reduction_needed = base_cost - breakeven_cost
        pct_reduction = (reduction_needed / base_cost) * 100
        print(f"{r['node']:15s}: Need base_cost ≤ {breakeven_cost:.6f} (reduce by {reduction_needed:.6f}, {pct_reduction:.1f}%)")
    else:
        print(f"{r['node']:15s}: Already profitable at current base_cost")

print("\n" + "="*100)

