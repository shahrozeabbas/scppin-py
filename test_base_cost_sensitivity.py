#!/usr/bin/env python
"""
Test PCST with different base_cost values to find where it starts finding larger modules.
"""

import pandas as pd
import networkx as nx
import numpy as np
from scppin.core.scoring import compute_node_scores
from scppin.core.bum_model import fit_bum
from scppin.core.pcst_solver import solve_pcst, prepare_edge_costs

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

print("="*100)
print("BASE COST SENSITIVITY ANALYSIS")
print("="*100)

print(f"\nPrize statistics:")
print(f"  Min: {min(prizes.values()):.6f}")
print(f"  Median: {np.median(list(prizes.values())):.6f}")
print(f"  Mean: {np.mean(list(prizes.values())):.6f}")
print(f"  Q25: {np.percentile(list(prizes.values()), 25):.6f}")
print(f"  Q10: {np.percentile(list(prizes.values()), 10):.6f}")

# Test different base_cost values
base_costs_to_test = [
    ("Current (median)", np.median(list(prizes.values()))),
    ("Q10 (10th percentile)", np.percentile(list(prizes.values()), 10)),
    ("Q25 (25th percentile)", np.percentile(list(prizes.values()), 25)),
    ("0.1 * median", 0.1 * np.median(list(prizes.values()))),
    ("0.05 * median", 0.05 * np.median(list(prizes.values()))),
    ("0.01 * median", 0.01 * np.median(list(prizes.values()))),
]

print(f"\n" + "="*100)
print("TESTING DIFFERENT BASE COSTS")
print("="*100)

results = []

for name, base_cost in base_costs_to_test:
    print(f"\n{name}: base_cost = {base_cost:.6f}")
    
    # Prepare edge costs
    edge_costs = prepare_edge_costs(
        network_filtered,
        base_cost,
        edge_weight_attr=None,
        c0=0.1
    )
    
    # Solve PCST
    try:
        solution_nodes = solve_pcst(
            network_filtered,
            prizes,
            edge_costs,
            num_clusters=1
        )
        
        solution_prize = sum(prizes.get(n, 0) for n in solution_nodes)
        
        # Count edges properly: find minimum spanning tree
        if len(solution_nodes) > 1:
            subgraph = network_filtered.subgraph(solution_nodes)
            solution_edges = subgraph.number_of_edges()
        else:
            solution_edges = 0
        
        solution_cost = solution_edges * base_cost
        solution_benefit = solution_prize - solution_cost
        
        print(f"  Nodes: {len(solution_nodes)}")
        print(f"  Edges: {solution_edges}")
        print(f"  Total prize: {solution_prize:.6f}")
        print(f"  Total cost: {solution_cost:.6f}")
        print(f"  Net benefit: {solution_benefit:.6f}")
        
        # Top node
        if solution_nodes:
            top_node = max(solution_nodes, key=lambda n: prizes.get(n, 0))
            top_prize = prizes.get(top_node, 0)
            print(f"  Top node: {top_node} (prize={top_prize:.6f})")
        
        results.append({
            'name': name,
            'base_cost': base_cost,
            'num_nodes': len(solution_nodes),
            'num_edges': solution_edges,
            'total_prize': solution_prize,
            'total_cost': solution_cost,
            'net_benefit': solution_benefit,
            'solution_nodes': solution_nodes
        })
        
    except Exception as e:
        print(f"  ERROR: {str(e)}")

# Summary
print(f"\n" + "="*100)
print("SUMMARY: Impact of Base Cost on Module Size")
print("="*100)

print(f"\n{'Base Cost':<25} {'Nodes':<8} {'Edges':<8} {'Prize':<12} {'Cost':<12} {'Net Benefit':<12}")
print("-" * 100)

for r in results:
    print(f"{r['base_cost']:<25.6f} {r['num_nodes']:<8} {r['num_edges']:<8} {r['total_prize']:<12.6f} {r['total_cost']:<12.6f} {r['net_benefit']:<12.6f}")

# Find threshold where module size increases
print(f"\n" + "="*100)
print("FINDINGS")
print("="*100)

size_changes = []
for i in range(1, len(results)):
    if results[i]['num_nodes'] != results[i-1]['num_nodes']:
        size_changes.append({
            'from': results[i-1]['num_nodes'],
            'to': results[i]['num_nodes'],
            'at_base_cost': results[i]['base_cost'],
            'between': f"{results[i-1]['base_cost']:.6f} and {results[i]['base_cost']:.6f}"
        })

if size_changes:
    print(f"\nModule size changes:")
    for change in size_changes:
        print(f"  {change['from']} → {change['to']} nodes at base_cost {change['between']}")
else:
    print(f"\nModule size remains constant across all tested base_cost values")
    print(f"Currently at: {results[0]['num_nodes']} nodes")

print("\n" + "="*100)

