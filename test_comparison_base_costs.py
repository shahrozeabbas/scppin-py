#!/usr/bin/env python
"""
Direct comparison: median base_cost vs min positive prize base_cost
"""

import pandas as pd
import networkx as nx
import numpy as np
from scppin.core.scoring import compute_node_scores
from scppin.core.pcst_solver import solve_pcst, prepare_edge_costs
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

# Calculate prizes
min_score = min(node_scores.values())
prizes = {node: score - min_score for node, score in node_scores.items()}
prizes_list = list(prizes.values())

print("="*120)
print("DIRECT COMPARISON: Different Base Cost Strategies")
print("="*120)

# Test different base_cost strategies
strategies = [
    ("Median", np.median(prizes_list)),
    ("Min positive prize", np.min([p for p in prizes_list if p > 0])),
    ("Q25 (25th percentile)", np.percentile(prizes_list, 25)),
    ("Q10 (10th percentile)", np.percentile(prizes_list, 10)),
    ("Very small (0.0001)", 0.0001),
]

results = []

for name, base_cost in strategies:
    print(f"\n{name} (base_cost={base_cost:.6f}):")
    
    edge_costs = prepare_edge_costs(network_filtered, base_cost, edge_weight_attr=None, c0=0.1)
    solution_nodes = solve_pcst(network_filtered, prizes, edge_costs, num_clusters=1)
    
    solution_prize = sum(prizes.get(n, 0) for n in solution_nodes)
    subgraph = network_filtered.subgraph(solution_nodes)
    solution_edges = subgraph.number_of_edges()
    solution_cost = solution_edges * base_cost
    solution_benefit = solution_prize - solution_cost
    
    print(f"  Nodes: {len(solution_nodes)}")
    print(f"  Edges: {solution_edges}")
    print(f"  Prize: {solution_prize:.6f}")
    print(f"  Cost: {solution_cost:.6f}")
    print(f"  Benefit: {solution_benefit:.6f}")
    
    results.append({
        'name': name,
        'base_cost': base_cost,
        'nodes': len(solution_nodes),
        'edges': solution_edges,
        'prize': solution_prize,
        'cost': solution_cost,
        'benefit': solution_benefit
    })

print(f"\n" + "="*120)
print("SUMMARY TABLE")
print("="*120)

print(f"\n{'Strategy':<25} {'Base Cost':<15} {'Nodes':<8} {'Edges':<8} {'Prize':<12} {'Benefit':<12}")
print("-" * 100)

for r in results:
    print(f"{r['name']:<25} {r['base_cost']:<15.6f} {r['nodes']:<8} {r['edges']:<8} {r['prize']:<12.6f} {r['benefit']:<12.6f}")

# Check if any strategy changed node count
unique_node_counts = set(r['nodes'] for r in results)
print(f"\n" + "="*120)
print(f"KEY FINDING:")
print(f"  Unique node counts across all strategies: {unique_node_counts}")
if len(unique_node_counts) == 1:
    print(f"  ❌ ALL strategies return {results[0]['nodes']} nodes")
    print(f"  ❌ Changing base_cost has NO EFFECT on module size")
else:
    print(f"  ✓ Different strategies produce different node counts")

print("\n" + "="*120)

