#!/usr/bin/env python
"""
Test if PCST behaves differently when all nodes have prizes vs only filtered nodes.

Maybe the issue is that we're only including nodes with FDR-filtered scores?
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

print("="*100)
print("TEST: FDR-FILTERED vs FDR=1.0 PRIZES")
print("="*100)

print(f"\nScoring results (FDR={fdr}):")
print(f"  Nodes with scores: {len(node_scores)}")
print(f"  Nodes without scores: {len(network_filtered.nodes()) - len(node_scores)}")

# Calculate prizes (FDR=0.1 only)
min_score = min(node_scores.values())
prizes_fdr = {node: score - min_score for node, score in node_scores.items()}
base_cost_fdr = np.median(list(prizes_fdr.values()))

print(f"\nPrizes (FDR=0.1 filtered):")
print(f"  Nodes with prizes: {len(prizes_fdr)}")
print(f"  Min prize: {min(prizes_fdr.values()):.6f}")
print(f"  Max prize: {max(prizes_fdr.values()):.6f}")
print(f"  Base cost (median): {base_cost_fdr:.6f}")

# Now compute scores with FDR=1.0
node_scores_fdr1 = compute_node_scores(
    network_filtered,
    pvalues_filtered,
    lambda_param,
    alpha,
    1.0,  # FDR=1.0
    missing_data_score=False
)

print(f"\nScoring results (FDR=1.0):")
print(f"  Nodes with scores: {len(node_scores_fdr1)}")

# Calculate prizes (FDR=1.0 - all nodes)
min_score_fdr1 = min(node_scores_fdr1.values())
prizes_fdr1 = {node: score - min_score_fdr1 for node, score in node_scores_fdr1.items()}
base_cost_fdr1 = np.median(list(prizes_fdr1.values()))

print(f"\nPrizes (FDR=1.0):")
print(f"  Nodes with prizes: {len(prizes_fdr1)}")
print(f"  Min prize: {min(prizes_fdr1.values()):.6f}")
print(f"  Max prize: {max(prizes_fdr1.values()):.6f}")
print(f"  Base cost (median): {base_cost_fdr1:.6f}")

# Test PCST with both
print(f"\n" + "="*100)
print("PCST RESULTS")
print("="*100)

# Test 1: FDR=0.1 (current)
print(f"\nTest 1: FDR=0.1 filtered prizes")
edge_costs_fdr = prepare_edge_costs(network_filtered, base_cost_fdr, edge_weight_attr=None, c0=0.1)
solution_nodes_fdr = solve_pcst(network_filtered, prizes_fdr, edge_costs_fdr, num_clusters=1)

solution_prize_fdr = sum(prizes_fdr.get(n, 0) for n in solution_nodes_fdr)
subgraph_fdr = network_filtered.subgraph(solution_nodes_fdr)
solution_edges_fdr = subgraph_fdr.number_of_edges()
solution_cost_fdr = solution_edges_fdr * base_cost_fdr
solution_benefit_fdr = solution_prize_fdr - solution_cost_fdr

print(f"  Nodes: {len(solution_nodes_fdr)}")
print(f"  Edges: {solution_edges_fdr}")
print(f"  Prize: {solution_prize_fdr:.6f}")
print(f"  Cost: {solution_cost_fdr:.6f}")
print(f"  Net benefit: {solution_benefit_fdr:.6f}")
print(f"  Top node: {max(solution_nodes_fdr, key=lambda n: prizes_fdr.get(n, 0))}")

# Test 2: FDR=1.0 (all nodes)
print(f"\nTest 2: FDR=1.0 all prizes")
edge_costs_fdr1 = prepare_edge_costs(network_filtered, base_cost_fdr1, edge_weight_attr=None, c0=0.1)
solution_nodes_fdr1 = solve_pcst(network_filtered, prizes_fdr1, edge_costs_fdr1, num_clusters=1)

solution_prize_fdr1 = sum(prizes_fdr1.get(n, 0) for n in solution_nodes_fdr1)
subgraph_fdr1 = network_filtered.subgraph(solution_nodes_fdr1)
solution_edges_fdr1 = subgraph_fdr1.number_of_edges()
solution_cost_fdr1 = solution_edges_fdr1 * base_cost_fdr1
solution_benefit_fdr1 = solution_prize_fdr1 - solution_cost_fdr1

print(f"  Nodes: {len(solution_nodes_fdr1)}")
print(f"  Edges: {solution_edges_fdr1}")
print(f"  Prize: {solution_prize_fdr1:.6f}")
print(f"  Cost: {solution_cost_fdr1:.6f}")
print(f"  Net benefit: {solution_benefit_fdr1:.6f}")
print(f"  Top node: {max(solution_nodes_fdr1, key=lambda n: prizes_fdr1.get(n, 0))}")

# Test 3: Force large solution by using very small base cost
print(f"\nTest 3: FDR=0.1 with base_cost=0.001 (aggressive)")
edge_costs_small = prepare_edge_costs(network_filtered, 0.001, edge_weight_attr=None, c0=0.1)
solution_nodes_small = solve_pcst(network_filtered, prizes_fdr, edge_costs_small, num_clusters=1)

solution_prize_small = sum(prizes_fdr.get(n, 0) for n in solution_nodes_small)
subgraph_small = network_filtered.subgraph(solution_nodes_small)
solution_edges_small = subgraph_small.number_of_edges()
solution_cost_small = solution_edges_small * 0.001
solution_benefit_small = solution_prize_small - solution_cost_small

print(f"  Nodes: {len(solution_nodes_small)}")
print(f"  Edges: {solution_edges_small}")
print(f"  Prize: {solution_prize_small:.6f}")
print(f"  Cost: {solution_cost_small:.6f}")
print(f"  Net benefit: {solution_benefit_small:.6f}")
if solution_nodes_small:
    print(f"  Top node: {max(solution_nodes_small, key=lambda n: prizes_fdr.get(n, 0))}")

print("\n" + "="*100)

