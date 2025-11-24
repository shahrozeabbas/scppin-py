#!/usr/bin/env python
"""
Analyze prize and cost distributions before PCST solving.
Shows where top prizes fall relative to edge costs.
"""

import pandas as pd
import networkx as nx
import numpy as np
from scppin.core.scoring import compute_node_scores
from scppin.core.pcst_solver import prepare_edge_costs
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

def main():
    print("="*80)
    print("PRIZE AND COST DISTRIBUTION ANALYSIS")
    print("="*80)
    
    # Load data
    pvalues = load_pvalues('tests/magma_gene_symbol_results.tsv')
    network = load_network_with_weights('tests/denoised_mg_ad_fava_network.tsv')
    
    genes_in_pvalues = set(pvalues.keys())
    genes_in_network = set(network.nodes())
    overlap = genes_in_pvalues & genes_in_network
    network_filtered = network.subgraph(overlap).copy()
    pvalues_filtered = {k: pvalues[k] for k in overlap}
    
    pvalues_array = np.array(list(pvalues_filtered.values()))
    lambda_param, alpha, success = fit_bum(pvalues_array)
    
    fdr = 0.01
    node_scores = compute_node_scores(
        network_filtered,
        pvalues_filtered,
        lambda_param,
        alpha,
        fdr,
        missing_data_score=False
    )
    
    min_score = min(node_scores.values())
    prizes = {node: score - min_score for node, score in node_scores.items()}
    base_cost = np.mean(list(prizes.values()))
    
    print(f"\n[1] BASIC PARAMETERS")
    print("-"*80)
    print(f"Network nodes: {network_filtered.number_of_nodes()}")
    print(f"Network edges: {network_filtered.number_of_edges()}")
    print(f"Nodes with scores: {len(node_scores)}")
    print(f"Min score: {min_score:.6f}")
    print(f"Base cost (mean prize): {base_cost:.6f}")
    
    # Prize distribution
    print(f"\n[2] PRIZE DISTRIBUTION")
    print("-"*80)
    
    prizes_array = np.array(list(prizes.values()))
    print(f"Prize statistics:")
    print(f"  Min: {np.min(prizes_array):.6f}")
    print(f"  Q1: {np.percentile(prizes_array, 25):.6f}")
    print(f"  Median: {np.median(prizes_array):.6f}")
    print(f"  Q3: {np.percentile(prizes_array, 75):.6f}")
    print(f"  Mean: {np.mean(prizes_array):.6f}")
    print(f"  Max: {np.max(prizes_array):.6f}")
    print(f"  Std: {np.std(prizes_array):.6f}")
    
    print(f"\nPrize quantiles:")
    for q in [50, 75, 90, 95, 99]:
        val = np.percentile(prizes_array, q)
        count = np.sum(prizes_array >= val)
        print(f"  Top {100-q}% (>= {val:.6f}): {count} nodes")
    
    # Top prizes
    print(f"\nTop 20 prizes:")
    sorted_by_prize = sorted(prizes.items(), key=lambda x: x[1], reverse=True)
    for i, (node, prize) in enumerate(sorted_by_prize[:20], 1):
        score = node_scores[node]
        pval = pvalues_filtered[node]
        ratio = prize / base_cost if base_cost > 0 else 0
        print(f"  {i:2d}. {node:15s}: prize={prize:8.4f} ({ratio:5.2f}x cost), score={score:8.4f}, p={pval:.6f}")
    
    # Edge cost distribution (uniform)
    print(f"\n[3] UNIFORM EDGE COSTS (no weights)")
    print("-"*80)
    
    edge_costs_uniform = prepare_edge_costs(
        network_filtered,
        base_cost,
        edge_weight_attr=None,
        c0=0.01
    )
    
    costs_uniform = np.array(list(edge_costs_uniform.values()))
    print(f"All uniform costs = base_cost = {base_cost:.6f}")
    print(f"Total edges: {len(costs_uniform)}")
    
    # Compare top prizes to cost
    print(f"\nCost-to-prize ratios for top nodes:")
    for i, (node, prize) in enumerate(sorted_by_prize[:10], 1):
        ratio = base_cost / prize if prize > 0 else float('inf')
        print(f"  Top {i}: {prize:.4f} prize, need {ratio:.2f}x this node's prize to break even on 1 edge")
    
    # Edge cost distribution (weighted)
    print(f"\n[4] WEIGHTED EDGE COSTS")
    print("-"*80)
    
    edge_costs_weighted = prepare_edge_costs(
        network_filtered,
        base_cost,
        edge_weight_attr='weight',
        c0=0.01
    )
    
    costs_weighted = np.array(list(edge_costs_weighted.values()))
    print(f"Edge cost statistics (with weights):")
    print(f"  Min: {np.min(costs_weighted):.6f}")
    print(f"  Q1: {np.percentile(costs_weighted, 25):.6f}")
    print(f"  Median: {np.median(costs_weighted):.6f}")
    print(f"  Q3: {np.percentile(costs_weighted, 75):.6f}")
    print(f"  Mean: {np.mean(costs_weighted):.6f}")
    print(f"  Max: {np.max(costs_weighted):.6f}")
    
    # How many edges are cheap vs expensive?
    print(f"\nEdge cost distribution:")
    cost_bins = [0.01, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2]
    for i in range(len(cost_bins)-1):
        count = np.sum((costs_weighted >= cost_bins[i]) & (costs_weighted < cost_bins[i+1]))
        pct = 100 * count / len(costs_weighted)
        print(f"  Cost in [{cost_bins[i]:.2f}, {cost_bins[i+1]:.2f}): {count:5d} edges ({pct:5.1f}%)")
    
    count = np.sum(costs_weighted >= cost_bins[-1])
    pct = 100 * count / len(costs_weighted)
    print(f"  Cost >= {cost_bins[-1]:.2f}: {count:5d} edges ({pct:5.1f}%)")
    
    # Critical insight: profitability analysis
    print(f"\n[5] PROFITABILITY ANALYSIS: Can We Build Modules?")
    print("-"*80)
    
    print(f"\nTo add N nodes to a module while connecting them with N-1 edges:")
    print(f"  Total prize needed: N × (average prize) = N × {np.mean(prizes_array):.4f}")
    print(f"  Total cost for edges: (N-1) × (cost per edge)")
    print(f"\nBreak-even analysis for uniform costs (edge_cost = {base_cost:.4f}):")
    
    for n_nodes in [2, 3, 4, 5, 10, 20]:
        top_n_prizes = [sorted_by_prize[i][1] for i in range(min(n_nodes, len(sorted_by_prize)))]
        avg_prize_n = np.mean(top_n_prizes)
        total_prize = avg_prize_n * n_nodes
        total_cost = (n_nodes - 1) * base_cost
        net = total_prize - total_cost
        print(f"  {n_nodes} nodes: prize={total_prize:.4f}, cost={total_cost:.4f}, net={net:+.4f} {'✓ PROFITABLE' if net > 0 else '✗ UNPROFITABLE'}")
    
    # Weighted costs analysis
    print(f"\nWith weighted edge costs (mean = {np.mean(costs_weighted):.4f}):")
    avg_cost_weighted = np.mean(costs_weighted)
    
    for n_nodes in [2, 3, 4, 5, 10, 20]:
        top_n_prizes = [sorted_by_prize[i][1] for i in range(min(n_nodes, len(sorted_by_prize)))]
        avg_prize_n = np.mean(top_n_prizes)
        total_prize = avg_prize_n * n_nodes
        total_cost = (n_nodes - 1) * avg_cost_weighted
        net = total_prize - total_cost
        print(f"  {n_nodes} nodes: prize={total_prize:.4f}, cost={total_cost:.4f}, net={net:+.4f} {'✓ PROFITABLE' if net > 0 else '✗ UNPROFITABLE'}")
    
    # Prize-to-cost comparison visualization
    print(f"\n[6] PRIZE vs COST COMPARISON")
    print("-"*80)
    
    print(f"\nTop prizes relative to cost:")
    print(f"{'Rank':<6} {'Prize':<10} {'Cost':<10} {'Prize/Cost':<12} {'Status':<15}")
    print(f"-" * 55)
    
    for i, (node, prize) in enumerate(sorted_by_prize[:15], 1):
        ratio = prize / base_cost if base_cost > 0 else 0
        status = "VALUABLE" if ratio >= 1.0 else "WEAK"
        cost_for_connection = base_cost
        print(f"{i:<6} {prize:<10.4f} {cost_for_connection:<10.4f} {ratio:<12.2f}x {status:<15}")
    
    print("\n" + "="*80)
    print("KEY INSIGHTS")
    print("="*80)
    print(f"""
1. Average node prize: {np.mean(prizes_array):.4f}
   - This is also the base_cost (uniform edge cost)
   - 1:1 ratio: each edge costs as much as average node's value

2. Top node prize: {sorted_by_prize[0][1]:.4f}
   - This is {sorted_by_prize[0][1] / base_cost:.2f}x the cost of an edge
   - Adding this top node + 1 edge: gain {sorted_by_prize[0][1]:.4f} - lose {base_cost:.4f} = {sorted_by_prize[0][1] - base_cost:+.4f}

3. To build a profitable 3-node module:
   - Need total prize >= 2 × edge cost = 2 × {base_cost:.4f} = {2 * base_cost:.4f}
   - Top 3 nodes total: {sum(sorted_by_prize[i][1] for i in range(3)):.4f}
   - This is WHY PCST selects ~3 nodes: diminishing returns beyond that

4. With weighted edges (mean cost: {np.mean(costs_weighted):.4f}):
   - Some edges are cheaper, some more expensive
   - On average still {np.mean(costs_weighted) / base_cost:.2f}x the base cost
""")

if __name__ == '__main__':
    main()

