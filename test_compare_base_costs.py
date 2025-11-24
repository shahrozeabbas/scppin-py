#!/usr/bin/env python
"""
Compare results between old and new base_cost calculations.
"""

import pandas as pd
import networkx as nx
import numpy as np
from scppin.core.scoring import compute_node_scores
from scppin.core.pcst_solver import prepare_edge_costs, solve_pcst
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
    print("COMPARING OLD vs NEW BASE_COST CALCULATIONS")
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
    
    # OLD base_cost
    base_cost_old = -min_score
    
    # NEW base_cost
    base_cost_new = np.mean(list(prizes.values()))
    
    print(f"\n[1] BASE COST VALUES")
    print("-"*80)
    print(f"Old (R impl):   base_cost = -min_score = {base_cost_old:.6f}")
    print(f"New (Strategy5): base_cost = mean_prize = {base_cost_new:.6f}")
    print(f"Ratio:          {base_cost_old / base_cost_new:.2f}x")
    
    # Test with uniform edge weights (no c0 issues)
    print(f"\n[2] SOLVING WITH UNIFORM EDGE COSTS (c0=0.01)")
    print("-"*80)
    
    for name, bc in [("OLD", base_cost_old), ("NEW", base_cost_new)]:
        edge_costs = prepare_edge_costs(
            network_filtered,
            bc,
            edge_weight_attr=None,
            c0=0.01
        )
        
        solution_nodes = solve_pcst(
            network_filtered,
            prizes,
            edge_costs
        )
        
        total_prize = sum(prizes[node] for node in solution_nodes)
        total_cost = sum(edge_costs[(u, v)] if (u, v) in edge_costs else edge_costs[(v, u)] 
                        for u, v in network_filtered.subgraph(solution_nodes).edges())
        net_value = total_prize - total_cost
        
        print(f"\n{name} approach:")
        print(f"  Module size: {len(solution_nodes)} nodes, {network_filtered.subgraph(solution_nodes).number_of_edges()} edges")
        print(f"  Total prize: {total_prize:.6f}")
        print(f"  Total cost:  {total_cost:.6f}")
        print(f"  Net value:   {net_value:.6f}")
        
        if len(solution_nodes) <= 10:
            print(f"  Nodes: {sorted(solution_nodes)}")
    
    # Test with weighted edge costs
    print(f"\n[3] SOLVING WITH WEIGHTED EDGE COSTS (c0=0.01)")
    print("-"*80)
    
    for name, bc in [("OLD", base_cost_old), ("NEW", base_cost_new)]:
        edge_costs = prepare_edge_costs(
            network_filtered,
            bc,
            edge_weight_attr='weight',
            c0=0.01
        )
        
        solution_nodes = solve_pcst(
            network_filtered,
            prizes,
            edge_costs
        )
        
        total_prize = sum(prizes[node] for node in solution_nodes)
        total_cost = sum(edge_costs[(u, v)] if (u, v) in edge_costs else edge_costs[(v, u)] 
                        for u, v in network_filtered.subgraph(solution_nodes).edges())
        net_value = total_prize - total_cost
        
        print(f"\n{name} approach:")
        print(f"  Module size: {len(solution_nodes)} nodes, {network_filtered.subgraph(solution_nodes).number_of_edges()} edges")
        print(f"  Total prize: {total_prize:.6f}")
        print(f"  Total cost:  {total_cost:.6f}")
        print(f"  Net value:   {net_value:.6f}")
        
        if len(solution_nodes) <= 10:
            print(f"  Nodes: {sorted(solution_nodes)}")
    
    # Analyze edge cost distribution
    print(f"\n[4] EDGE COST DISTRIBUTION COMPARISON")
    print("-"*80)
    
    edge_costs_old = prepare_edge_costs(network_filtered, base_cost_old, edge_weight_attr=None, c0=0.01)
    edge_costs_new = prepare_edge_costs(network_filtered, base_cost_new, edge_weight_attr=None, c0=0.01)
    
    print(f"\nUniform costs (no weights):")
    print(f"  Old: {np.mean(list(edge_costs_old.values())):.6f} (all edges = base_cost)")
    print(f"  New: {np.mean(list(edge_costs_new.values())):.6f} (all edges = base_cost)")
    
    edge_costs_old_weighted = prepare_edge_costs(network_filtered, base_cost_old, edge_weight_attr='weight', c0=0.01)
    edge_costs_new_weighted = prepare_edge_costs(network_filtered, base_cost_new, edge_weight_attr='weight', c0=0.01)
    
    print(f"\nWeighted costs (with edge weights):")
    print(f"  Old: min={np.min(list(edge_costs_old_weighted.values())):.6f}, max={np.max(list(edge_costs_old_weighted.values())):.6f}, mean={np.mean(list(edge_costs_old_weighted.values())):.6f}")
    print(f"  New: min={np.min(list(edge_costs_new_weighted.values())):.6f}, max={np.max(list(edge_costs_new_weighted.values())):.6f}, mean={np.mean(list(edge_costs_new_weighted.values())):.6f}")
    
    print("\n" + "="*80)

if __name__ == '__main__':
    main()

