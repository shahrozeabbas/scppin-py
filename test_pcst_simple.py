#!/usr/bin/env python
"""
Test pcst_fast with a simple manual example to understand the API.
"""

import numpy as np
import pcst_fast

print("="*80)
print("SIMPLE pcst_fast TEST")
print("="*80)

# Create a simple graph: 5 nodes, fully connected
# Nodes: 0, 1, 2, 3, 4
# We want to collect nodes with high prizes

nodes = np.array([0, 1, 2, 3, 4])
prizes = np.array([0.0, 10.0, 5.0, 2.0, 1.0])  # Node 1 has highest prize

# Create edges (all pairs)
edges = []
costs = []
cost_default = 1.0

for i in range(5):
    for j in range(i+1, 5):
        edges.append([i, j])
        costs.append(cost_default)

edges = np.array(edges, dtype=int)
costs = np.array(costs, dtype=float)

print(f"\nSimple graph:")
print(f"  Nodes: {nodes}")
print(f"  Prizes: {prizes}")
print(f"  Edges: {len(edges)}")
print(f"  Edge cost: {cost_default}")

# Call pcst_fast
print(f"\nCalling pcst_fast...")
vertices, edges_in_solution = pcst_fast.pcst_fast(
    edges,
    prizes,
    costs,
    -1,  # root
    1,   # num_clusters
    "gw",
    1    # verbosity
)

print(f"\nReturn values:")
print(f"  vertices: {vertices}")
print(f"  edges_in_solution: {edges_in_solution}")

# Interpret
print(f"\nInterpretation:")
print(f"  vertices should contain node indices? len={len(vertices)}")
print(f"  edges_in_solution should contain edge indices? len={len(edges_in_solution)}")

# Try to reconstruct
solution_nodes = set()
for edge_idx in edges_in_solution:
    if 0 <= edge_idx < len(edges):
        u, v = edges[edge_idx]
        solution_nodes.add(u)
        solution_nodes.add(v)

for v_idx in vertices:
    if 0 <= v_idx < len(nodes):
        solution_nodes.add(v_idx)

print(f"\nReconstructed solution: {sorted(solution_nodes)}")
print(f"  Prizes: {[prizes[i] for i in sorted(solution_nodes)]}")
print(f"  Total prize: {sum(prizes[i] for i in sorted(solution_nodes))}")
print(f"  Num edges in solution: {len([e for e in edges_in_solution if 0 <= e < len(edges)])}")

print("\n" + "="*80)

