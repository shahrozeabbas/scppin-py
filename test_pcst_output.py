"""
Test script to investigate pcst_fast output format and behavior.
"""

import warnings
warnings.filterwarnings('ignore')
import pandas as pd
import numpy as np
import networkx as nx
from scppin.core import filter_network_by_pvalues
from scppin.core.scoring import shift_scores_for_pcst, compute_node_scores
from scppin.core import fit_bum
from scppin.core.pcst_solver import prepare_edge_costs
import pcst_fast

# Load data
print("Loading data...")
network_df = pd.read_csv('tests/denoised_mg_ad_fava_network.tsv', sep='\t')
pvalues_df = pd.read_csv('tests/magma_gene_symbol_results.tsv', sep='\t')
pvalues = dict(zip(pvalues_df['GENE'], pvalues_df['P']))

# Build network
edges_list = [(row['Protein_1'], row['Protein_2']) for _, row in network_df.iterrows()]
network = nx.Graph(edges_list)
network_filtered = filter_network_by_pvalues(network, pvalues, missing_data_score=False)

print(f"Network: {network_filtered.number_of_nodes()} nodes, {network_filtered.number_of_edges()} edges")

# Compute scores and prizes
pvalues_array = np.array(list(pvalues.values()))
lambda_param, alpha, success = fit_bum(pvalues_array)
node_scores = compute_node_scores(network_filtered, pvalues, lambda_param, alpha, 0.01, missing_data_score=False)
node_prizes = shift_scores_for_pcst(node_scores)

# Calculate base_cost (using mean prize)
prizes_array = np.array(list(node_prizes.values()))
base_cost = float(np.mean(prizes_array))
print(f"\nBase cost (mean prize): {base_cost:.6f}")

# Prepare for pcst_fast
nodes = list(network_filtered.nodes())
node_to_idx = {node: idx for idx, node in enumerate(nodes)}
idx_to_node = {idx: node for node, idx in node_to_idx.items()}

prizes = np.zeros(len(nodes))
for node, prize in node_prizes.items():
    if node in node_to_idx:
        prizes[node_to_idx[node]] = prize

edges = []
costs = []
for u, v in network_filtered.edges():
    u_idx = node_to_idx[u]
    v_idx = node_to_idx[v]
    edges.append((u_idx, v_idx))
    costs.append(base_cost)

edges = np.array(edges, dtype=int)
costs = np.array(costs, dtype=float)

print(f"\nPrepared for pcst_fast:")
print(f"  Edges shape: {edges.shape}, dtype: {edges.dtype}")
print(f"  Prizes shape: {prizes.shape}, dtype: {prizes.dtype}, min={prizes.min():.6f}, max={prizes.max():.6f}, mean={prizes.mean():.6f}")
print(f"  Costs shape: {costs.shape}, dtype: {costs.dtype}, min={costs.min():.6f}, max={costs.max():.6f}, mean={costs.mean():.6f}")
print(f"  Non-zero prizes: {np.sum(prizes > 0)}")

# Call pcst_fast
print("\n" + "="*80)
print("CALLING pcst_fast")
print("="*80)
vertices, edges_in_solution = pcst_fast.pcst_fast(
    edges,
    prizes,
    costs,
    -1,  # root (let solver choose)
    1,   # num_clusters
    'gw', # pruning
    0    # verbosity
)

print(f"\npcst_fast returned:")
print(f"  vertices type: {type(vertices)}, shape: {vertices.shape}, dtype: {vertices.dtype}")
print(f"  edges_in_solution type: {type(edges_in_solution)}, shape: {edges_in_solution.shape}, dtype: {edges_in_solution.dtype}")

# Analyze vertices array
print(f"\nVertices array analysis:")
print(f"  Length: {len(vertices)}")
print(f"  Unique values: {len(np.unique(vertices))}")
print(f"  Unique values (first 20): {np.unique(vertices)[:20]}")
print(f"  Min: {vertices.min()}, Max: {vertices.max()}")
print(f"  First 50 values: {vertices[:50]}")
print(f"  Last 50 values: {vertices[-50:]}")

# Check if vertices are valid indices
valid_vertices = vertices[(vertices >= 0) & (vertices < len(nodes))]
print(f"  Valid vertex indices (0 <= idx < {len(nodes)}): {len(valid_vertices)}")
if len(valid_vertices) > 0:
    print(f"  Valid vertices (first 20): {valid_vertices[:20]}")

# Analyze edges_in_solution array
print(f"\nEdges_in_solution array analysis:")
print(f"  Length: {len(edges_in_solution)}")
print(f"  Unique values: {len(np.unique(edges_in_solution))}")
print(f"  Unique values (first 20): {np.unique(edges_in_solution)[:20]}")
print(f"  Min: {edges_in_solution.min()}, Max: {edges_in_solution.max()}")
print(f"  First 50 values: {edges_in_solution[:50]}")
print(f"  Last 50 values: {edges_in_solution[-50:]}")

# Check if edge indices are valid
valid_edges = edges_in_solution[(edges_in_solution >= 0) & (edges_in_solution < len(edges))]
print(f"  Valid edge indices (0 <= idx < {len(edges)}): {len(valid_edges)}")
if len(valid_edges) > 0:
    print(f"  Valid edge indices (first 20): {valid_edges[:20]}")

# Try to reconstruct solution using current method
print("\n" + "="*80)
print("RECONSTRUCTING SOLUTION (current method)")
print("="*80)
solution_vertex_indices = set()

for edge_idx in edges_in_solution:
    edge_idx_int = int(edge_idx)
    if 0 <= edge_idx_int < len(edges):
        u_idx, v_idx = edges[edge_idx_int]
        solution_vertex_indices.add(int(u_idx))
        solution_vertex_indices.add(int(v_idx))

print(f"Vertices reconstructed from edges: {len(solution_vertex_indices)}")
if len(solution_vertex_indices) > 0:
    print(f"  Vertex indices: {sorted(list(solution_vertex_indices))[:20]}")
    print(f"  Node names: {[idx_to_node[idx] for idx in sorted(list(solution_vertex_indices))[:10]]}")

# Try alternative: use vertices array directly
print("\n" + "="*80)
print("ALTERNATIVE: Using vertices array directly")
print("="*80)
vertices_from_array = set()
for v_idx in vertices:
    v_idx_int = int(v_idx)
    if 0 <= v_idx_int < len(prizes) and prizes[v_idx_int] > 0:
        vertices_from_array.add(v_idx_int)

print(f"Vertices from array (with positive prizes): {len(vertices_from_array)}")
if len(vertices_from_array) > 0:
    print(f"  Vertex indices: {sorted(list(vertices_from_array))[:20]}")
    print(f"  Node names: {[idx_to_node[idx] for idx in sorted(list(vertices_from_array))[:10]]}")

# Try alternative: use all unique vertices
print("\n" + "="*80)
print("ALTERNATIVE: Using all unique vertices")
print("="*80)
unique_vertices = np.unique(vertices)
vertices_from_unique = set()
for v_idx in unique_vertices:
    v_idx_int = int(v_idx)
    if 0 <= v_idx_int < len(prizes):
        vertices_from_unique.add(v_idx_int)

print(f"Unique vertices (all): {len(vertices_from_unique)}")
if len(vertices_from_unique) > 0:
    print(f"  Vertex indices: {sorted(list(vertices_from_unique))[:20]}")
    print(f"  Node names: {[idx_to_node[idx] for idx in sorted(list(vertices_from_unique))[:10]]}")

# Check what edges connect the vertices from array
if len(vertices_from_array) > 0:
    print("\n" + "="*80)
    print("CHECKING EDGES BETWEEN VERTICES FROM ARRAY")
    print("="*80)
    vertices_list = list(vertices_from_array)
    edges_between = []
    for i, u_idx in enumerate(vertices_list):
        for v_idx in vertices_list[i+1:]:
            # Check if edge exists
            edge_found = False
            for edge_idx, (eu, ev) in enumerate(edges):
                if (eu == u_idx and ev == v_idx) or (eu == v_idx and ev == u_idx):
                    edges_between.append((edge_idx, u_idx, v_idx))
                    edge_found = True
                    break
    print(f"Edges found between vertices: {len(edges_between)}")
    if len(edges_between) > 0:
        print(f"  Sample edges: {edges_between[:10]}")

# Test with a simple case
print("\n" + "="*80)
print("TESTING WITH SIMPLE CASE")
print("="*80)
# Create a simple 5-node star graph
test_edges = np.array([
    [0, 1],
    [0, 2],
    [0, 3],
    [0, 4]
], dtype=int)
test_prizes = np.array([10.0, 5.0, 5.0, 5.0, 5.0], dtype=float)
test_costs = np.array([1.0, 1.0, 1.0, 1.0], dtype=float)

print(f"Test case:")
print(f"  Edges: {test_edges}")
print(f"  Prizes: {test_prizes}")
print(f"  Costs: {test_costs}")

test_vertices, test_edges_sol = pcst_fast.pcst_fast(
    test_edges,
    test_prizes,
    test_costs,
    -1,
    1,
    'gw',
    0
)

print(f"\nSimple test results:")
print(f"  vertices: {test_vertices}, unique: {np.unique(test_vertices)}")
print(f"  edges_in_solution: {test_edges_sol}, unique: {np.unique(test_edges_sol)}")

# Reconstruct from edges
test_solution = set()
for edge_idx in test_edges_sol:
    if 0 <= edge_idx < len(test_edges):
        u, v = test_edges[edge_idx]
        test_solution.add(u)
        test_solution.add(v)
print(f"  Solution from edges: {sorted(test_solution)}")

# Reconstruct from vertices
test_solution2 = set(test_vertices)
print(f"  Solution from vertices: {sorted(test_solution2)}")

