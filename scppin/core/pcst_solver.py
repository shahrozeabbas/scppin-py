"""Prize-Collecting Steiner Tree solver integration."""

import numpy as np
import networkx as nx
from typing import Dict, Tuple, Optional, List
import warnings

try:
    import pcst_fast
    PCST_AVAILABLE = True
except ImportError:
    PCST_AVAILABLE = False
    warnings.warn(
        "pcst_fast library not found. Install with: pip install pcst-fast\n"
        "PCST solver will not be available."
    )


def prepare_edge_costs(
    network: nx.Graph,
    base_cost: float,
    edge_weight_attr: Optional[str] = None,
    edge_weight_scale: float = 1.0,
    c0: Optional[float] = None
) -> Dict[Tuple[int, int], float]:
    """
    Prepare edge costs for PCST solver using author's recommended formula.
    
    Formula (when edge weights present): cost = base_cost * (1 - weight * scale) + c0
    
    Where:
    - base_cost: typically -minimumScore (positive value)
    - weight: edge weight in [0, 1]
    - scale: scaling factor for weight influence
    - c0: minimum cost to prevent zeros (only used when edge weights are present)
    
    Higher edge weights → Lower costs → More likely to include edge
    
    When no edge weights are set, matches R implementation: cost = base_cost
    
    Parameters
    ----------
    network : nx.Graph
        Network with optional edge weights
    base_cost : float
        Base cost for edges (typically -min_node_score, positive)
    edge_weight_attr : str, optional
        Name of edge attribute containing weights [0, 1]
    edge_weight_scale : float, optional
        Scaling factor for edge weight influence (default: 1.0)
    c0 : float, optional
        Minimum cost to prevent zeros when using edge weights (default: 0.1 * base_cost).
        Ignored when edge_weight_attr is None.
        
    Returns
    -------
    Dict[Tuple[int, int], float]
        Dictionary mapping edge tuples to costs
        
    References
    ----------
    Formula recommended by Florian Klimm (scPPIN author) in response to
    GitHub Issue #10: https://github.com/floklimm/scPPIN/issues/10
    """
    edge_costs = {}
    
    if edge_weight_attr is None:
        # Uniform edge costs - match R implementation exactly (no c0 added)
        for u, v in network.edges():
            edge_costs[(u, v)] = base_cost
    else:
        # Edge-weighted costs using author's formula
        # Set default c0 only when using edge weights
        if c0 is None:
            c0 = 0.1 * base_cost  # Default: 10% of base cost
        
        # Check if weights exist
        weights = []
        for u, v in network.edges():
            if edge_weight_attr in network[u][v]:
                weights.append(network[u][v][edge_weight_attr])
        
        if not weights:
            warnings.warn(f"No edges have attribute '{edge_weight_attr}'. "
                         "Using uniform costs.")
            for u, v in network.edges():
                edge_costs[(u, v)] = base_cost + c0
            return edge_costs
        
        # Check if already normalized to [0, 1]
        min_w = min(weights)
        max_w = max(weights)
        needs_normalization = (min_w < 0 or max_w > 1)
        
        for u, v in network.edges():
            if edge_weight_attr in network[u][v]:
                weight = network[u][v][edge_weight_attr]
                
                # Normalize if needed
                if needs_normalization and max_w > min_w:
                    weight_norm = (weight - min_w) / (max_w - min_w)
                elif needs_normalization:
                    weight_norm = 0.5  # All weights same
                else:
                    weight_norm = weight
                
                # Apply author's formula
                # cost = base_cost * (1 - weight * scale) + c0
                cost = base_cost * (1 - weight_norm * edge_weight_scale) + c0
                edge_costs[(u, v)] = cost
            else:
                # Edge doesn't have weight attribute, use base cost
                edge_costs[(u, v)] = base_cost + c0
    
    return edge_costs


def solve_pcst(
    network: nx.Graph,
    node_prizes: Dict[str, float],
    edge_costs: Optional[Dict[Tuple[str, str], float]] = None,
    root: Optional[str] = None,
    num_clusters: int = 1,
    pruning: str = 'gw',
    verbosity: int = 0
) -> List[str]:
    """
    Solve Prize-Collecting Steiner Tree problem.
    
    Parameters
    ----------
    network : nx.Graph
        Network graph
    node_prizes : Dict[str, float]
        Dictionary mapping node names to prizes (must be non-negative)
    edge_costs : Dict[Tuple[str, str], float], optional
        Dictionary mapping edge tuples to costs. If None, uses uniform costs.
    root : str, optional
        Root node (if None, solver chooses automatically)
    num_clusters : int, optional
        Number of connected components to return (default: 1)
    pruning : str, optional
        Pruning method: 'gw' (Goemans-Williamson) or 'strong' (default: 'gw')
    verbosity : int, optional
        Verbosity level (default: 0)
        
    Returns
    -------
    List[str]
        List of node names in the optimal solution
        
    Raises
    ------
    ImportError
        If pcst_fast is not installed
    ValueError
        If node prizes are negative or network is invalid
    """
    if not PCST_AVAILABLE:
        raise ImportError(
            "pcst_fast library is required for PCST solving.\n"
            "Install with: pip install pcst-fast"
        )
    
    # Validate inputs
    if not node_prizes:
        raise ValueError("node_prizes cannot be empty")
    
    if any(prize < 0 for prize in node_prizes.values()):
        raise ValueError("All node prizes must be non-negative")
    
    # Create node mapping (string names to integer indices)
    nodes = list(network.nodes())
    node_to_idx = {node: idx for idx, node in enumerate(nodes)}
    idx_to_node = {idx: node for node, idx in node_to_idx.items()}
    
    # Prepare prizes array
    prizes = np.zeros(len(nodes))
    for node, prize in node_prizes.items():
        if node in node_to_idx:
            prizes[node_to_idx[node]] = prize
    
    # Prepare edges and costs
    edges = []
    costs = []
    
    for u, v in network.edges():
        u_idx = node_to_idx[u]
        v_idx = node_to_idx[v]
        edges.append((u_idx, v_idx))
        
        if edge_costs is not None and (u, v) in edge_costs:
            cost = edge_costs[(u, v)]
        elif edge_costs is not None and (v, u) in edge_costs:
            cost = edge_costs[(v, u)]
        else:
            # Default cost (should not happen if edge_costs properly prepared)
            cost = 0.0
        
        costs.append(cost)
    
    edges = np.array(edges, dtype=int)
    costs = np.array(costs, dtype=float)
    
    # Prepare root
    if root is not None:
        if root not in node_to_idx:
            raise ValueError(f"Root node '{root}' not in network")
        root_idx = node_to_idx[root]
    else:
        root_idx = -1  # Let solver choose
    
    # Solve PCST
    vertices, edges_in_solution = pcst_fast.pcst_fast(
        edges,
        prizes,
        costs,
        root_idx,
        num_clusters,
        pruning,
        verbosity
    )
    
    # IMPORTANT: There appears to be an issue with pcst_fast where both vertices and 
    # edges_in_solution arrays contain repeated single values rather than arrays of 
    # different indices. The edges_in_solution array seems more reliable.
    # We reconstruct the solution from edges_in_solution to get all nodes in the solution.
    solution_vertex_indices = set()
    
    # Get unique edge indices (pcst_fast returns repeated values)
    unique_edge_indices = np.unique(edges_in_solution)
    
    # Add vertices from edges in solution
    for edge_idx in unique_edge_indices:
        edge_idx_int = int(edge_idx)
        if 0 <= edge_idx_int < len(edges):
            u_idx, v_idx = edges[edge_idx_int]
            solution_vertex_indices.add(int(u_idx))
            solution_vertex_indices.add(int(v_idx))
    
    # Also check vertices array for any isolated nodes (nodes with no edges in solution)
    # The vertices array may contain additional nodes not connected by edges
    unique_vertex_indices = np.unique(vertices)
    for v_idx in unique_vertex_indices:
        v_idx_int = int(v_idx)
        if 0 <= v_idx_int < len(nodes):
            # Only add if this node has a non-zero prize (isolated high-value nodes)
            if prizes[v_idx_int] > 0:
                solution_vertex_indices.add(v_idx_int)
    
    # Convert back to node names
    solution_nodes = [idx_to_node[idx] for idx in sorted(solution_vertex_indices)]
    
    return solution_nodes


def detect_functional_module_core(
    network: nx.Graph,
    node_scores: Dict[str, float],
    edge_weight_attr: Optional[str] = None,
    edge_weight_scale: float = 1.0,
    c0: Optional[float] = None
) -> nx.Graph:
    """
    Core function to detect functional module using PCST.
    
    This function:
    1. Shifts node scores to make them non-negative (prizes)
    2. Prepares edge costs (uniform or weighted)
    3. Solves PCST
    4. Returns subgraph
    
    Parameters
    ----------
    network : nx.Graph
        Protein-protein interaction network
    node_scores : Dict[str, float]
        Node scores (can be negative)
    edge_weight_attr : str, optional
        Edge attribute for weights
    edge_weight_scale : float, optional
        Scaling factor for edge weights
    c0 : float, optional
        Minimum edge cost to prevent zeros
        
    Returns
    -------
    nx.Graph
        Subgraph representing the functional module
    """
    from .scoring import shift_scores_for_pcst, get_minimum_score
    
    # Shift scores to make them non-negative (prizes)
    node_prizes = shift_scores_for_pcst(node_scores)
    
    # Calculate base_cost as mean prize to match prize scale
    # This ensures costs and prizes are on similar scales
    if node_prizes:
        prizes_array = np.array(list(node_prizes.values()))
        base_cost = float(np.mean(prizes_array))
    else:
        # Fallback to original method if no prizes
        min_score = get_minimum_score(node_scores)
        base_cost = -min_score
    
    # Prepare edge costs
    # Create mapping with node objects (not strings)
    node_mapping = {node: node for node in network.nodes()}
    edge_costs_mapped = prepare_edge_costs(
        network,
        base_cost,
        edge_weight_attr,
        edge_weight_scale,
        c0
    )
    
    # Solve PCST
    solution_nodes = solve_pcst(
        network,
        node_prizes,
        edge_costs_mapped
    )
    
    # Handle empty solution
    if not solution_nodes:
        warnings.warn(
            "PCST solver returned empty solution. This may indicate:\n"
            "1) All node scores are negative after shifting (all prizes are zero)\n"
            "2) Edge costs are too high relative to node prizes\n"
            "3) Network structure prevents module formation\n"
            "4) FDR threshold may be too strict\n"
            f"Network had {network.number_of_nodes()} nodes, "
            f"{len(node_scores)} had scores."
        )
        # Return empty graph with metadata
        empty_graph = nx.Graph()
        empty_graph.graph['empty_solution'] = True
        empty_graph.graph['reason'] = 'pcst_returned_empty'
        return empty_graph
    
    # Create subgraph
    subgraph = network.subgraph(solution_nodes).copy()
    
    # Add original scores as attributes
    for node in subgraph.nodes():
        if node in node_scores:
            subgraph.nodes[node]['score'] = node_scores[node]
            subgraph.nodes[node]['prize'] = node_prizes[node]
    
    return subgraph

