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
    c0: float = 0.1
) -> Dict[Tuple[int, int], float]:
    """
    Prepare edge costs for PCST solver using author's recommended formula.
    
    Formula (when edge weights present): cost = base_cost * (1 - weight) + c0
    
    Where:
    - base_cost: typically -minimumScore (positive value)
    - weight: edge weight in [0, 1] (normalized)
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
    c0 : float, optional
        Minimum cost to prevent zeros when using edge weights (default: 0.1).
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
    # No edge weights: use R implementation (base_cost only)
    if edge_weight_attr is None:
        return {edge: base_cost for edge in network.edges()}
    
    # Edge weights provided: calculate weighted costs
    
    # Collect weights from edges
    weights_dict = {}
    for u, v in network.edges():
        if edge_weight_attr in network[u][v]:
            weights_dict[(u, v)] = network[u][v][edge_weight_attr]
    
    # No weights found: fall back to R implementation
    if not weights_dict:
        warnings.warn(f"No edges have attribute '{edge_weight_attr}'. "
                     "Using uniform costs (R implementation).")
        return {edge: base_cost for edge in network.edges()}
    
    # Normalize weights to [0, 1] if needed
    weights = np.array(list(weights_dict.values()))
    min_w, max_w = weights.min(), weights.max()
    
    if max_w > min_w:
        scale_factor = 1.0 / (max_w - min_w)
        offset = -min_w * scale_factor
    else:
        scale_factor = 0.0
        offset = 0.5
    
    # Calculate weighted costs: cost = base_cost * (1 - weight) + c0
    edge_costs = {}
    for edge, weight in weights_dict.items():
        weight_norm = weight * scale_factor + offset if max_w > min_w else 0.5
        edge_costs[edge] = base_cost * (1 - weight_norm) + c0
    
    # Fill edges without weights (use base_cost, matching R for unweighted edges)
    for u, v in network.edges():
        if (u, v) not in edge_costs and (v, u) not in edge_costs:
            edge_costs[(u, v)] = base_cost
    
    return edge_costs


def solve_pcst(
    network: nx.Graph,
    node_prizes: Dict[str, float],
    edge_costs: Dict[Tuple[str, str], float],
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
    edge_costs : Dict[Tuple[str, str], float]
        Dictionary mapping edge tuples to costs.
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
    
    # Prepare edges and costs - pre-allocate arrays
    num_edges = network.number_of_edges()
    edges = np.zeros((num_edges, 2), dtype=int)
    costs = np.zeros(num_edges, dtype=float)
    
    for i, (u, v) in enumerate(network.edges()):
        edges[i] = [node_to_idx[u], node_to_idx[v]]
        cost = edge_costs.get((u, v)) or edge_costs.get((v, u), 0.0)
        costs[i] = cost
    
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
    
    # Reconstruct solution from edges and vertices
    solution_indices = set()
    
    # Add from edges_in_solution
    for edge_idx in edges_in_solution:
        if 0 <= edge_idx < len(edges):
            u, v = edges[edge_idx]
            solution_indices.add(u)
            solution_indices.add(v)
    
    # Add from vertices
    for v_idx in vertices:
        if 0 <= v_idx < len(nodes):
            solution_indices.add(v_idx)
    
    return [idx_to_node[idx] for idx in solution_indices]


def detect_functional_module_core(
    network: nx.Graph,
    node_scores: Dict[str, float],
    edge_weight_attr: Optional[str] = None,
    c0: float = 0.1
) -> nx.Graph:
    """
    Core function to detect functional module using PCST.
    
    This function:
    1. Calculates prizes from node scores (prize = score - min_score)
    2. Prepares edge costs (base_cost = mean of prizes for balanced cost/prize ratio)
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
    c0 : float, optional
        Minimum edge cost to prevent zeros
        
    Returns
    -------
    nx.Graph
        Subgraph representing the functional module
    """
    if not node_scores:
        empty_graph = nx.Graph()
        empty_graph.graph['empty_solution'] = True
        empty_graph.graph['reason'] = 'no_node_scores'
        return empty_graph
    
    # Calculate min_score and prizes (shifted scores: prize = score - min_score)
    min_score = min(node_scores.values())
    prizes = {node: score - min_score for node, score in node_scores.items()}
    
    # Calculate base_cost as mean of prizes for balanced cost/prize ratio
    # This creates a 1:1 relationship: edge cost ≈ typical node prize
    base_cost = np.mean(list(prizes.values()))
    
    # Prepare edge costs
    edge_costs = prepare_edge_costs(
        network,
        base_cost,
        edge_weight_attr,
        c0
    )
    
    # Solve PCST
    solution_nodes = solve_pcst(
        network,
        prizes,
        edge_costs
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
        empty_graph = nx.Graph()
        empty_graph.graph['empty_solution'] = True
        empty_graph.graph['reason'] = 'pcst_returned_empty'
        return empty_graph
    
    # Create subgraph
    subgraph = network.subgraph(solution_nodes).copy()
    
    # Add scores and prizes as attributes
    for node in subgraph.nodes():
        if node in node_scores:
            subgraph.nodes[node]['score'] = node_scores[node]
            subgraph.nodes[node]['prize'] = prizes[node]
    
    return subgraph

