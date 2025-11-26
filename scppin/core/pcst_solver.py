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

from .edge_weights import normalize_edge_weights


def prepare_edge_costs(
    network: nx.Graph,
    base_cost: float,
    edge_weight_attr: Optional[str] = None,
    c0: float = 0.01,
    normalization: str = 'minmax'
) -> Dict[Tuple[str, str], float]:
    """
    Prepare edge costs for PCST solver from edge weights.
    
    This function normalizes edge weights (if provided) and converts them to costs.
    Formula (when edge weights present): cost = base_cost * (1 - weight) + c0
    
    Where:
    - base_cost: absolute value of minimum node score (positive value)
    - weight: edge weight (normalized to [0, 1] range internally)
    - c0: minimum cost to prevent zeros (only used when edge weights are present)
    
    Higher edge weights → Lower costs → More likely to include edge
    
    When no edge weights are set, matches R implementation: cost = base_cost
    
    Parameters
    ----------
    network : nx.Graph
        Network with optional edge weights (will be normalized internally if provided)
    base_cost : float
        Base cost for edges (absolute value of minimum node score, positive)
    edge_weight_attr : str, optional
        Name of edge attribute containing weights. If provided, weights will be
        normalized using the specified normalization method.
    c0 : float, optional
        Minimum cost to prevent zeros when using edge weights (default: 0.01).
        Ignored when edge_weight_attr is None.
    normalization : str, optional
        Normalization method for edge weights: 'minmax', 'log1p', or 'power'
        (default: 'minmax'). Only used when edge_weight_attr is provided.
        
    Returns
    -------
    Dict[Tuple[str, str], float]
        Dictionary mapping edge tuples to costs
        
    References
    ----------
    Formula recommended by Florian Klimm (scPPIN author) in response to
    GitHub Issue #10: https://github.com/floklimm/scPPIN/issues/10
    """
    
    # No edge weights: use R implementation (base_cost only)
    if edge_weight_attr is None:
        return {edge: base_cost for edge in network.edges()}
    
    # Edge weights provided: normalize and calculate weighted costs
    # Check if weights exist
    has_weights = any(
        edge_weight_attr in network[u][v]
        for u, v in network.edges()
    )
    
    if not has_weights:
        warnings.warn(f"No edges have attribute '{edge_weight_attr}'. "
                     "Using uniform costs (R implementation).")
        return {edge: base_cost for edge in network.edges()}
    
    # Normalize edge weights in-place
    normalize_edge_weights(
        network,
        edge_weight_attr,
        method=normalization,
        output_attr=f'{edge_weight_attr}_norm'
    )
    
    # Use normalized weights for cost calculation
    norm_attr = f'{edge_weight_attr}_norm'
    edge_costs = {}
    
    # Single pass: collect normalized weights and compute costs
    for u, v in network.edges():
        if norm_attr in network[u][v]:
            weight = network[u][v][norm_attr]
            cost = (base_cost * (1 - weight)) + c0
            edge_costs[(u, v)] = cost
        else:
            # Edge without normalized weight (shouldn't happen, but handle gracefully)
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
    idx_to_node = {v: k for k, v in node_to_idx.items()}
    
    # Prepare prizes array
    prizes = np.zeros(len(nodes))
    for node, prize in node_prizes.items():
        if node in node_to_idx:
            prizes[node_to_idx[node]] = prize
    
    # Prepare edges and costs - pre-allocate arrays
    # Use int64 to match pcst_fast expectations (C++ code expects int64_t)
    num_edges = network.number_of_edges()
    edges = np.zeros((num_edges, 2), dtype=np.int64)
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
    pcst_nodes, pcst_edges = pcst_fast.pcst_fast(
        edges,
        prizes,
        costs,
        root_idx,
        num_clusters,
        pruning,
        verbosity
    )
    
    module_edges = []

    for edge in pcst_edges:
        if 0 <= edge < len(edges):
            u, v = edges[edge]
            module_edges.append((idx_to_node[u], idx_to_node[v]))
    
    if module_edges:
        subgraph = network.edge_subgraph(module_edges)
    else:
        subgraph = nx.Graph()
        subgraph.graph['empty_solution'] = True
        subgraph.graph['reason'] = 'pcst_returned_empty'
    
    return subgraph


def detect_functional_module_core(
    network: nx.Graph,
    node_scores: Dict[str, float],
    edge_weight_attr: Optional[str] = None,
    c0: float = 0.01,
    normalization: str = 'minmax',
    num_clusters: int = 1,
    pruning: str = 'gw',
    use_max_prize_root: bool = False
) -> nx.Graph:
    """
    Core function to detect functional module using PCST optimization.
    
    This function:
    1. Calculates prizes from node scores (prize = score - min_score)
    2. Prepares edge costs from edge weights (base_cost = abs(min_score))
    3. Optionally selects root node (highest prize if use_max_prize_root=True)
    4. Solves PCST
    5. Returns subgraph
    
    Edge weights are normalized internally if provided.
    
    Parameters
    ----------
    network : nx.Graph
        Protein-protein interaction network
    node_scores : Dict[str, float]
        Node scores (can be negative)
    edge_weight_attr : str, optional
        Edge attribute for weights. If provided, weights will be normalized
        using the specified normalization method.
    c0 : float, optional
        Minimum edge cost to prevent zeros (default: 0.01)
    normalization : str, optional
        Normalization method for edge weights: 'minmax', 'log1p', or 'power'
        (default: 'minmax'). Only used when edge_weight_attr is provided.
    num_clusters : int, optional
        Number of connected components to return (default: 1)
    pruning : str, optional
        Pruning method: 'gw' (Goemans-Williamson) or 'strong' (default: 'gw')
    use_max_prize_root : bool, optional
        If True, use the node with highest prize as root. If False, let PCST
        algorithm choose root automatically (default: False)
        
    Returns
    -------
    nx.Graph
        Subgraph representing the functional module(s)
        
    Examples
    --------
    >>> # Use default PCST
    >>> module = detect_functional_module_core(network, node_scores)
    >>> 
    >>> # Use PCST with custom parameters
    >>> module = detect_functional_module_core(
    ...     network, node_scores,
    ...     num_clusters=1, pruning='gw', use_max_prize_root=True
    ... )
    """
    if not node_scores:
        module_subgraph = nx.Graph()
        module_subgraph.graph['empty_solution'] = True
        module_subgraph.graph['reason'] = 'no_node_scores'
        return module_subgraph
    
    # Calculate prizes from node scores (prize = score - min_score)
    min_score = min(node_scores.values())
    prizes = {node: score - min_score for node, score in node_scores.items()}
    
    # Calculate base_cost as absolute value of minimum node score
    # This encourages connectivity by reducing edge costs relative to prizes
    base_cost = abs(min_score) if min_score else np.quantile(list(prizes.values()), 0.1)
    
    # Prepare edge costs using centralized function
    edge_costs = prepare_edge_costs(
        network, base_cost, c0=c0,
        edge_weight_attr=edge_weight_attr,
        normalization=normalization
    )
    
    # Determine root node if requested
    root_node = None
    
    if use_max_prize_root and prizes:
        # Find node with highest prize
        root_node = max(prizes.items(), key=lambda x: x[1])[0]
    
    # Solve PCST
    module_subgraph = solve_pcst(
        network,
        prizes,
        edge_costs,
        root=root_node,
        num_clusters=num_clusters,
        pruning=pruning
    )
    
    # Handle empty solution
    if module_subgraph.graph.get('empty_solution', False):
        warnings.warn(
            "PCST solver returned empty solution. This may indicate:\n"
            "1) All node scores are negative after shifting (all prizes are zero)\n"
            "2) Edge costs are too high relative to node prizes\n"
            "3) Network structure prevents module formation\n"
            "4) FDR threshold may be too strict\n"
            f"Network had {network.number_of_nodes()} nodes, "
            f"{len(node_scores)} had scores."
        )
        return module_subgraph
    
    for node in module_subgraph.nodes():
        if node in node_scores:
            module_subgraph.nodes[node]['score'] = node_scores[node]
        if node in prizes:
            module_subgraph.nodes[node]['prize'] = prizes[node]
    
    return module_subgraph
    

