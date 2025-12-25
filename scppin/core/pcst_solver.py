"""Prize-Collecting Steiner Tree solver integration."""

import numpy as np
import igraph as ig
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
    network: ig.Graph,
    base_cost: float,
    edge_weight_attr: Optional[str] = None,
    c0: float = 0.01,
    normalization: Optional[str] = 'minmax'
) -> Dict[Tuple[str, str], float]:
    """
    Prepare edge costs for PCST solver from edge weights.
    
    This function normalizes edge weights (if provided) and converts them to costs.
    Formula (when edge weights present): cost = base_cost * (1 - weight) + c0
    
    Where:
    - base_cost: absolute value of minimum node score (positive value)
    - weight: edge weight (normalized to [0, 1] range internally, or used directly if normalization=None)
    - c0: minimum cost to prevent zeros (only used when edge weights are present)
    
    Higher edge weights → Lower costs → More likely to include edge
    
    When no edge weights are set, matches R implementation: cost = base_cost
    
    Parameters
    ----------
    network : ig.Graph
        Network with optional edge weights (will be normalized internally if provided and normalization is not None)
    base_cost : float
        Base cost for edges (absolute value of minimum node score, positive)
    edge_weight_attr : str, optional
        Name of edge attribute containing weights. If provided, weights will be
        normalized using the specified normalization method (unless normalization=None).
    c0 : float, optional
        Minimum cost to prevent zeros when using edge weights (default: 0.01).
        Ignored when edge_weight_attr is None.
    normalization : str, optional
        Normalization method for edge weights: 'minmax', 'log1p', 'power', or None
        (default: 'minmax'). If None, uses weights directly without normalization
        (assumes weights are already in [0, 1] range). Only used when edge_weight_attr is provided.
        
    Returns
    -------
    Dict[Tuple[str, str], float]
        Dictionary mapping edge tuples to costs
        
    References
    ----------
    Formula recommended by Florian Klimm (scPPIN author) in response to
    GitHub Issue #10: https://github.com/floklimm/scPPIN/issues/10
    """
    
    # Get edge list with node names
    edge_list = network.get_edgelist()
    node_names = network.vs['name']
    
    # No edge weights: use R implementation (base_cost only)
    if edge_weight_attr is None:
        return {
            (node_names[u], node_names[v]): base_cost
            for u, v in edge_list
        }
    
    # Edge weights provided: normalize and calculate weighted costs
    # Check if weights exist
    try:
        weights_raw = network.es[edge_weight_attr]
        has_weights = any(w is not None for w in weights_raw)
    except (KeyError, TypeError):
        has_weights = False
    
    if not has_weights:
        warnings.warn(
            f"Edge attribute '{edge_weight_attr}' not found. Using uniform edge costs.",
            UserWarning,
            stacklevel=2
        )
        return {
            (node_names[u], node_names[v]): base_cost
            for u, v in edge_list
        }
    
    # If normalization is None, use weights directly without normalization
    if normalization is None:
        norm_weights_raw = weights_raw
    else:
        # Normalize edge weights in-place
        normalize_edge_weights(
            network,
            edge_weight_attr,
            method=normalization,
            output_attr=f'{edge_weight_attr}_norm'
        )
        norm_attr = f'{edge_weight_attr}_norm'
        norm_weights_raw = network.es[norm_attr]
    
    # Calculate costs - handle None values
    edge_costs = {}
    for (u, v), weight in zip(edge_list, norm_weights_raw):
        if weight is not None:
            cost = (base_cost * (1 - weight)) + c0
        else:
            cost = base_cost
        edge_costs[(node_names[u], node_names[v])] = float(cost)
    
    return edge_costs


def solve_pcst(
    network: ig.Graph,
    node_prizes: Dict[str, float],
    edge_costs: Dict[Tuple[str, str], float],
    root: Optional[str] = None,
    num_clusters: int = 1,
    pruning: str = 'gw',
    verbosity: int = 0
) -> ig.Graph:
    """
    Solve Prize-Collecting Steiner Tree problem.
    
    Parameters
    ----------
    network : ig.Graph
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
    ig.Graph
        Induced subgraph from solution vertices. Includes all edges between
        solution nodes, with edge attribute 'in_solution' (bool) marking which
        edges were selected by the PCST algorithm.
        
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
    
    # Create name to index mapping
    name_to_idx = {network.vs[i]['name']: i for i in range(network.vcount())}
    idx_to_name = {i: network.vs[i]['name'] for i in range(network.vcount())}
    
    # Prepare prizes array
    prizes = np.zeros(network.vcount())
    for node, prize in node_prizes.items():
        if node in name_to_idx:
            prizes[name_to_idx[node]] = prize
    
    # Get edges directly as integer pairs (igraph already uses integers)
    edges = np.array(network.get_edgelist(), dtype=np.int64)
    
    # Get costs - vectorized lookup
    node_names = network.vs['name']
    costs = np.array([
        edge_costs.get((node_names[u], node_names[v])) or 
        edge_costs.get((node_names[v], node_names[u]), 0.0)
        for u, v in edges
    ])
    
    # Prepare root
    if root is not None:
        if root not in name_to_idx:
            raise ValueError(f"Root node '{root}' not in network")
        root_idx = name_to_idx[root]
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
    
    if len(pcst_nodes) > 0:
        solution_vertices = [int(v) for v in pcst_nodes if 0 <= v < network.vcount()]
        if solution_vertices:
            subgraph = network.induced_subgraph(solution_vertices)
            
            # Build solution edges set using names
            node_names = network.vs['name']
            solution_edges = {
                frozenset((node_names[edges[eid][0]], node_names[edges[eid][1]]))
                for eid in pcst_edges
                if 0 <= eid < len(edges)
            }
            
            # Batch mark edges (single attribute access)
            subgraph_names = subgraph.vs['name']
            subgraph.es['in_solution'] = [
                frozenset((subgraph_names[e.source], subgraph_names[e.target])) in solution_edges
                for e in subgraph.es
            ]
        else:
            subgraph = ig.Graph()
    else:
        subgraph = ig.Graph()
    
    return subgraph


def detect_functional_module_core(
    network: ig.Graph,
    node_scores: Dict[str, float],
    edge_weight_attr: Optional[str] = None,
    c0: float = 0.01,
    normalization: str = 'minmax',
    num_clusters: int = 1,
    pruning: str = 'gw',
    use_max_prize_root: bool = False
) -> ig.Graph:
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
    network : ig.Graph
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
    ig.Graph
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
        module_subgraph = ig.Graph()
        # Store metadata (igraph doesn't have direct graph attributes, use a workaround)
        # We'll check module_subgraph.vcount() == 0 to detect empty solutions
        return module_subgraph
    
    # Calculate prizes from node scores (prize = score - min_score)
    min_score = min(node_scores.values())
    prizes = {node: score - min_score for node, score in node_scores.items()}
    
    # Calculate base_cost as absolute value of maximum prize
    # This encourages connectivity by reducing edge costs relative to prizes
    
    # base_cost = abs(min_score) if min_score else np.quantile(list(prizes.values()), 0.1)
    # base_cost = abs(np.median(list(prizes.values())))
    base_cost = abs(max(prizes.values()))
    
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
    if module_subgraph.vcount() == 0:
        warnings.warn(
            f'PCST returned empty solution ({network.vcount()} nodes, '
            f'{len(node_scores)} scored). Try relaxing FDR or checking edge costs.'
        )
        return module_subgraph
    
    # Add node scores and prizes as attributes
    score_list = [None] * module_subgraph.vcount()
    prize_list = [None] * module_subgraph.vcount()
    
    for v in module_subgraph.vs:
        node_name = v['name']
        if node_name in node_scores:
            score_list[v.index] = node_scores[node_name]
        if node_name in prizes:
            prize_list[v.index] = prizes[node_name]
    
    module_subgraph.vs['score'] = score_list
    module_subgraph.vs['prize'] = prize_list
    
    return module_subgraph
    

