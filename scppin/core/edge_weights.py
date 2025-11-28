"""Edge weight computation and normalization utilities."""

import numpy as np
import igraph as ig
from typing import Optional, Dict, Tuple
import warnings


def normalize_edge_weights(
    network: ig.Graph,
    edge_attr: str,
    method: Optional[str] = 'minmax',
    output_attr: Optional[str] = None
) -> ig.Graph:
    """
    Normalize edge weights to [0, 1] range.
    
    Parameters
    ----------
    network : ig.Graph
        Network with edge weights
    edge_attr : str
        Name of edge attribute to normalize
    method : str, optional
        Normalization method: 'minmax', 'log1p', 'power', or None (default: 'minmax').
        If None, copies weights without normalization (assumes already normalized).
        'power': Normalizes to [0, 1] then raises to power 6 to emphasize stronger edges
    output_attr : str, optional
        Name for normalized attribute. If None, uses f'{edge_attr}_norm'
        
    Returns
    -------
    ig.Graph
        Network with normalized edge weights (modifies in place and returns)
    """
    if output_attr is None:
        output_attr = f'{edge_attr}_norm'
    
    # Extract all edge weights - batch operation
    try:
        weights_raw = network.es[edge_attr]
        # Collect edges with weights and their indices
        edges_with_weights = [(i, w) for i, w in enumerate(weights_raw) if w is not None]
        if not edges_with_weights:
            warnings.warn(f"No edges have attribute '{edge_attr}'. Skipping normalization.")
            return network
        
        # Extract weights for normalization
        weights = np.array([w for _, w in edges_with_weights])
    except (KeyError, TypeError):
        warnings.warn(f"No edges have attribute '{edge_attr}'. Skipping normalization.")
        return network
    
    # If method is None, copy weights without normalization
    if method is None:
        normalized = weights
    # Normalize based on method
    elif method == 'minmax':
        min_w = np.min(weights)
        max_w = np.max(weights)
        
        if max_w > min_w:
            normalized = (weights - min_w) / (max_w - min_w)
        else:
            # All weights are the same
            normalized = np.full_like(weights, 0.5)
            warnings.warn("All edge weights are identical. Setting normalized values to 0.5")
    
    elif method == 'log1p':
        # Apply log1p transformation, then normalize to [0, 1]
        log1p_weights = np.log1p(weights)
        min_log = np.min(log1p_weights)
        max_log = np.max(log1p_weights)
        
        if max_log > min_log:
            normalized = (log1p_weights - min_log) / (max_log - min_log)
        else:
            # All log1p-transformed weights are the same
            normalized = np.full_like(log1p_weights, 0.5)
            warnings.warn("All log1p-transformed weights are identical. Setting normalized values to 0.5")
    
    elif method == 'power':
        # Power normalization: raise normalized weights to power 6
        # This emphasizes stronger edges while suppressing weaker ones
        min_w = np.min(weights)
        max_w = np.max(weights)
        
        # Normalize to [0, 1] first
        # normalized_base = (weights - min_w) / (max_w - min_w)
        # Apply power transformation (exponent = 6)
        normalized = np.power(weights, 6)

    else:
        raise ValueError(f"Unknown normalization method: {method}. "
                        "Use 'minmax', 'log1p', 'power', or None.")
    
    # Set normalized weights - only for edges that had weights
    # Initialize output with None
    normalized_list = [None] * network.ecount()
    for idx, (edge_idx, _) in enumerate(edges_with_weights):
        normalized_list[edge_idx] = float(normalized[idx])
    
    network.es[output_attr] = normalized_list
    
    return network


def add_edge_weights_to_network(
    network: ig.Graph,
    edge_weights: Dict[Tuple[str, str], float],
    attr_name: str = 'weight',
    symmetric: bool = True
) -> ig.Graph:
    """
    Add edge weights from a dictionary to the network.
    
    Parameters
    ----------
    network : ig.Graph
        Network to modify
    edge_weights : Dict[Tuple[str, str], float]
        Dictionary mapping (node1, node2) tuples to weights
    attr_name : str, optional
        Name for edge attribute (default: 'weight')
    symmetric : bool, optional
        If True, also set weight for (node2, node1) (default: True)
        
    Returns
    -------
    ig.Graph
        Network with edge weights added (modifies in place and returns)
    """
    # Create name to index mapping
    name_to_idx = {network.vs[i]['name']: i for i in range(network.vcount())}
    
    # Initialize weights with None
    weights_list = [None] * network.ecount()
    
    for (u, v), weight in edge_weights.items():
        u_idx = name_to_idx.get(u)
        v_idx = name_to_idx.get(v)
        
        if u_idx is not None and v_idx is not None:
            try:
                eid = network.get_eid(u_idx, v_idx, directed=False, error=False)
                if eid != -1:
                    weights_list[eid] = weight
            except:
                pass
    
    # Batch assign weights
    network.es[attr_name] = weights_list
    
    return network

