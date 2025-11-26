"""Edge weight computation and normalization utilities."""

import numpy as np
import networkx as nx
from typing import Optional, Dict, Tuple
import warnings


def normalize_edge_weights(
    network: nx.Graph,
    edge_attr: str,
    method: str = 'minmax',
    output_attr: Optional[str] = None
) -> nx.Graph:
    """
    Normalize edge weights to [0, 1] range.
    
    Parameters
    ----------
    network : nx.Graph
        Network with edge weights
    edge_attr : str
        Name of edge attribute to normalize
    method : str, optional
        Normalization method: 'minmax', 'log1p', or 'power' (default: 'minmax').
        'power': Normalizes to [0, 1] then raises to power 6 to emphasize stronger edges
    output_attr : str, optional
        Name for normalized attribute. If None, uses f'{edge_attr}_norm'
        
    Returns
    -------
    nx.Graph
        Network with normalized edge weights (modifies in place and returns)
    """
    if output_attr is None:
        output_attr = f'{edge_attr}_norm'
    
    # Extract all edge weights
    weights = []
    edges_with_weights = []
    
    for u, v in network.edges():
        if edge_attr in network[u][v]:
            weight = network[u][v][edge_attr]
            weights.append(weight)
            edges_with_weights.append((u, v))
    
    if not weights:
        warnings.warn(f"No edges have attribute '{edge_attr}'. Skipping normalization.")
        return network
    
    weights = np.array(weights)
    
    # Normalize based on method
    if method == 'minmax':
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
                        "Use 'minmax', 'log1p', or 'power'.")
    
    # Set normalized weights
    for (u, v), norm_weight in zip(edges_with_weights, normalized):
        network[u][v][output_attr] = float(norm_weight)
    
    return network


def add_edge_weights_to_network(
    network: nx.Graph,
    edge_weights: Dict[Tuple[str, str], float],
    attr_name: str = 'weight',
    symmetric: bool = True
) -> nx.Graph:
    """
    Add edge weights from a dictionary to the network.
    
    Parameters
    ----------
    network : nx.Graph
        Network to modify
    edge_weights : Dict[Tuple[str, str], float]
        Dictionary mapping (node1, node2) tuples to weights
    attr_name : str, optional
        Name for edge attribute (default: 'weight')
    symmetric : bool, optional
        If True, also set weight for (node2, node1) (default: True)
        
    Returns
    -------
    nx.Graph
        Network with edge weights added (modifies in place and returns)
    """
    for (u, v), weight in edge_weights.items():
        if network.has_edge(u, v):
            network[u][v][attr_name] = weight
        
        if symmetric and network.has_edge(v, u):
            network[v][u][attr_name] = weight
    
    return network

