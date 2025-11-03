"""Edge weight computation and normalization utilities."""

import numpy as np
import networkx as nx
from typing import Optional, Dict, Tuple
from scipy.stats import pearsonr
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
        Normalization method: 'minmax', 'zscore', or 'rank' (default: 'minmax')
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
    
    elif method == 'zscore':
        mean_w = np.mean(weights)
        std_w = np.std(weights)
        
        if std_w > 0:
            normalized = (weights - mean_w) / std_w
            # Clip to reasonable range and rescale to [0, 1]
            normalized = np.clip(normalized, -3, 3)
            normalized = (normalized + 3) / 6
        else:
            normalized = np.full_like(weights, 0.5)
            warnings.warn("Zero standard deviation. Setting normalized values to 0.5")
    
    elif method == 'rank':
        # Rank-based normalization
        ranks = np.argsort(np.argsort(weights))
        normalized = ranks / (len(weights) - 1) if len(weights) > 1 else np.array([0.5])
    
    else:
        raise ValueError(f"Unknown normalization method: {method}. "
                        "Use 'minmax', 'zscore', or 'rank'.")
    
    # Set normalized weights
    for (u, v), norm_weight in zip(edges_with_weights, normalized):
        network[u][v][output_attr] = float(norm_weight)
    
    return network


def compute_edge_weights_from_expression(
    network: nx.Graph,
    expression_matrix: np.ndarray,
    gene_names: np.ndarray,
    use_abs: bool = True,
    min_correlation: float = 0.0,
    output_attr: str = 'pearson_corr'
) -> nx.Graph:
    """
    Compute Pearson correlation between genes as edge weights.
    
    Only computes correlations for edges that exist in the network (efficient).
    
    Parameters
    ----------
    network : nx.Graph
        Protein-protein interaction network
    expression_matrix : np.ndarray
        Expression matrix (cells x genes)
    gene_names : np.ndarray
        Array of gene names corresponding to columns in expression_matrix
    use_abs : bool, optional
        Use absolute correlation (default: True)
    min_correlation : float, optional
        Minimum correlation to set (default: 0.0)
    output_attr : str, optional
        Name for edge attribute (default: 'pearson_corr')
        
    Returns
    -------
    nx.Graph
        Network with correlation edge weights (modifies in place and returns)
    """
    # Create gene name to index mapping
    gene_to_idx = {gene: idx for idx, gene in enumerate(gene_names)}
    
    # Compute correlations only for edges in the network
    edges_computed = 0
    edges_skipped = 0
    
    for u, v in network.edges():
        u_str = str(u)
        v_str = str(v)
        
        # Check if both genes are in expression data
        if u_str in gene_to_idx and v_str in gene_to_idx:
            idx_u = gene_to_idx[u_str]
            idx_v = gene_to_idx[v_str]
            
            # Get expression vectors
            expr_u = expression_matrix[:, idx_u]
            expr_v = expression_matrix[:, idx_v]
            
            # Check for valid data
            if np.std(expr_u) > 0 and np.std(expr_v) > 0:
                # Compute Pearson correlation
                corr, _ = pearsonr(expr_u, expr_v)
                
                if use_abs:
                    corr = abs(corr)
                
                # Apply minimum threshold
                corr = max(corr, min_correlation)
                
                network[u][v][output_attr] = float(corr)
                edges_computed += 1
            else:
                # One or both genes have zero variance
                network[u][v][output_attr] = min_correlation
                edges_skipped += 1
        else:
            # One or both genes not in expression data
            edges_skipped += 1
    
    if edges_computed == 0:
        warnings.warn("No edge correlations could be computed. "
                     "Check that gene names match between network and expression data.")
    
    return network


def compute_edge_weights_from_scanpy(
    network: nx.Graph,
    adata,  # AnnData type
    cluster_key: Optional[str] = None,
    cluster_id: Optional[str] = None,
    groupby: Optional[str] = None,
    group: Optional[str] = None,
    layer: Optional[str] = None,
    use_abs: bool = True,
    min_correlation: float = 0.0,
    output_attr: str = 'pearson_corr'
) -> nx.Graph:
    """
    Compute edge weights from scanpy AnnData object.
    
    Parameters
    ----------
    network : nx.Graph
        Protein-protein interaction network
    adata : AnnData
        Scanpy AnnData object with expression data
    cluster_key : str, optional
        DEPRECATED: Use groupby instead. Key in adata.obs for grouping labels.
    cluster_id : str, optional
        DEPRECATED: Use group instead. Specific group to use.
    groupby : str, optional
        Key in adata.obs for grouping labels. If provided, compute correlation
        within specific group only.
    group : str, optional
        Specific group to use (requires groupby)
    layer : str, optional
        Which expression layer to use (None = .X)
    use_abs : bool, optional
        Use absolute correlation (default: True)
    min_correlation : float, optional
        Minimum correlation to set (default: 0.0)
    output_attr : str, optional
        Name for edge attribute (default: 'pearson_corr')
        
    Returns
    -------
    nx.Graph
        Network with correlation edge weights
    """
    # Handle deprecated parameters
    if cluster_key is not None or cluster_id is not None:
        import warnings
        warnings.warn(
            "cluster_key and cluster_id are deprecated. Use groupby and group instead.",
            DeprecationWarning,
            stacklevel=2
        )
        if groupby is None:
            groupby = cluster_key
        if group is None:
            group = cluster_id
    
    # Get expression matrix
    if layer is None:
        expr_matrix = adata.X
    else:
        expr_matrix = adata.layers[layer]
    
    # Convert to dense if sparse
    if hasattr(expr_matrix, 'toarray'):
        expr_matrix = expr_matrix.toarray()
    
    # Filter to specific group if requested
    if groupby is not None and group is not None:
        if groupby not in adata.obs:
            raise ValueError(f"Grouping key '{groupby}' not found in adata.obs")
        
        mask = adata.obs[groupby] == group
        expr_matrix = expr_matrix[mask, :]
        
        if expr_matrix.shape[0] == 0:
            raise ValueError(f"No cells found for group '{group}'")
    
    # Get gene names
    gene_names = adata.var_names.values
    
    # Compute correlations
    return compute_edge_weights_from_expression(
        network=network,
        expression_matrix=expr_matrix,
        gene_names=gene_names,
        use_abs=use_abs,
        min_correlation=min_correlation,
        output_attr=output_attr
    )


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

