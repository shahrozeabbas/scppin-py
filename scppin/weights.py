"""Edge weight computation from expression data."""

import numpy as np
import pandas as pd
import networkx as nx
from typing import Union, List, Tuple, Optional, Dict
from scipy.stats import pearsonr, spearmanr
from pathlib import Path


def _compute_edge_weights(
    edges: Union[str, List[Tuple], nx.Graph],
    adata,
    cluster_key: Optional[str] = None,
    cluster_id: Optional[str] = None,
    method: str = 'pearson'
) -> Dict[Tuple[str, str], float]:
    """
    Compute edge weights from expression data.
    
    Computes correlation between genes and returns normalized weights in [0, 1].
    
    Implementation:
    - ALWAYS uses absolute value of correlation (positive weights)
    - ALWAYS normalizes to [0, 1] using min-max scaling
    - Only computes for edges in provided edge list
    
    Parameters
    ----------
    edges : str, list of tuples, or nx.Graph
        Edge list to compute weights for. Can be:
        - Path to edge list file
        - List of tuples: [('GENE1', 'GENE2'), ...]
        - NetworkX graph (will extract edges)
    adata : AnnData
        Expression data from scanpy
    cluster_key : str, optional
        Compute correlation within specific cluster
    cluster_id : str, optional
        Specific cluster to use (requires cluster_key)
    method : str, optional
        Correlation method: 'pearson' or 'spearman' (default: 'pearson')
        
    Returns
    -------
    Dict[Tuple[str, str], float]
        Edge weights normalized to [0, 1] range.
        Keys are (gene1, gene2) tuples, values are weights.
        
    Examples
    --------
    >>> # Compute weights for edge list
    >>> weights = scppin.compute_edge_weights('edges.csv', adata)
    >>> network = scppin.build_graph('edges.csv', weights=weights)
    
    >>> # Compute within specific cluster
    >>> weights = scppin.compute_edge_weights(
    ...     'edges.csv', adata,
    ...     cluster_key='louvain',
    ...     cluster_id='0'
    ... )
    
    >>> # From existing network
    >>> network = scppin.build_graph('edges.csv')
    >>> weights = scppin.compute_edge_weights(network, adata)
    >>> # Add weights to network
    >>> for (u, v), w in weights.items():
    ...     if network.has_edge(u, v):
    ...         network[u][v]['weight'] = w
    """
    # Extract edge list
    if isinstance(edges, nx.Graph):
        edge_list = list(edges.edges())
    elif isinstance(edges, (str, Path)):
        # Load from file
        from .graph import _load_edges_from_file
        edges_df = _load_edges_from_file(edges)
        edge_list = [(str(row['source']), str(row['target'])) 
                     for _, row in edges_df.iterrows()]
    elif isinstance(edges, list):
        edge_list = [(str(u), str(v)) for u, v in edges]
    else:
        raise TypeError(f"edges must be str, list, or nx.Graph, got {type(edges)}")
    
    # Get expression matrix
    if cluster_key is not None and cluster_id is not None:
        # Filter to specific cluster
        if cluster_key not in adata.obs:
            raise ValueError(f"Cluster key '{cluster_key}' not found in adata.obs")
        
        mask = adata.obs[cluster_key] == cluster_id
        expr_matrix = adata.X[mask, :]
        
        if expr_matrix.shape[0] == 0:
            raise ValueError(f"No cells found for cluster '{cluster_id}'")
    else:
        expr_matrix = adata.X
    
    # Convert to dense if sparse
    if hasattr(expr_matrix, 'toarray'):
        expr_matrix = expr_matrix.toarray()
    
    # Get gene names
    gene_names = adata.var_names.values
    gene_to_idx = {gene: idx for idx, gene in enumerate(gene_names)}
    
    # Compute correlations for each edge
    correlations = {}
    
    for u, v in edge_list:
        if u in gene_to_idx and v in gene_to_idx:
            idx_u = gene_to_idx[u]
            idx_v = gene_to_idx[v]
            
            # Get expression vectors
            expr_u = expr_matrix[:, idx_u]
            expr_v = expr_matrix[:, idx_v]
            
            # Check for valid data (non-zero variance)
            if np.std(expr_u) > 0 and np.std(expr_v) > 0:
                # Compute correlation
                if method == 'pearson':
                    corr, _ = pearsonr(expr_u, expr_v)
                elif method == 'spearman':
                    corr, _ = spearmanr(expr_u, expr_v)
                else:
                    raise ValueError(f"Unknown method: {method}. Use 'pearson' or 'spearman'.")
                
                # ALWAYS use absolute value
                corr = abs(corr)
                
                correlations[(u, v)] = corr
            else:
                # Zero variance - set to 0
                correlations[(u, v)] = 0.0
    
    if not correlations:
        raise ValueError("No correlations could be computed. "
                        "Check that gene names match between network and expression data.")
    
    # ALWAYS normalize to [0, 1]
    correlations = _normalize_weights(correlations)
    
    return correlations


def _normalize_weights(weights: Dict[Tuple[str, str], float]) -> Dict[Tuple[str, str], float]:
    """Normalize weights to [0, 1] using min-max scaling."""
    if not weights:
        return weights
    
    values = list(weights.values())
    min_val = min(values)
    max_val = max(values)
    
    if max_val > min_val:
        # Min-max normalization
        normalized = {
            edge: (weight - min_val) / (max_val - min_val)
            for edge, weight in weights.items()
        }
    else:
        # All weights are the same - set to 0.5
        normalized = {edge: 0.5 for edge in weights.keys()}
    
    return normalized

