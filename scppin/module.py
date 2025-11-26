"""Core module detection function."""

import networkx as nx
import numpy as np
from typing import Dict, Optional
import warnings

from .core import (
    fit_bum,
    compute_node_scores,
    filter_network_by_pvalues,
    simplify_network,
    validate_network
)
from .core.pcst_solver import detect_functional_module_core


def _validate_pvalues(pvalues: Dict[str, float]) -> np.ndarray:
    """
    Validate p-values and convert to numpy array.
    
    Parameters
    ----------
    pvalues : Dict[str, float]
        Dictionary mapping gene names to p-values
        
    Returns
    -------
    np.ndarray
        Validated p-values as numpy array
        
    Raises
    ------
    ValueError
        If p-values are empty, invalid range, or contain NaN/inf
    """
    # Early empty check (fail fast)
    if not pvalues:
        raise ValueError("pvalues dictionary is empty")
    
    # Convert to numpy array and validate using vectorized operations
    pvalues_array = np.array(list(pvalues.values()), dtype=float)
    
    # Vectorized validation for range check
    invalid_range = (pvalues_array <= 0) | (pvalues_array > 1)
    if np.any(invalid_range):
        # Find invalid genes only when needed for error message
        invalid_genes = [gene for gene, pval in pvalues.items() 
                        if pval <= 0 or pval > 1]
        invalid_preview = invalid_genes[:5]
        preview_str = ', '.join(f"{g}: {pvalues[g]}" for g in invalid_preview)
        if len(invalid_genes) > 5:
            preview_str += f", ... ({len(invalid_genes)} total)"
        raise ValueError(f"P-values must be in (0, 1]. Invalid: {preview_str}")
    
    # Vectorized validation for NaN/inf check
    invalid_nan = np.isnan(pvalues_array) | np.isinf(pvalues_array)
    if np.any(invalid_nan):
        # Find invalid genes only when needed for error message
        invalid_genes = [gene for gene, pval in pvalues.items() 
                        if np.isnan(pval) or np.isinf(pval)]
        invalid_preview = invalid_genes[:5]
        preview_str = ', '.join(f"{g}: {pvalues[g]}" for g in invalid_preview)
        if len(invalid_genes) > 5:
            preview_str += f", ... ({len(invalid_genes)} total)"
        raise ValueError(f"P-values contain NaN or inf. Invalid: {preview_str}")
    
    return pvalues_array


def _detect_module(
    network: nx.Graph,
    pvalues: Dict[str, float],
    fdr: float = 0.01,
    edge_weight_attr: Optional[str] = None,
    c0: float = 0.01,
    normalization: str = 'minmax',
    missing_data_score: bool = False,
    simplify: bool = True,
    validate: bool = True,
    num_clusters: int = 1,
    pruning: str = 'gw',
    use_max_prize_root: bool = False
) -> nx.Graph:
    """
    Detect functional module in protein-protein interaction network.
    
    This is the core function for detecting functional modules by integrating
    p-values from differential expression analysis with a protein-protein
    interaction network.
    
    Parameters
    ----------
    network : nx.Graph
        Protein-protein interaction network (PPIN)
    pvalues : Dict[str, float]
        Dictionary mapping gene names to p-values from differential expression.
        P-values must be in (0, 1] (zeros not allowed).
    fdr : float, optional
        False discovery rate threshold (default: 0.01)
    edge_weight_attr : str, optional
        Name of edge attribute containing weights (e.g., 'weight', 'confidence').
        If None, uses uniform edge costs matching R implementation behavior.
        Weights should be in [0, 1]. Higher weights indicate stronger preference
        for including the edge. (default: None)
    c0 : float, optional
        Minimum edge cost to prevent zeros. If None, uses 0.1 * base_cost.
        This parameter prevents edges from having zero cost, which can cause
        numerical issues. (default: None)
    normalization : str, optional
        Normalization method for edge weights: 'minmax', 'log1p', or 'power'
        (default: 'minmax'). Only used when edge_weight_attr is provided.
    missing_data_score : bool, optional
        If True, include genes without expression data in the analysis
        (default: False)
    simplify : bool, optional
        Remove self-loops and parallel edges (default: True)
    validate : bool, optional
        Validate network structure (default: True)
    num_clusters : int, optional
        Number of connected components to return (default: 1)
    pruning : str, optional
        Pruning method: 'gw' (Goemans-Williamson) or 'strong' (default: 'gw')
    use_max_prize_root : bool, optional
        If True, use the node with highest prize as root (default: False)
        
    Returns
    -------
    nx.Graph
        Subgraph representing the functional module. Node attributes include:
        - 'score': Original node score
        - 'prize': Shifted score used in PCST
        
    Raises
    ------
    ValueError
        If p-values are invalid or network is empty
        
    Examples
    --------
    >>> import scppin
    >>> 
    >>> from scppin import scPPIN
    >>> analyzer = scPPIN()
    >>> analyzer.load_network('edges.csv')
    >>> analyzer.set_node_weights({'GENE1': 1e-4, 'GENE2': 5e-3})
    >>> module = analyzer.detect_module(fdr=0.01)
    >>> 
    >>> # With edge weights
    >>> weights = {('GENE1', 'GENE2'): 0.9}
    >>> analyzer.set_edge_weights(weights=weights)
    >>> module = analyzer.detect_module(
    ...     fdr=0.01,
    ...     edge_weight_attr='weight',
    ...     c0=0.1
    ... )
    """
    # Validate p-values and convert to array
    pvalues_array = _validate_pvalues(pvalues)
    
    # Simplify network if requested (in-place for efficiency)
    if simplify:
        network = simplify_network(network, copy=False)
    
    # Validate network
    if validate:
        validate_network(network)
    
    # Filter network to genes with p-values (unless missing_data_score=True)
    network_filtered = filter_network_by_pvalues(
        network, pvalues, missing_data_score
    )
    
    if network_filtered.number_of_nodes() == 0:
        raise ValueError("No genes in network have p-values")
    
    # Fit BUM model
    lambda_param, alpha, success = fit_bum(pvalues_array)
    
    if not success:
        warnings.warn("BUM model fitting did not converge. Results may be unreliable.")
    
    # Compute node scores
    node_scores = compute_node_scores(
        network_filtered,
        pvalues,
        lambda_param,
        alpha,
        fdr,
        missing_data_score=missing_data_score
    )
    
    if not node_scores:
        raise ValueError("No node scores computed. Check that gene names match "
                        "between network and p-values.")
    
    # Detect functional module using PCST
    module = detect_functional_module_core(
        network_filtered,
        node_scores,
        edge_weight_attr=edge_weight_attr,
        c0=c0,
        normalization=normalization,
        num_clusters=num_clusters,
        pruning=pruning,
        use_max_prize_root=use_max_prize_root
    )
    
    return module
