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


def _detect_module(
    network: nx.Graph,
    pvalues: Dict[str, float],
    fdr: float = 0.01,
    edge_weight_attr: Optional[str] = None,
    edge_weight_scale: float = 1.0,
    c0: Optional[float] = None,
    missing_data_score: bool = False,
    simplify: bool = True,
    validate: bool = True
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
    edge_weight_scale : float, optional
        Scaling factor for edge weight influence (0-1 recommended).
        Higher values give more weight to edge confidence. (default: 1.0)
    c0 : float, optional
        Minimum edge cost to prevent zeros. If None, uses 0.1 * base_cost.
        This parameter prevents edges from having zero cost, which can cause
        numerical issues. (default: None)
    missing_data_score : bool, optional
        If True, include genes without expression data in the analysis
        (default: False)
    simplify : bool, optional
        Remove self-loops and parallel edges (default: True)
    validate : bool, optional
        Validate network structure (default: True)
        
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
    >>> # Build network
    >>> network = scppin.build_graph('edges.csv')
    >>> 
    >>> # Get p-values
    >>> pvalues = {'GENE1': 0.001, 'GENE2': 0.05, 'GENE3': 0.0001}
    >>> 
    >>> # Detect module
    >>> module = scppin.detect_module(network, pvalues, fdr=0.01)
    >>> 
    >>> # With edge weights
    >>> network = scppin.build_graph('edges.csv', weights='confidence')
    >>> module = scppin.detect_module(
    ...     network, pvalues, fdr=0.01,
    ...     edge_weight_attr='confidence',
    ...     edge_weight_scale=0.5,
    ...     c0=0.1
    ... )
    
    >>> # Complete workflow with scanpy
    >>> import scanpy as sc
    >>> 
    >>> # Extract p-values
    >>> pvalues = scppin.extract_pvalues(adata, 'louvain', '0')
    >>> 
    >>> # Compute edge weights
    >>> weights = scppin.compute_edge_weights('edges.csv', adata)
    >>> network = scppin.build_graph('edges.csv', weights=weights)
    >>> 
    >>> # Detect module
    >>> module = scppin.detect_module(network, pvalues, fdr=0.01)
    """
    # Validate p-values
    pvalues = dict(pvalues)  # Make a copy
    
    for gene, pval in pvalues.items():
        if pval <= 0 or pval > 1:
            raise ValueError(f"P-value for {gene} is {pval}, must be in (0, 1]")
        if np.isnan(pval) or np.isinf(pval):
            raise ValueError(f"P-value for {gene} is NaN or inf")
    
    if len(pvalues) == 0:
        raise ValueError("pvalues dictionary is empty")
    
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
    pvalues_array = np.array(list(pvalues.values()))
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
        edge_weight_scale=edge_weight_scale,
        c0=c0
    )
    
    return module

