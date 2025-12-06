"""Node scoring functions for converting p-values to network weights."""

import numpy as np
import igraph as ig
from typing import Dict, Optional
import warnings


def node_score_function(
    pvalues: np.ndarray,
    alpha: float,
    tau: float
) -> np.ndarray:
    """
    Compute node scores from p-values using BUM parameters.
    
    score = (α - 1) * (log(p) - log(τ))
    
    Parameters
    ----------
    pvalues : np.ndarray
        P-values for genes
    alpha : float
        BUM shape parameter
    tau : float
        Threshold from FDR
        
    Returns
    -------
    np.ndarray
        Node scores (higher = more significant)
    """
    # Avoid log(0) by adding small epsilon
    epsilon = 1e-300
    pvalues_safe = np.maximum(pvalues, epsilon)
    tau_safe = max(tau, epsilon)
    
    scores = (alpha - 1) * (np.log(pvalues_safe) - np.log(tau_safe))
    
    return scores


def compute_node_scores(
    network: ig.Graph,
    pvalues: Dict[str, float],
    lambda_param: float,
    alpha: float,
    fdr: float
) -> Dict[str, float]:
    """
    Compute node scores for all nodes in the network.
    
    Parameters
    ----------
    network : ig.Graph
        Protein-protein interaction network (should be filtered to genes with p-values)
    pvalues : Dict[str, float]
        Dictionary mapping gene names to p-values
    lambda_param : float
        BUM mixing parameter
    alpha : float
        BUM shape parameter
    fdr : float
        False discovery rate threshold
        
    Returns
    -------
    Dict[str, float]
        Dictionary mapping node names to scores
    """
    from .bum_model import compute_tau_threshold
    
    # Compute tau threshold
    tau = compute_tau_threshold(lambda_param, alpha, fdr)
    
    # Vectorized computation: collect all nodes and p-values
    node_names = network.vs['name']
    
    # Build p-values array (NaN for missing)
    pvalues_array = np.array([pvalues.get(name, np.nan) for name in node_names])
    
    # Identify valid p-values
    valid_mask = ~np.isnan(pvalues_array) & ~np.isinf(pvalues_array)
    
    # Initialize scores array with NaN
    scores_array = np.full(len(node_names), np.nan, dtype=float)
    
    # Compute scores for valid p-values (vectorized)
    if np.any(valid_mask):
        scores_array[valid_mask] = node_score_function(
            pvalues_array[valid_mask], alpha, tau
        )
    
    # Build result dictionary (skip NaN scores)
    node_scores = {}
    for name, score in zip(node_names, scores_array):
        if not np.isnan(score):
            node_scores[name] = float(score)
    
    return node_scores
