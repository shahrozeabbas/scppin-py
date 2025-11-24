"""Node scoring functions for converting p-values to network weights."""

import numpy as np
import networkx as nx
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
    network: nx.Graph,
    pvalues: Dict[str, float],
    lambda_param: float,
    alpha: float,
    fdr: float,
    missing_data_score: bool = False,
    missing_penalty: float = -1.0
) -> Dict[str, float]:
    """
    Compute node scores for all nodes in the network.
    
    Parameters
    ----------
    network : nx.Graph
        Protein-protein interaction network
    pvalues : Dict[str, float]
        Dictionary mapping gene names to p-values
    lambda_param : float
        BUM mixing parameter
    alpha : float
        BUM shape parameter
    fdr : float
        False discovery rate threshold
    missing_data_score : bool, optional
        If True, assign score to nodes without p-values (default: False)
    missing_penalty : float, optional
        Score for nodes without p-values when missing_data_score=True (default: -1.0)
        
    Returns
    -------
    Dict[str, float]
        Dictionary mapping node names to scores
    """
    from .bum_model import compute_tau_threshold
    
    # Compute tau threshold
    tau = compute_tau_threshold(lambda_param, alpha, fdr)
    
    # Vectorized computation: collect all nodes and p-values
    nodes_list = list(network.nodes())
    node_names = [str(node) for node in nodes_list]
    num_nodes = len(node_names)
    
    # Build p-values array (NaN for missing)
    pvalues_array = np.array([pvalues.get(name, np.nan) for name in node_names])
    
    # Identify valid p-values
    valid_mask = ~np.isnan(pvalues_array) & ~np.isinf(pvalues_array)
    has_pvalue_mask = np.array([name in pvalues for name in node_names])
    
    # Initialize scores array
    if missing_data_score:
        scores_array = np.full(num_nodes, missing_penalty, dtype=float)
    else:
        scores_array = np.full(num_nodes, np.nan, dtype=float)
    
    # Compute scores for valid p-values (vectorized)
    if np.any(valid_mask):
        scores_array[valid_mask] = node_score_function(
            pvalues_array[valid_mask], alpha, tau
        )
    
    # Handle NaN/inf p-values when missing_data_score=True
    if missing_data_score:
        nan_inf_mask = has_pvalue_mask & ~valid_mask
        if np.any(nan_inf_mask):
            scores_array[nan_inf_mask] = missing_penalty
    
    # Build result dictionary (skip NaN scores if missing_data_score=False)
    node_scores = {}
    for name, score in zip(node_names, scores_array):
        if missing_data_score or not np.isnan(score):
            node_scores[name] = float(score)
    
    return node_scores


def add_node_scores_to_network(
    network: nx.Graph,
    node_scores: Dict[str, float],
    attr_name: str = 'score'
) -> nx.Graph:
    """
    Add node scores as node attributes to the network.
    
    Parameters
    ----------
    network : nx.Graph
        Network to modify
    node_scores : Dict[str, float]
        Dictionary mapping node names to scores
    attr_name : str, optional
        Name of the node attribute (default: 'score')
        
    Returns
    -------
    nx.Graph
        Network with scores added (modifies in place and returns)
    """
    nx.set_node_attributes(network, node_scores, attr_name)
    return network


def get_minimum_score(node_scores: Dict[str, float]) -> float:
    """
    Get the minimum score (used for edge cost calculation).
    
    Parameters
    ----------
    node_scores : Dict[str, float]
        Node scores
        
    Returns
    -------
    float
        Minimum score value
    """
    if not node_scores:
        return 0.0
    
    return min(node_scores.values())

