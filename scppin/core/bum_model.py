"""Beta-Uniform Mixture (BUM) model for p-value fitting.

This module implements the BUM model to distinguish true signals from noise
in p-value distributions from differential expression analysis.
"""

import numpy as np
from scipy.optimize import minimize
from typing import Tuple, Optional
import warnings


def bum_density(p: np.ndarray, lambda_param: float, alpha: float) -> np.ndarray:
    """
    Beta-Uniform Mixture density function.
    
    f(p) = λ + (1-λ)αp^(α-1)
    
    Parameters
    ----------
    p : np.ndarray
        P-values in (0, 1]
    lambda_param : float
        Mixing parameter in (0, 1)
    alpha : float
        Beta distribution shape parameter in (0, 1)
        
    Returns
    -------
    np.ndarray
        Density values
    """
    return lambda_param + (1 - lambda_param) * alpha * np.power(p, alpha - 1)


def bum_density_cumulative(p: np.ndarray, lambda_param: float, alpha: float) -> np.ndarray:
    """
    Cumulative distribution function for BUM model.
    
    F(p) = λp + (1-λ)p^α
    
    Parameters
    ----------
    p : np.ndarray
        P-values in (0, 1]
    lambda_param : float
        Mixing parameter in (0, 1)
    alpha : float
        Beta distribution shape parameter in (0, 1)
        
    Returns
    -------
    np.ndarray
        Cumulative probability values
    """
    return lambda_param * p + (1 - lambda_param) * np.power(p, alpha)


def _negative_log_likelihood(params: np.ndarray, pvalues: np.ndarray) -> float:
    """
    Negative log-likelihood for BUM model (for minimization).
    
    Parameters
    ----------
    params : np.ndarray
        [lambda_param, alpha]
    pvalues : np.ndarray
        P-values to fit
        
    Returns
    -------
    float
        Negative log-likelihood value
    """
    lambda_param, alpha = params
    
    # Add small epsilon to avoid log(0)
    epsilon = 1e-10
    density = bum_density(pvalues, lambda_param, alpha)
    density = np.maximum(density, epsilon)
    
    return -np.sum(np.log(density))


def fit_bum(
    pvalues: np.ndarray,
    lambda_init: float = 0.5,
    alpha_init: float = 0.5,
    bounds: Optional[Tuple[Tuple[float, float], Tuple[float, float]]] = None
) -> Tuple[float, float, bool]:
    """
    Fit Beta-Uniform Mixture model to p-values using maximum likelihood.
    
    Parameters
    ----------
    pvalues : np.ndarray
        P-values in (0, 1]. Must not contain zeros.
    lambda_init : float, optional
        Initial value for lambda parameter (default: 0.5)
    alpha_init : float, optional
        Initial value for alpha parameter (default: 0.5)
    bounds : tuple, optional
        Bounds for (lambda, alpha). Default: ((1e-9, 1-1e-9), (1e-9, 1-1e-9))
        
    Returns
    -------
    lambda_param : float
        Fitted mixing parameter
    alpha : float
        Fitted shape parameter
    success : bool
        Whether optimization converged
    """
    # Convert to array (validation done upstream in set_node_weights)
    pvalues = np.asarray(pvalues)
    
    if len(pvalues) < 10:
        warnings.warn("Very few p-values (<10) for BUM fitting. Results may be unreliable.")
    
    # Set bounds
    if bounds is None:
        epsilon = 1e-9
        bounds = ((epsilon, 1 - epsilon), (epsilon, 1 - epsilon))
    
    # Initial parameters
    x0 = np.array([lambda_init, alpha_init])
    
    # Optimize
    result = minimize(
        _negative_log_likelihood,
        x0,
        args=(pvalues,),
        method='L-BFGS-B',
        bounds=bounds
    )
    
    lambda_param, alpha = result.x
    
    if not result.success:
        warnings.warn(f"BUM fitting did not converge: {result.message}")
    
    return lambda_param, alpha, result.success


def compute_tau_threshold(
    lambda_param: float,
    alpha: float,
    fdr: float
) -> float:
    """
    Compute tau threshold from BUM parameters and FDR.
    
    The threshold tau is computed such that the expected FDR is controlled
    at the specified level.
    
    Parameters
    ----------
    lambda_param : float
        BUM mixing parameter
    alpha : float
        BUM shape parameter
    fdr : float
        False discovery rate threshold
        
    Returns
    -------
    float
        Tau threshold value
    """
    # Compute pi_hat (density at p=1)
    pi_hat = lambda_param + (1 - lambda_param) * alpha
    
    # Compute tau from FDR
    # tau = ((pi_hat - fdr*lambda) / (fdr*(1-lambda)))^(1/(alpha-1))
    numerator = pi_hat - (fdr * lambda_param)
    denominator = fdr * (1 - lambda_param)
    
    if denominator <= 0 or numerator <= 0:
        warnings.warn("Invalid parameters for tau computation. Using FDR as threshold.")
        return fdr
    
    tau = np.power(numerator / denominator, 1.0 / (alpha - 1))
    
    # Ensure tau is in valid range
    tau = np.clip(tau, 0, 1)
    
    return tau

