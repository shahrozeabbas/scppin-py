"""Tests for Beta-Uniform Mixture model."""

import pytest
import numpy as np
from scppin.core.bum_model import (
    bum_density,
    bum_density_cumulative,
    fit_bum,
    compute_tau_threshold
)


class TestBUMDensity:
    """Test BUM density functions."""
    
    def test_bum_density_basic(self):
        """Test basic BUM density calculation."""
        p = np.array([0.1, 0.5, 0.9])
        lambda_param = 0.5
        alpha = 0.5
        
        density = bum_density(p, lambda_param, alpha)
        
        assert len(density) == 3
        assert np.all(density > 0)
        assert np.all(np.isfinite(density))
    
    def test_bum_density_edge_cases(self):
        """Test BUM density at edge cases."""
        # At p=1
        density = bum_density(np.array([1.0]), 0.5, 0.5)
        assert np.isfinite(density[0])
        
        # Very small p
        density = bum_density(np.array([0.001]), 0.5, 0.5)
        assert np.isfinite(density[0])
    
    def test_bum_density_cumulative(self):
        """Test cumulative BUM density."""
        p = np.array([0.1, 0.5, 0.9])
        lambda_param = 0.5
        alpha = 0.5
        
        cdf = bum_density_cumulative(p, lambda_param, alpha)
        
        assert len(cdf) == 3
        assert np.all(cdf >= 0)
        assert np.all(cdf <= 1)
        # CDF should be monotonically increasing
        assert cdf[0] < cdf[1] < cdf[2]


class TestFitBUM:
    """Test BUM model fitting."""
    
    def test_fit_bum_uniform(self):
        """Test fitting to uniform distribution."""
        # Generate uniform p-values
        np.random.seed(42)
        pvalues = np.random.uniform(0.01, 1.0, 1000)
        
        lambda_param, alpha, success = fit_bum(pvalues)
        
        assert success
        assert 0 < lambda_param < 1
        assert 0 < alpha < 1
        # For uniform, lambda should be close to 1
        assert lambda_param > 0.8
    
    def test_fit_bum_beta(self):
        """Test fitting to beta-like distribution."""
        # Generate beta-distributed p-values
        np.random.seed(42)
        pvalues = np.random.beta(0.5, 5, 1000)
        pvalues = np.clip(pvalues, 0.001, 1.0)  # Ensure in valid range
        
        lambda_param, alpha, success = fit_bum(pvalues)
        
        assert success
        assert 0 < lambda_param < 1
        assert 0 < alpha < 1
    
    def test_fit_bum_invalid_pvalues(self):
        """Test that invalid p-values raise errors."""
        # P-values with zeros
        with pytest.raises(ValueError, match="must be in the interval"):
            fit_bum(np.array([0.0, 0.5, 1.0]))
        
        # P-values > 1
        with pytest.raises(ValueError, match="must be in the interval"):
            fit_bum(np.array([0.5, 1.5]))
        
        # P-values with NaN
        with pytest.raises(ValueError, match="NaN or inf"):
            fit_bum(np.array([0.5, np.nan, 0.8]))
    
    def test_fit_bum_few_values(self):
        """Test fitting with few p-values (should warn)."""
        pvalues = np.array([0.1, 0.2, 0.3])
        
        with pytest.warns(UserWarning, match="Very few p-values"):
            lambda_param, alpha, success = fit_bum(pvalues)
        
        assert 0 < lambda_param < 1
        assert 0 < alpha < 1


class TestComputeTau:
    """Test tau threshold computation."""
    
    def test_compute_tau_basic(self):
        """Test basic tau computation."""
        lambda_param = 0.5
        alpha = 0.5
        fdr = 0.01
        
        tau = compute_tau_threshold(lambda_param, alpha, fdr)
        
        assert 0 <= tau <= 1
        assert np.isfinite(tau)
    
    def test_compute_tau_different_fdr(self):
        """Test that tau decreases with stricter FDR."""
        lambda_param = 0.5
        alpha = 0.5
        
        tau_001 = compute_tau_threshold(lambda_param, alpha, 0.01)
        tau_005 = compute_tau_threshold(lambda_param, alpha, 0.05)
        
        # Stricter FDR should give smaller tau
        assert tau_001 < tau_005
    
    def test_compute_tau_edge_cases(self):
        """Test tau computation at edge cases."""
        # Very small FDR
        tau = compute_tau_threshold(0.5, 0.5, 0.001)
        assert 0 <= tau <= 1
        
        # Large FDR
        tau = compute_tau_threshold(0.5, 0.5, 0.1)
        assert 0 <= tau <= 1


if __name__ == '__main__':
    pytest.main([__file__, '-v'])

