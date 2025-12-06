"""Tests for node scoring functions."""

import pytest
import numpy as np
import igraph as ig
from scppin.core.scoring import (
    node_score_function,
    compute_node_scores,
)


class TestNodeScoreFunction:
    """Test node score computation."""
    
    def test_node_score_basic(self):
        """Test basic node score calculation."""
        pvalues = np.array([0.001, 0.01, 0.1])
        alpha = 0.5
        tau = 0.05
        
        scores = node_score_function(pvalues, alpha, tau)
        
        assert len(scores) == 3
        assert np.all(np.isfinite(scores))
        # More significant p-values should have higher scores
        assert scores[0] > scores[1] > scores[2]
    
    def test_node_score_zero_pvalue(self):
        """Test that very small p-values are handled."""
        pvalues = np.array([1e-300, 0.01])
        alpha = 0.5
        tau = 0.05
        
        scores = node_score_function(pvalues, alpha, tau)
        
        assert np.all(np.isfinite(scores))


class TestComputeNodeScores:
    """Test node score computation for networks."""
    
    def test_compute_node_scores_basic(self):
        """Test computing scores for a simple network."""
        # Create simple network
        network = ig.Graph.TupleList([('A', 'B'), ('B', 'C')], directed=False, vertex_name_attr='name')
        
        pvalues = {'A': 0.001, 'B': 0.01, 'C': 0.1}
        
        scores = compute_node_scores(
            network, pvalues,
            lambda_param=0.5, alpha=0.5, fdr=0.01
        )
        
        assert len(scores) == 3
        assert 'A' in scores
        assert 'B' in scores
        assert 'C' in scores
        assert scores['A'] > scores['B'] > scores['C']
    
    def test_compute_node_scores_missing_data(self):
        """Test handling of nodes without p-values."""
        network = ig.Graph.TupleList([('A', 'B'), ('B', 'C'), ('C', 'D')], directed=False, vertex_name_attr='name')
        
        pvalues = {'A': 0.001, 'B': 0.01}  # C and D missing
        
        # Only nodes with p-values get scores
        scores = compute_node_scores(
            network, pvalues,
            lambda_param=0.5, alpha=0.5, fdr=0.01
        )
        assert len(scores) == 2
        assert 'C' not in scores
        assert 'D' not in scores


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
