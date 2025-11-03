"""Tests for node scoring functions."""

import pytest
import numpy as np
import networkx as nx
from scppin.core.scoring import (
    node_score_function,
    compute_node_scores,
    shift_scores_for_pcst,
    get_minimum_score
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
        network = nx.Graph()
        network.add_nodes_from(['A', 'B', 'C'])
        network.add_edges_from([('A', 'B'), ('B', 'C')])
        
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
        """Test handling of missing data."""
        network = nx.Graph()
        network.add_nodes_from(['A', 'B', 'C', 'D'])
        
        pvalues = {'A': 0.001, 'B': 0.01}  # C and D missing
        
        # Without missing_data_score
        scores = compute_node_scores(
            network, pvalues,
            lambda_param=0.5, alpha=0.5, fdr=0.01,
            missing_data_score=False
        )
        assert len(scores) == 2
        assert 'C' not in scores
        assert 'D' not in scores
        
        # With missing_data_score
        scores = compute_node_scores(
            network, pvalues,
            lambda_param=0.5, alpha=0.5, fdr=0.01,
            missing_data_score=True,
            missing_penalty=-1.0
        )
        assert len(scores) == 4
        assert scores['C'] == -1.0
        assert scores['D'] == -1.0


class TestShiftScores:
    """Test score shifting for PCST."""
    
    def test_shift_scores_basic(self):
        """Test basic score shifting."""
        scores = {'A': 10, 'B': 5, 'C': -5}
        
        shifted = shift_scores_for_pcst(scores)
        
        assert len(shifted) == 3
        # Minimum should be 0
        assert min(shifted.values()) == 0
        # Relative differences preserved
        assert shifted['A'] - shifted['B'] == scores['A'] - scores['B']
    
    def test_shift_scores_all_positive(self):
        """Test shifting when all scores are positive."""
        scores = {'A': 10, 'B': 5, 'C': 1}
        
        shifted = shift_scores_for_pcst(scores)
        
        assert min(shifted.values()) == 0
        assert shifted['C'] == 0
    
    def test_shift_scores_empty(self):
        """Test shifting empty dictionary."""
        shifted = shift_scores_for_pcst({})
        assert shifted == {}


class TestGetMinimumScore:
    """Test minimum score extraction."""
    
    def test_get_minimum_score_basic(self):
        """Test getting minimum score."""
        scores = {'A': 10, 'B': 5, 'C': -5}
        
        min_score = get_minimum_score(scores)
        
        assert min_score == -5
    
    def test_get_minimum_score_empty(self):
        """Test getting minimum from empty dictionary."""
        min_score = get_minimum_score({})
        assert min_score == 0.0


if __name__ == '__main__':
    pytest.main([__file__, '-v'])

