"""Tests for edge weight functionality."""

import pytest
import numpy as np
import networkx as nx
from scppin.core.edge_weights import (
    normalize_edge_weights,
    compute_edge_weights_from_expression,
    add_edge_weights_to_network
)


class TestNormalizeEdgeWeights:
    """Test edge weight normalization."""
    
    def test_normalize_minmax(self):
        """Test min-max normalization."""
        network = nx.Graph()
        network.add_edge('A', 'B', weight=0.1)
        network.add_edge('B', 'C', weight=0.5)
        network.add_edge('C', 'D', weight=0.9)
        
        network = normalize_edge_weights(network, 'weight', method='minmax')
        
        # Check normalized values
        assert network['A']['B']['weight_norm'] == 0.0
        assert network['C']['D']['weight_norm'] == 1.0
        assert 0 < network['B']['C']['weight_norm'] < 1
    
    def test_normalize_all_same(self):
        """Test normalization when all weights are the same."""
        network = nx.Graph()
        network.add_edge('A', 'B', weight=0.5)
        network.add_edge('B', 'C', weight=0.5)
        
        with pytest.warns(UserWarning, match="identical"):
            network = normalize_edge_weights(network, 'weight', method='minmax')
        
        # Should set to 0.5
        assert network['A']['B']['weight_norm'] == 0.5
        assert network['B']['C']['weight_norm'] == 0.5
    
    def test_normalize_missing_attribute(self):
        """Test normalization with missing attribute."""
        network = nx.Graph()
        network.add_edge('A', 'B')  # No weight attribute
        
        with pytest.warns(UserWarning, match="No edges have attribute"):
            network = normalize_edge_weights(network, 'weight', method='minmax')
    
    def test_normalize_zscore(self):
        """Test z-score normalization."""
        network = nx.Graph()
        network.add_edge('A', 'B', weight=0.1)
        network.add_edge('B', 'C', weight=0.5)
        network.add_edge('C', 'D', weight=0.9)
        
        network = normalize_edge_weights(network, 'weight', method='zscore')
        
        # All values should be in [0, 1] after z-score normalization
        for u, v in network.edges():
            assert 0 <= network[u][v]['weight_norm'] <= 1


class TestComputeEdgeWeightsFromExpression:
    """Test computing edge weights from expression data."""
    
    def test_compute_correlation_basic(self):
        """Test basic correlation computation."""
        # Create network
        network = nx.Graph()
        network.add_edge('GENE1', 'GENE2')
        network.add_edge('GENE2', 'GENE3')
        
        # Create expression matrix (cells x genes)
        np.random.seed(42)
        expr_matrix = np.random.randn(100, 3)
        gene_names = np.array(['GENE1', 'GENE2', 'GENE3'])
        
        network = compute_edge_weights_from_expression(
            network, expr_matrix, gene_names
        )
        
        # Check that correlations were added
        assert 'pearson_corr' in network['GENE1']['GENE2']
        assert 'pearson_corr' in network['GENE2']['GENE3']
        
        # Correlations should be in [0, 1] (absolute values)
        for u, v in network.edges():
            corr = network[u][v]['pearson_corr']
            assert 0 <= corr <= 1
    
    def test_compute_correlation_perfect(self):
        """Test correlation with perfectly correlated genes."""
        network = nx.Graph()
        network.add_edge('GENE1', 'GENE2')
        
        # Create perfectly correlated expression
        expr_matrix = np.array([[1, 1], [2, 2], [3, 3], [4, 4]])
        gene_names = np.array(['GENE1', 'GENE2'])
        
        network = compute_edge_weights_from_expression(
            network, expr_matrix, gene_names
        )
        
        # Should be 1.0 (perfect correlation)
        assert np.isclose(network['GENE1']['GENE2']['pearson_corr'], 1.0)
    
    def test_compute_correlation_missing_genes(self):
        """Test correlation when some genes are missing."""
        network = nx.Graph()
        network.add_edge('GENE1', 'GENE2')
        network.add_edge('GENE2', 'GENE3')
        network.add_edge('GENE3', 'GENE4')  # GENE4 not in expression
        
        expr_matrix = np.random.randn(100, 3)
        gene_names = np.array(['GENE1', 'GENE2', 'GENE3'])
        
        with pytest.warns(UserWarning, match="No edge correlations"):
            network = compute_edge_weights_from_expression(
                network, expr_matrix, gene_names
            )


class TestAddEdgeWeights:
    """Test adding edge weights from dictionary."""
    
    def test_add_edge_weights_basic(self):
        """Test adding edge weights."""
        network = nx.Graph()
        network.add_edge('A', 'B')
        network.add_edge('B', 'C')
        
        edge_weights = {
            ('A', 'B'): 0.8,
            ('B', 'C'): 0.6
        }
        
        network = add_edge_weights_to_network(network, edge_weights)
        
        assert network['A']['B']['weight'] == 0.8
        assert network['B']['C']['weight'] == 0.6
    
    def test_add_edge_weights_symmetric(self):
        """Test symmetric edge weight addition."""
        network = nx.Graph()
        network.add_edge('A', 'B')
        
        edge_weights = {('A', 'B'): 0.8}
        
        network = add_edge_weights_to_network(
            network, edge_weights, symmetric=True
        )
        
        # Both directions should have weight
        assert network['A']['B']['weight'] == 0.8


if __name__ == '__main__':
    pytest.main([__file__, '-v'])

