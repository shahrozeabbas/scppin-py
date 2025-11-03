"""Tests for network utilities."""

import pytest
import networkx as nx
import tempfile
import os
from scppin.core.network_utils import (
    simplify_network,
    filter_network,
    filter_network_by_pvalues,
    get_largest_connected_component,
    network_statistics,
    validate_network
)


class TestSimplifyNetwork:
    """Test network simplification."""
    
    def test_remove_self_loops(self):
        """Test removal of self-loops."""
        network = nx.Graph()
        network.add_edge('A', 'A')  # Self-loop
        network.add_edge('A', 'B')
        
        simplified = simplify_network(network)
        
        assert not simplified.has_edge('A', 'A')
        assert simplified.has_edge('A', 'B')
    
    def test_no_changes_needed(self):
        """Test when network is already simple."""
        network = nx.Graph()
        network.add_edge('A', 'B')
        network.add_edge('B', 'C')
        
        simplified = simplify_network(network)
        
        assert simplified.number_of_nodes() == 3
        assert simplified.number_of_edges() == 2


class TestFilterNetwork:
    """Test network filtering."""
    
    def test_filter_keep_genes(self):
        """Test filtering to keep specific genes."""
        network = nx.Graph()
        network.add_nodes_from(['A', 'B', 'C', 'D'])
        network.add_edges_from([('A', 'B'), ('B', 'C'), ('C', 'D')])
        
        genes_to_keep = {'A', 'B', 'C'}
        
        filtered = filter_network(network, genes_to_keep, keep_only_genes=True)
        
        assert filtered.number_of_nodes() == 3
        assert 'D' not in filtered.nodes()
    
    def test_filter_remove_genes(self):
        """Test filtering to remove specific genes."""
        network = nx.Graph()
        network.add_nodes_from(['A', 'B', 'C', 'D'])
        network.add_edges_from([('A', 'B'), ('B', 'C'), ('C', 'D')])
        
        genes_to_remove = {'D'}
        
        filtered = filter_network(network, genes_to_remove, keep_only_genes=False)
        
        assert filtered.number_of_nodes() == 3
        assert 'D' not in filtered.nodes()
    
    def test_filter_empty_result(self):
        """Test filtering that results in empty network."""
        network = nx.Graph()
        network.add_nodes_from(['A', 'B'])
        
        with pytest.warns(UserWarning, match="No nodes remain"):
            filtered = filter_network(network, {'C', 'D'}, keep_only_genes=True)
        
        assert filtered.number_of_nodes() == 0


class TestFilterNetworkByPvalues:
    """Test filtering by p-values."""
    
    def test_filter_with_pvalues(self):
        """Test filtering to genes with p-values."""
        network = nx.Graph()
        network.add_nodes_from(['A', 'B', 'C', 'D'])
        
        pvalues = {'A': 0.01, 'B': 0.05}
        
        filtered = filter_network_by_pvalues(
            network, pvalues, missing_data_score=False
        )
        
        assert filtered.number_of_nodes() == 2
        assert 'A' in filtered.nodes()
        assert 'B' in filtered.nodes()
    
    def test_filter_with_missing_data(self):
        """Test filtering with missing_data_score=True."""
        network = nx.Graph()
        network.add_nodes_from(['A', 'B', 'C', 'D'])
        
        pvalues = {'A': 0.01, 'B': 0.05}
        
        filtered = filter_network_by_pvalues(
            network, pvalues, missing_data_score=True
        )
        
        # Should keep all nodes
        assert filtered.number_of_nodes() == 4


class TestLargestConnectedComponent:
    """Test largest connected component extraction."""
    
    def test_single_component(self):
        """Test when network is already connected."""
        network = nx.Graph()
        network.add_edges_from([('A', 'B'), ('B', 'C')])
        
        lcc = get_largest_connected_component(network)
        
        assert lcc.number_of_nodes() == 3
    
    def test_multiple_components(self):
        """Test extracting largest from multiple components."""
        network = nx.Graph()
        # Component 1: A-B-C
        network.add_edges_from([('A', 'B'), ('B', 'C')])
        # Component 2: D-E
        network.add_edges_from([('D', 'E')])
        
        with pytest.warns(UserWarning, match="connected components"):
            lcc = get_largest_connected_component(network)
        
        assert lcc.number_of_nodes() == 3
        assert 'A' in lcc.nodes()
        assert 'D' not in lcc.nodes()


class TestNetworkStatistics:
    """Test network statistics computation."""
    
    def test_statistics_basic(self):
        """Test basic statistics computation."""
        network = nx.Graph()
        network.add_edges_from([('A', 'B'), ('B', 'C'), ('C', 'D')])
        
        stats = network_statistics(network)
        
        assert stats['num_nodes'] == 4
        assert stats['num_edges'] == 3
        assert stats['num_components'] == 1
        assert 'avg_degree' in stats
        assert 'density' in stats
    
    def test_statistics_empty(self):
        """Test statistics for empty network."""
        network = nx.Graph()
        
        stats = network_statistics(network)
        
        assert stats['num_nodes'] == 0
        assert stats['num_edges'] == 0


class TestValidateNetwork:
    """Test network validation."""
    
    def test_validate_good_network(self):
        """Test validation of good network."""
        network = nx.Graph()
        network.add_edges_from([('A', 'B'), ('B', 'C')])
        
        assert validate_network(network) is True
    
    def test_validate_empty_network(self):
        """Test validation of empty network."""
        network = nx.Graph()
        
        with pytest.raises(ValueError, match="no nodes"):
            validate_network(network)
    
    def test_validate_no_edges(self):
        """Test validation of network with no edges."""
        network = nx.Graph()
        network.add_node('A')
        
        with pytest.raises(ValueError, match="no edges"):
            validate_network(network)
    
    def test_validate_disconnected(self):
        """Test validation warns about disconnected network."""
        network = nx.Graph()
        network.add_edge('A', 'B')
        network.add_edge('C', 'D')
        
        with pytest.warns(UserWarning, match="connected components"):
            validate_network(network)
    
    def test_validate_self_loops(self):
        """Test validation warns about self-loops."""
        network = nx.Graph()
        network.add_edge('A', 'A')
        network.add_edge('A', 'B')
        
        with pytest.warns(UserWarning, match="self-loops"):
            validate_network(network)


if __name__ == '__main__':
    pytest.main([__file__, '-v'])

