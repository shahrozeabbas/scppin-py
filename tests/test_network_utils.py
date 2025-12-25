"""Tests for network utilities."""

import pytest
import igraph as ig
from scppin.core.network_utils import (
    simplify_network,
    filter_network,
    get_largest_connected_component,
    network_statistics,
    validate_network
)


class TestSimplifyNetwork:
    """Test network simplification."""
    
    def test_remove_self_loops(self):
        """Test removal of self-loops."""
        # Create graph with self-loop
        network = ig.Graph.TupleList([('A', 'A'), ('A', 'B')], directed=False, vertex_name_attr='name')
        
        simplified = simplify_network(network)
        
        # Check no self-loops
        assert simplified.get_eid('A', 'A', error=False) == -1
        assert simplified.get_eid('A', 'B', error=False) != -1
    
    def test_no_changes_needed(self):
        """Test when network is already simple."""
        network = ig.Graph.TupleList([('A', 'B'), ('B', 'C')], directed=False, vertex_name_attr='name')
        
        simplified = simplify_network(network)
        
        assert simplified.vcount() == 3
        assert simplified.ecount() == 2


class TestFilterNetwork:
    """Test network filtering."""
    
    def test_filter_keep_genes(self):
        """Test filtering to keep specific genes."""
        network = ig.Graph.TupleList([('A', 'B'), ('B', 'C'), ('C', 'D')], directed=False, vertex_name_attr='name')
        
        genes_to_keep = {'A', 'B', 'C'}
        
        filtered = filter_network(network, genes_to_keep, keep_only_genes=True)
        
        assert filtered.vcount() == 3
        node_names = filtered.vs['name']
        assert 'D' not in node_names
    
    def test_filter_remove_genes(self):
        """Test filtering to remove specific genes."""
        network = ig.Graph.TupleList([('A', 'B'), ('B', 'C'), ('C', 'D')], directed=False, vertex_name_attr='name')
        
        genes_to_remove = {'D'}
        
        filtered = filter_network(network, genes_to_remove, keep_only_genes=False)
        
        assert filtered.vcount() == 3
        node_names = filtered.vs['name']
        assert 'D' not in node_names
    
    def test_filter_empty_result(self):
        """Test filtering that results in empty network."""
        network = ig.Graph.TupleList([('A', 'B')], directed=False, vertex_name_attr='name')
        
        with pytest.warns(UserWarning, match="No nodes remain"):
            filtered = filter_network(network, {'C', 'D'}, keep_only_genes=True)
        
        assert filtered.vcount() == 0


class TestFilterNetworkByPvalues:
    """Test filtering by p-values."""
    
    def test_filter_with_pvalues(self):
        """Test filtering to genes with p-values."""
        network = ig.Graph.TupleList([('A', 'B'), ('B', 'C'), ('C', 'D')], directed=False, vertex_name_attr='name')
        
        pvalues = {'A': 0.01, 'B': 0.05}
        
        filtered = filter_network(network, set(pvalues.keys()))
        
        assert filtered.vcount() == 2
        node_names = filtered.vs['name']
        assert 'A' in node_names
        assert 'B' in node_names


class TestLargestConnectedComponent:
    """Test largest connected component extraction."""
    
    def test_single_component(self):
        """Test when network is already connected."""
        network = ig.Graph.TupleList([('A', 'B'), ('B', 'C')], directed=False, vertex_name_attr='name')
        
        lcc = get_largest_connected_component(network)
        
        assert lcc.vcount() == 3
    
    def test_multiple_components(self):
        """Test extracting largest from multiple components."""
        # Create disconnected graph: A-B-C and D-E
        network = ig.Graph.TupleList([('A', 'B'), ('B', 'C'), ('D', 'E')], directed=False, vertex_name_attr='name')
        
        with pytest.warns(UserWarning, match="largest component"):
            lcc = get_largest_connected_component(network)
        
        assert lcc.vcount() == 3
        node_names = lcc.vs['name']
        assert 'A' in node_names
        assert 'D' not in node_names


class TestNetworkStatistics:
    """Test network statistics computation."""
    
    def test_statistics_basic(self):
        """Test basic statistics computation."""
        network = ig.Graph.TupleList([('A', 'B'), ('B', 'C'), ('C', 'D')], directed=False, vertex_name_attr='name')
        
        stats = network_statistics(network)
        
        assert stats['num_nodes'] == 4
        assert stats['num_edges'] == 3
        assert stats['num_components'] == 1
        assert 'avg_degree' in stats
        assert 'density' in stats
    
    def test_statistics_empty(self):
        """Test statistics for empty network."""
        network = ig.Graph()
        
        stats = network_statistics(network)
        
        assert stats['num_nodes'] == 0
        assert stats['num_edges'] == 0


class TestValidateNetwork:
    """Test network validation."""
    
    def test_validate_good_network(self):
        """Test validation of good network."""
        network = ig.Graph.TupleList([('A', 'B'), ('B', 'C')], directed=False, vertex_name_attr='name')
        
        validated = validate_network(network)
        assert validated.vcount() == 3
        assert validated.ecount() == 2
    
    def test_validate_empty_network(self):
        """Test validation of empty network."""
        network = ig.Graph()
        
        with pytest.raises(ValueError, match="no nodes"):
            validate_network(network)
    
    def test_validate_no_edges(self):
        """Test validation of network with no edges."""
        network = ig.Graph()
        network.add_vertex('A')
        network.vs[0]['name'] = 'A'
        
        with pytest.raises(ValueError, match="no edges"):
            validate_network(network)
    
    def test_validate_disconnected(self):
        """Test validation warns about disconnected network."""
        network = ig.Graph.TupleList([('A', 'B'), ('C', 'D')], directed=False, vertex_name_attr='name')
        
        with pytest.warns(UserWarning, match="largest component"):
            validate_network(network)
    
    def test_validate_self_loops(self):
        """Test validation warns about self-loops."""
        network = ig.Graph.TupleList([('A', 'A'), ('A', 'B')], directed=False, vertex_name_attr='name')
        
        with pytest.warns(UserWarning, match="self-loops"):
            validated = validate_network(network)
        
        # Should still return the network
        assert validated.vcount() == 2


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
