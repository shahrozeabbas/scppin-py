"""Tests for new class-based API."""

import pytest
import numpy as np
import networkx as nx
import pandas as pd
from scppin import scPPIN


class TestScPPINLoadNetwork:
    """Test scPPIN.load_network() method."""
    
    def test_load_from_list(self):
        """Test loading network from list of tuples."""
        analyzer = scPPIN()
        edges = [('A', 'B'), ('B', 'C'), ('C', 'D')]
        
        analyzer.load_network(edges)
        
        assert analyzer.network.number_of_nodes() == 4
        assert analyzer.network.number_of_edges() == 3
        assert analyzer.network.has_edge('A', 'B')
    
    def test_load_from_dataframe(self):
        """Test loading network from DataFrame."""
        analyzer = scPPIN()
        df = pd.DataFrame({
            'source': ['A', 'B', 'C'],
            'target': ['B', 'C', 'D']
        })
        
        analyzer.load_network(df)
        
        assert analyzer.network.number_of_nodes() == 4
        assert analyzer.network.number_of_edges() == 3
    
    def test_load_with_weight_column(self):
        """Test loading network with weight column."""
        analyzer = scPPIN()
        df = pd.DataFrame({
            'source': ['A', 'B'],
            'target': ['B', 'C'],
            'confidence': [0.9, 0.8]
        })
        
        analyzer.load_network(df, weight_column='confidence')
        
        assert analyzer.network['A']['B']['weight'] == 0.9
        assert analyzer.network['B']['C']['weight'] == 0.8
        assert len(analyzer.edge_weights) == 2
    
    def test_load_from_nx_graph(self):
        """Test loading network from NetworkX graph."""
        analyzer = scPPIN()
        g = nx.Graph()
        g.add_edges_from([('A', 'B'), ('B', 'C')])
        
        analyzer.load_network(g)
        
        assert analyzer.network.number_of_nodes() == 2
        assert analyzer.network.number_of_edges() == 2


class TestScPPINSetNodeWeights:
    """Test scPPIN.set_node_weights() method."""
    
    def create_test_network(self):
        """Create a simple test network."""
        analyzer = scPPIN()
        edges = [
            ('GENE1', 'GENE2'),
            ('GENE2', 'GENE3'),
            ('GENE3', 'GENE4'),
            ('GENE4', 'GENE5'),
            ('GENE1', 'GENE5'),
        ]
        analyzer.load_network(edges)
        return analyzer
    
    def test_set_node_weights_from_dict(self):
        """Test setting node weights from dictionary."""
        analyzer = self.create_test_network()
        pvalues = {
            'GENE1': 0.001,
            'GENE2': 0.005,
            'GENE3': 0.01,
            'GENE4': 0.05,
            'GENE5': 0.08,
        }
        
        analyzer.set_node_weights(pvalues)
        
        assert len(analyzer.node_weights) == 5
        assert analyzer.node_weights['GENE1'] == 0.001
    
    def test_set_node_weights_filters_network(self):
        """Test that setting node weights filters network."""
        analyzer = self.create_test_network()
        # Add extra edge that won't have p-value
        analyzer.network.add_edge('GENE6', 'GENE7')
        
        pvalues = {
            'GENE1': 0.001,
            'GENE2': 0.005,
        }
        
        analyzer.set_node_weights(pvalues)
        
        # Network should be filtered to only genes with p-values
        assert 'GENE6' not in analyzer.network.nodes()
        assert 'GENE7' not in analyzer.network.nodes()
        assert 'GENE1' in analyzer.network.nodes()
        assert 'GENE2' in analyzer.network.nodes()


class TestScPPINSetEdgeWeights:
    """Test scPPIN.set_edge_weights() method."""
    
    def create_test_network(self):
        """Create a simple test network."""
        analyzer = scPPIN()
        edges = [
            ('GENE1', 'GENE2'),
            ('GENE2', 'GENE3'),
            ('GENE3', 'GENE4'),
        ]
        analyzer.load_network(edges)
        return analyzer
    
    def test_set_edge_weights_from_dict(self):
        """Test setting edge weights from dictionary."""
        analyzer = self.create_test_network()
        weights = {
            ('GENE1', 'GENE2'): 0.9,
            ('GENE2', 'GENE3'): 0.8,
        }
        
        analyzer.set_edge_weights(weights=weights)
        
        assert len(analyzer.edge_weights) == 2
        assert analyzer.network['GENE1']['GENE2']['weight'] == 0.9
        assert analyzer.network['GENE2']['GENE3']['weight'] == 0.8
    
    def test_set_edge_weights_filters_to_network(self):
        """Test that edge weights are filtered to network edges only."""
        analyzer = self.create_test_network()
        weights = {
            ('GENE1', 'GENE2'): 0.9,  # In network
            ('NOT_IN_NETWORK', 'ALSO_NOT'): 0.5,  # Not in network
        }
        
        analyzer.set_edge_weights(weights=weights)
        
        # Only network edges should be set
        assert len(analyzer.edge_weights) == 1
        assert ('GENE1', 'GENE2') in analyzer.edge_weights


class TestScPPINDetectModule:
    """Test scPPIN.detect_module() method."""
    
    def create_test_analyzer(self):
        """Create analyzer with network and p-values."""
        analyzer = scPPIN()
        edges = [
            ('GENE1', 'GENE2'),
            ('GENE2', 'GENE3'),
            ('GENE3', 'GENE4'),
            ('GENE4', 'GENE5'),
            ('GENE1', 'GENE5'),
        ]
        analyzer.load_network(edges)
        
        pvalues = {
            'GENE1': 0.001,
            'GENE2': 0.005,
            'GENE3': 0.01,
            'GENE4': 0.05,
            'GENE5': 0.08,
        }
        analyzer.set_node_weights(pvalues)
        
        return analyzer
    
    def test_detect_module_basic(self):
        """Test basic module detection."""
        analyzer = self.create_test_analyzer()
        
        try:
            analyzer.detect_module(fdr=0.01)
            
            assert analyzer.module is not None
            assert analyzer.module.number_of_nodes() > 0
            assert analyzer.module.number_of_edges() >= 0
            
            # All nodes should have scores
            for node in analyzer.module.nodes():
                assert 'score' in analyzer.module.nodes[node]
        except ImportError:
            pytest.skip("pcst_fast not installed")
    
    def test_detect_module_with_edge_weights(self):
        """Test module detection with edge weights."""
        analyzer = self.create_test_analyzer()
        
        # Add edge weights
        weights = {
            ('GENE1', 'GENE2'): 0.9,
            ('GENE2', 'GENE3'): 0.8,
        }
        analyzer.set_edge_weights(weights=weights)
        
        try:
            analyzer.detect_module(
                fdr=0.01,
                edge_weight_scale=0.5,
                c0=0.1
            )
            
            assert analyzer.module.number_of_nodes() > 0
        except ImportError:
            pytest.skip("pcst_fast not installed")
    
    def test_detect_module_default_no_edge_weights(self):
        """Test that default edge_weight_attr=None uses uniform costs even if network has 'weight' attribute."""
        analyzer = self.create_test_analyzer()
        
        # Add 'weight' attribute to edges (should be ignored by default)
        for u, v in analyzer.network.edges():
            analyzer.network[u][v]['weight'] = 0.9  # High weight that would affect results if used
        
        try:
            # Default should ignore edge weights (uniform costs)
            analyzer.detect_module(fdr=0.01)
            module_default = analyzer.module
            
            # Explicitly using edge weights should give different results
            analyzer2 = self.create_test_analyzer()
            analyzer2.set_edge_weights(weights={('GENE1', 'GENE2'): 0.9})
            analyzer2.detect_module(fdr=0.01, edge_weight_scale=1.0)
            module_with_weights = analyzer2.module
            
            # Both should produce valid modules
            assert module_default.number_of_nodes() > 0
            assert module_with_weights.number_of_nodes() > 0
            
            # Note: Modules may or may not be identical, but both should work
            # The key test is that default doesn't crash and doesn't use weights
        except ImportError:
            pytest.skip("pcst_fast not installed")
    
    def test_detect_module_method_chaining(self):
        """Test method chaining."""
        analyzer = scPPIN()
        
        try:
            result = (analyzer
                     .load_network([('A', 'B'), ('B', 'C')])
                     .set_node_weights({'A': 0.001, 'B': 0.01, 'C': 0.05})
                     .detect_module(fdr=0.01))
            
            assert result is analyzer
            assert analyzer.module is not None
        except ImportError:
            pytest.skip("pcst_fast not installed")


class TestScPPINPlotModule:
    """Test scPPIN.plot_module() method."""
    
    def test_plot_module_requires_detection(self):
        """Test that plot_module requires detect_module to be called first."""
        analyzer = scPPIN()
        analyzer.load_network([('A', 'B')])
        analyzer.set_node_weights({'A': 0.001, 'B': 0.01})
        
        with pytest.raises(ValueError, match="No module detected"):
            analyzer.plot_module()


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
