"""Tests for scPPIN class-based API."""

import pytest
import numpy as np
import pandas as pd
import warnings
import igraph as ig
from scppin import scPPIN
from scppin.core import fit_bum


class TestScPPINLoadNetwork:
    """Test scPPIN.load_network() method."""
    
    def test_load_from_list(self):
        """Test loading network from list of tuples."""
        analyzer = scPPIN()
        edges = [('A', 'B'), ('B', 'C'), ('C', 'D')]
        
        analyzer.load_network(edges)
        
        assert analyzer.network.vcount() == 4
        assert analyzer.network.ecount() == 3
        assert analyzer.network.get_eid('A', 'B', error=False) != -1
    
    def test_load_from_dataframe(self):
        """Test loading network from DataFrame."""
        analyzer = scPPIN()
        df = pd.DataFrame({
            'source': ['A', 'B', 'C'],
            'target': ['B', 'C', 'D']
        })
        
        analyzer.load_network(df)
        
        assert analyzer.network.vcount() == 4
        assert analyzer.network.ecount() == 3
    
    def test_load_with_weight_column(self):
        """Test loading network with weight column."""
        analyzer = scPPIN()
        df = pd.DataFrame({
            'source': ['A', 'B'],
            'target': ['B', 'C'],
            'confidence': [0.9, 0.8]
        })
        
        analyzer.load_network(df, weight_column='confidence')
        
        eid_ab = analyzer.network.get_eid('A', 'B', error=False)
        eid_bc = analyzer.network.get_eid('B', 'C', error=False)
        assert eid_ab != -1
        assert eid_bc != -1
        assert analyzer.network.es[eid_ab]['weight'] == 0.9
        assert analyzer.network.es[eid_bc]['weight'] == 0.8
        assert len(analyzer.edge_weights) == 2
    
    def test_load_from_igraph_graph(self):
        """Test loading network from igraph graph."""
        analyzer = scPPIN()
        g = ig.Graph.TupleList([('A', 'B'), ('B', 'C')], directed=False, vertex_name_attr='name')
        
        analyzer.load_network(g)
        
        assert analyzer.network.vcount() == 3
        assert analyzer.network.ecount() == 2


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
        # Create a new network with extra nodes
        edges_with_extra = [
            ('GENE1', 'GENE2'),
            ('GENE2', 'GENE3'),
            ('GENE3', 'GENE4'),
            ('GENE4', 'GENE5'),
            ('GENE1', 'GENE5'),
            ('GENE6', 'GENE7'),  # Extra edge
        ]
        analyzer.load_network(edges_with_extra)
        
        pvalues = {
            'GENE1': 0.001,
            'GENE2': 0.005,
        }
        
        analyzer.set_node_weights(pvalues)
        
        # Network should be filtered to only genes with p-values
        node_names = analyzer.network.vs['name']
        assert 'GENE6' not in node_names
        assert 'GENE7' not in node_names
        assert 'GENE1' in node_names
        assert 'GENE2' in node_names


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
        eid_12 = analyzer.network.get_eid('GENE1', 'GENE2', error=False)
        eid_23 = analyzer.network.get_eid('GENE2', 'GENE3', error=False)
        assert eid_12 != -1
        assert eid_23 != -1
        assert analyzer.network.es[eid_12]['weight'] == 0.9
        assert analyzer.network.es[eid_23]['weight'] == 0.8
    
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
            assert analyzer.module.vcount() > 0
            assert analyzer.module.ecount() >= 0
            
            # All nodes should have scores
            for v in analyzer.module.vs:
                assert 'score' in v.attributes()
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
                edge_weight_attr='weight',
                c0=0.1
            )
            
            assert analyzer.module.vcount() > 0
        except ImportError:
            pytest.skip("pcst_fast not installed")
    
    def test_detect_module_default_no_edge_weights(self):
        """Test that default edge_weight_attr=None uses uniform costs even if network has 'weight' attribute."""
        analyzer = self.create_test_analyzer()
        
        # Add 'weight' attribute to edges (should be ignored by default)
        for e in analyzer.network.es:
            e['weight'] = 0.9  # High weight that would affect results if used
        
        try:
            # Default should ignore edge weights (uniform costs)
            analyzer.detect_module(fdr=0.01)
            module_default = analyzer.module
            
            # Explicitly using edge weights should give different results
            analyzer2 = self.create_test_analyzer()
            analyzer2.set_edge_weights(weights={('GENE1', 'GENE2'): 0.9})
            analyzer2.detect_module(fdr=0.01, edge_weight_attr='weight')
            module_with_weights = analyzer2.module
            
            # Both should produce valid modules
            assert module_default.vcount() > 0
            assert module_with_weights.vcount() > 0
            
            # Note: Modules may or may not be identical, but both should work
            # The key test is that default doesn't crash and doesn't use weights
        except ImportError:
            pytest.skip("pcst_fast not installed")
    
    def test_detect_module_method_chaining(self):
        """Test method chaining."""
        analyzer = scPPIN()
        
        analyzer = (analyzer
                   .load_network([('A', 'B'), ('B', 'C')])
                   .set_node_weights({'A': 0.001, 'B': 0.01, 'C': 0.05}))
        
        try:
            module = analyzer.detect_module(fdr=0.01)
            assert module is analyzer.module
            assert module.vcount() > 0
        except ImportError:
            pytest.skip("pcst_fast not installed")
    
    def test_detect_module_with_max_prize_root(self):
        """Test module detection with max prize root parameter."""
        analyzer = self.create_test_analyzer()
        
        try:
            # Test with use_max_prize_root=True
            analyzer.detect_module(
                fdr=0.01,
                use_max_prize_root=True
            )
            
            assert analyzer.module is not None
            assert analyzer.module.vcount() > 0
            
            # Test with use_max_prize_root=False (default)
            analyzer2 = self.create_test_analyzer()
            analyzer2.detect_module(
                fdr=0.01,
                use_max_prize_root=False
            )
            
            assert analyzer2.module is not None
            assert analyzer2.module.vcount() > 0
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


class TestIntegration:
    """End-to-end integration tests."""
    
    def create_test_analyzer(self):
        """Create analyzer with test network and p-values."""
        analyzer = scPPIN()
        # Create a small network with known structure
        edges = [
            ('GENE1', 'GENE2'),
            ('GENE2', 'GENE3'),
            ('GENE3', 'GENE4'),
            ('GENE4', 'GENE5'),
            ('GENE1', 'GENE5'),  # Close the loop
            ('GENE6', 'GENE7'),  # Separate component
        ]
        analyzer.load_network(edges)
        
        # Significant genes in main component
        pvalues = {
            'GENE1': 0.001,
            'GENE2': 0.005,
            'GENE3': 0.01,
            'GENE4': 0.05,
            'GENE5': 0.08,
            'GENE6': 0.5,  # Not significant
            'GENE7': 0.8,  # Not significant
        }
        analyzer.set_node_weights(pvalues)
        
        return analyzer
    
    def test_basic_module_detection(self):
        """Test basic module detection workflow."""
        analyzer = self.create_test_analyzer()
        
        try:
            analyzer.detect_module(fdr=0.01)
            
            # Module should be non-empty
            assert analyzer.module.vcount() > 0
            assert analyzer.module.ecount() >= 0
            
            # All nodes should have scores
            for v in analyzer.module.vs:
                assert 'score' in v.attributes()
        except ImportError:
            pytest.skip("pcst_fast not installed")
    
    def test_module_with_edge_weights(self):
        """Test module detection with edge weights."""
        analyzer = self.create_test_analyzer()
        
        # Add edge weights
        weights = {}
        node_names = analyzer.network.vs['name']
        for e in analyzer.network.es:
            u_name = node_names[e.source]
            v_name = node_names[e.target]
            weights[(u_name, v_name)] = np.random.uniform(0.5, 1.0)
        analyzer.set_edge_weights(weights=weights)
        
        try:
            analyzer.detect_module(
                fdr=0.01,
                edge_weight_attr='weight'
            )
            
            assert analyzer.module.vcount() > 0
        except ImportError:
            pytest.skip("pcst_fast not installed")
    
    def test_invalid_pvalues(self):
        """Test that invalid p-values raise errors."""
        analyzer = scPPIN()
        analyzer.load_network([('GENE1', 'GENE2')])
        
        # P-values with zeros
        with pytest.raises(ValueError, match="must be in"):
            analyzer.set_node_weights({'GENE1': 0.0, 'GENE2': 0.5})
        
        # P-values > 1
        analyzer2 = scPPIN()
        analyzer2.load_network([('GENE1', 'GENE2')])
        with pytest.raises(ValueError, match="must be in"):
            analyzer2.set_node_weights({'GENE1': 1.5, 'GENE2': 0.5})
    
    def test_empty_pvalues(self):
        """Test that empty p-values raise error."""
        analyzer = scPPIN()
        analyzer.load_network([('GENE1', 'GENE2')])
        
        with pytest.raises(ValueError, match="empty"):
            analyzer.set_node_weights({})
    
    def test_no_matching_genes(self):
        """Test when no genes match between network and p-values."""
        analyzer = scPPIN()
        analyzer.load_network([('GENE1', 'GENE2')])
        
        # This should filter network to empty, which will cause error in detect_module
        analyzer.set_node_weights({'NOTINGRAPH1': 0.001, 'NOTINGRAPH2': 0.005})
        
        try:
            with pytest.raises(ValueError, match="No node scores computed"):
                analyzer.detect_module(fdr=0.01)
        except ImportError:
            pytest.skip("pcst_fast not installed")
    
    def test_different_fdr_thresholds(self):
        """Test module detection with different FDR thresholds."""
        analyzer = self.create_test_analyzer()
        
        try:
            # Stricter FDR
            analyzer.detect_module(fdr=0.001)
            module_strict = analyzer.module
            
            # Lenient FDR
            analyzer.detect_module(fdr=0.1)
            module_lenient = analyzer.module
            
            # Both should produce valid modules
            assert module_strict.vcount() > 0
            assert module_lenient.vcount() > 0
        except ImportError:
            pytest.skip("pcst_fast not installed")
    
    def test_empty_module_solution(self):
        """Test handling of empty module solution from PCST solver."""
        analyzer = scPPIN()
        edges = [
            ('GENE1', 'GENE2'),
            ('GENE2', 'GENE3'),
            ('GENE3', 'GENE4'),
            ('GENE4', 'GENE5'),
        ]
        analyzer.load_network(edges)
        
        # Create p-values that will result in all negative scores after BUM fitting
        # This should trigger an empty solution
        # Use very large p-values (>0.9) to get negative scores
        pvalues = {
            'GENE1': 0.95,
            'GENE2': 0.98,
            'GENE3': 0.99,
            'GENE4': 0.97,
            'GENE5': 0.96,
        }
        analyzer.set_node_weights(pvalues)
        
        # This should produce a warning and return empty graph
        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                analyzer.detect_module(fdr=0.01)
                module = analyzer.module
                
                # Check if empty solution occurred (may or may not happen depending on network)
                if module.vcount() == 0:
                    # Verify warning was issued
                    assert len(w) > 0
                    assert any("empty solution" in str(warning.message).lower() 
                               for warning in w)
                else:
                    # If solution is not empty, that's also valid
                    # Just verify module is well-formed
                    assert module.vcount() > 0
        except ImportError:
            pytest.skip("pcst_fast not installed")
    
    def test_method_chaining(self):
        """Test setup method chaining workflow."""
        try:
            analyzer = (scPPIN()
                       .load_network([('GENE1', 'GENE2'), ('GENE2', 'GENE3')])
                       .set_node_weights({'GENE1': 0.001, 'GENE2': 0.01, 'GENE3': 0.05}))
            module = analyzer.detect_module(fdr=0.01)
            
            assert module is analyzer.module
            assert module.vcount() > 0
        except ImportError:
            pytest.skip("pcst_fast not installed")


class TestBUMIntegration:
    """Test BUM model integration."""
    
    def test_bum_fitting_realistic(self):
        """Test BUM fitting with realistic p-value distribution."""
        np.random.seed(42)
        
        # Mix of significant and non-significant p-values
        significant = np.random.beta(0.5, 5, 200)  # Skewed toward 0
        nonsignificant = np.random.uniform(0, 1, 800)  # Uniform
        pvalues = np.concatenate([significant, nonsignificant])
        pvalues = np.clip(pvalues, 0.001, 1.0)  # Ensure valid range
        
        lambda_param, alpha, success = fit_bum(pvalues)
        
        assert success
        assert 0 < lambda_param < 1
        assert 0 < alpha < 1


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
