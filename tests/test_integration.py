"""Integration tests for scppin."""

import pytest
import numpy as np
import networkx as nx
import warnings
from scppin import scPPIN
from scppin.core import fit_bum


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
            assert analyzer.module.number_of_nodes() > 0
            assert analyzer.module.number_of_edges() >= 0
            
            # All nodes should have scores
            for node in analyzer.module.nodes():
                assert 'score' in analyzer.module.nodes[node]
        except ImportError:
            pytest.skip("pcst_fast not installed")
    
    def test_module_with_edge_weights(self):
        """Test module detection with edge weights."""
        analyzer = self.create_test_analyzer()
        
        # Add edge weights
        weights = {}
        for u, v in analyzer.network.edges():
            weights[(u, v)] = np.random.uniform(0.5, 1.0)
        analyzer.set_edge_weights(weights=weights)
        
        try:
            analyzer.detect_module(
                fdr=0.01,
                edge_weight_scale=0.5
            )
            
            assert analyzer.module.number_of_nodes() > 0
        except ImportError:
            pytest.skip("pcst_fast not installed")
    
    def test_module_with_missing_data(self):
        """Test module detection with missing data."""
        analyzer = scPPIN()
        edges = [
            ('GENE1', 'GENE2'),
            ('GENE2', 'GENE3'),
            ('GENE3', 'GENE4'),
            ('GENE4', 'GENE5'),
        ]
        analyzer.load_network(edges)
        
        pvalues = {
            'GENE1': 0.001,
            'GENE2': 0.005,
            # GENE3, GENE4, GENE5 missing
        }
        analyzer.set_node_weights(pvalues)
        
        # Should work with missing_data_score=True
        try:
            analyzer.detect_module(
                fdr=0.01,
                missing_data_score=True
            )
            
            # Module might include genes without p-values
            assert analyzer.module.number_of_nodes() >= 2
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
            with pytest.raises(ValueError, match="No genes in network"):
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
            assert module_strict.number_of_nodes() > 0
            assert module_lenient.number_of_nodes() > 0
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
        
        # This should produce a warning and return empty graph with metadata
        try:
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                analyzer.detect_module(fdr=0.01)
                module = analyzer.module
                
                # Check if empty solution occurred (may or may not happen depending on network)
                if module.number_of_nodes() == 0:
                    # Verify warning was issued
                    assert len(w) > 0
                    assert any("empty solution" in str(warning.message).lower() 
                               for warning in w)
                    
                    # Verify metadata
                    assert module.graph.get('empty_solution', False) == True
                    assert module.graph.get('reason') == 'pcst_returned_empty'
                else:
                    # If solution is not empty, that's also valid
                    # Just verify module is well-formed
                    assert module.number_of_nodes() > 0
        except ImportError:
            pytest.skip("pcst_fast not installed")
    
    def test_method_chaining(self):
        """Test method chaining workflow."""
        try:
            analyzer = (scPPIN()
                       .load_network([('GENE1', 'GENE2'), ('GENE2', 'GENE3')])
                       .set_node_weights({'GENE1': 0.001, 'GENE2': 0.01, 'GENE3': 0.05})
                       .detect_module(fdr=0.01))
            
            assert analyzer.module is not None
            assert analyzer.module.number_of_nodes() > 0
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
        
        # For this mixture, lambda should be around 0.8
        assert 0.6 < lambda_param < 0.95


if __name__ == '__main__':
    pytest.main([__file__, '-v'])
