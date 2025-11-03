"""
Edge Weights Example

This example demonstrates edge weight functionality using the class-based API.
"""

import sys
sys.path.insert(0, '..')

import numpy as np
import matplotlib.pyplot as plt
from scppin import scPPIN


def create_example_with_edge_weights():
    """Create example network with edge confidence scores."""
    print("Creating example network with edge weights...")
    
    # High confidence edges (from databases like STRING)
    high_conf_edges = [
        ('TP53', 'MDM2', 0.95),
        ('TP53', 'CDKN1A', 0.90),
        ('MDM2', 'CDKN1A', 0.85),
        ('BRCA1', 'BRCA2', 0.92),
        ('BRCA1', 'RAD51', 0.88),
    ]
    
    # Medium confidence edges
    med_conf_edges = [
        ('TP53', 'BRCA1', 0.65),
        ('CDKN1A', 'BRCA1', 0.60),
        ('BRCA2', 'RAD51', 0.70),
    ]
    
    # Low confidence edges
    low_conf_edges = [
        ('MDM2', 'BRCA2', 0.35),
        ('RAD51', 'TP53', 0.40),
    ]
    
    all_edges = [(u, v) for u, v, _ in high_conf_edges + med_conf_edges + low_conf_edges]
    weights = {(u, v): conf for u, v, conf in high_conf_edges + med_conf_edges + low_conf_edges}
    
    print(f"Network: {len(set(sum(all_edges, ())))} nodes, {len(all_edges)} edges")
    
    return all_edges, weights


def create_pvalues():
    """Create example p-values."""
    # Simulate differential expression
    pvalues = {
        # TP53 pathway - highly significant
        'TP53': 0.0001,
        'MDM2': 0.0005,
        'CDKN1A': 0.001,
        
        # BRCA pathway - moderately significant
        'BRCA1': 0.01,
        'BRCA2': 0.02,
        'RAD51': 0.03,
    }
    
    return pvalues


def main():
    """Run edge weights example."""
    print("=" * 70)
    print("Edge Weights Example (Class-Based API)")
    print("=" * 70)
    
    edges, weights_dict = create_example_with_edge_weights()
    pvalues = create_pvalues()
    
    fdr = 0.01
    
    # Example 1: Without edge weights (uniform costs)
    print("\n" + "=" * 70)
    print("Example 1: Without edge weights (baseline)")
    print("=" * 70)
    
    try:
        analyzer_uniform = scPPIN()
        analyzer_uniform.load_network(edges)
        analyzer_uniform.set_node_weights(pvalues)
        analyzer_uniform.detect_module(fdr=fdr)
        
        print(f"Module without edge weights:")
        print(f"  Nodes: {analyzer_uniform.module.number_of_nodes()}")
        print(f"  Edges: {analyzer_uniform.module.number_of_edges()}")
        print(f"  Genes: {sorted(list(analyzer_uniform.module.nodes()))}")
        
    except Exception as e:
        print(f"Error: {e}")
        print("Note: Requires pcst-fast package")
        return
    
    # Example 2: With edge weights (confidence-based)
    print("\n" + "=" * 70)
    print("Example 2: With edge weights (confidence-based)")
    print("=" * 70)
    
    analyzer_weighted = scPPIN()
    analyzer_weighted.load_network(edges)
    analyzer_weighted.set_node_weights(pvalues)
    analyzer_weighted.set_edge_weights(weights=weights_dict)
    analyzer_weighted.detect_module(fdr=fdr, edge_weight_scale=0.5)
    
    print(f"Module with edge weights (scale=0.5):")
    print(f"  Nodes: {analyzer_weighted.module.number_of_nodes()}")
    print(f"  Edges: {analyzer_weighted.module.number_of_edges()}")
    print(f"  Genes: {sorted(list(analyzer_weighted.module.nodes()))}")
    
    # Example 3: Strong edge weight influence
    print("\n" + "=" * 70)
    print("Example 3: Strong edge weight influence")
    print("=" * 70)
    
    analyzer_strong = scPPIN()
    analyzer_strong.load_network(edges)
    analyzer_strong.set_node_weights(pvalues)
    analyzer_strong.set_edge_weights(weights=weights_dict)
    analyzer_strong.detect_module(fdr=fdr, edge_weight_scale=1.0)
    
    print(f"Module with strong edge weights (scale=1.0):")
    print(f"  Nodes: {analyzer_strong.module.number_of_nodes()}")
    print(f"  Edges: {analyzer_strong.module.number_of_edges()}")
    print(f"  Genes: {sorted(list(analyzer_strong.module.nodes()))}")
    
    # Visualize comparison
    print("\n" + "=" * 70)
    print("Creating comparison visualization...")
    print("=" * 70)
    
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    # Plot 1: Without edge weights
    analyzer_uniform.plot_module(
        fdr=fdr,
        title='Without Edge Weights',
        ax=axes[0],
        show_colorbar=False
    )
    
    # Plot 2: With moderate edge weights
    analyzer_weighted.plot_module(
        fdr=fdr,
        title='With Edge Weights (scale=0.5)',
        ax=axes[1],
        show_colorbar=False
    )
    
    # Plot 3: With strong edge weights
    analyzer_strong.plot_module(
        fdr=fdr,
        title='With Edge Weights (scale=1.0)',
        ax=axes[2],
        show_colorbar=False
    )
    
    plt.tight_layout()
    plt.savefig('edge_weights_comparison.png', dpi=150, bbox_inches='tight')
    print("Saved comparison to 'edge_weights_comparison.png'")
    
    # Summary
    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)
    print("\nEdge weights allow you to:")
    print("  1. Incorporate interaction confidence from databases (STRING, BioGRID)")
    print("  2. Use co-expression data to weight edges")
    print("  3. Prioritize high-confidence interactions in module detection")
    print("\nHigher edge weights → Lower costs → More likely to include in module")
    print("=" * 70)


if __name__ == '__main__':
    main()
