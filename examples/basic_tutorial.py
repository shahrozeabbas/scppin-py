"""
Basic Tutorial for scPPIN

This tutorial demonstrates the basic usage of scPPIN for detecting
functional modules in protein-protein interaction networks using
the class-based API.
"""

import numpy as np
import matplotlib.pyplot as plt

# Import scppin
import sys
sys.path.insert(0, '..')
from scppin import scPPIN


def main():
    """Run the basic tutorial."""
    print("=" * 60)
    print("scPPIN Basic Tutorial")
    print("=" * 60)
    
    # Create analyzer instance
    analyzer = scPPIN()
    
    # Step 1: Load network
    print("\n1. Loading network...")
    edges = [
        ('APP', 'APOE'),
        ('APP', 'PSEN1'),
        ('APP', 'PSEN2'),
        ('APOE', 'CLU'),
        ('APOE', 'ABCA1'),
        ('PSEN1', 'PSEN2'),
        ('MAPT', 'GSK3B'),
        ('MAPT', 'CDK5'),
        ('GSK3B', 'CDK5'),
        ('TNF', 'IL6'),
        ('TNF', 'IL1B'),
        ('IL6', 'IL1B'),
    ]
    
    analyzer.load_network(edges)
    print(f"   Network: {analyzer.network.number_of_nodes()} nodes, "
          f"{analyzer.network.number_of_edges()} edges")
    
    # Step 2: Set node weights (p-values)
    print("\n2. Setting node weights (p-values)...")
    pvalues = {
        # Significant genes (cluster 1)
        'APP': 0.0001,
        'APOE': 0.0005,
        'PSEN1': 0.001,
        'PSEN2': 0.005,
        'CLU': 0.01,
        'ABCA1': 0.02,
        
        # Moderately significant
        'MAPT': 0.05,
        'GSK3B': 0.08,
        'CDK5': 0.1,
        
        # Not significant
        'TNF': 0.5,
        'IL6': 0.7,
        'IL1B': 0.9,
    }
    
    analyzer.set_node_weights(pvalues)
    print(f"   Node weights for {len(analyzer.node_weights)} genes")
    
    # Step 3: Detect functional module
    print("\n" + "=" * 60)
    print("3. Detecting functional module...")
    print("=" * 60)
    
    fdr = 0.01
    
    try:
        analyzer.detect_module(fdr=fdr)
        
        print(f"\nModule detected!")
        print(f"  Nodes: {analyzer.module.number_of_nodes()}")
        print(f"  Edges: {analyzer.module.number_of_edges()}")
        print(f"  Genes in module: {list(analyzer.module.nodes())}")
        
        # Print node scores
        print("\nNode scores:")
        for node in analyzer.module.nodes():
            score = analyzer.module.nodes[node].get('score', 'N/A')
            print(f"  {node}: {score:.4f}" if isinstance(score, float) else f"  {node}: {score}")
        
    except Exception as e:
        print(f"\nError detecting module: {e}")
        print("Note: This example requires pcst_fast to be installed.")
        print("Install with: pip install pcst-fast")
        return
    
    # Step 4: Visualize (optional)
    print("\n" + "=" * 60)
    print("4. Visualizing module...")
    print("=" * 60)
    
    try:
        fig, axes = plt.subplots(1, 2, figsize=(14, 6))
        
        # Plot original network
        ax1 = axes[0]
        import networkx as nx
        pos = nx.spring_layout(analyzer.network, seed=42)
        nx.draw(analyzer.network, pos, with_labels=True, node_color='lightblue',
                node_size=500, font_size=8, ax=ax1)
        ax1.set_title('Original Network')
        ax1.axis('off')
        
        # Plot detected module
        ax2 = axes[1]
        analyzer.plot_module(fdr=fdr, ax=ax2)
        
        plt.tight_layout()
        plt.savefig('basic_tutorial_result.png', dpi=150, bbox_inches='tight')
        print("Visualization saved to 'basic_tutorial_result.png'")
        
    except Exception as e:
        print(f"Visualization skipped: {e}")
    
    # Step 5: Example with edge weights
    print("\n" + "=" * 60)
    print("5. Example with edge weights...")
    print("=" * 60)
    
    # Create new analyzer for edge weights example
    analyzer2 = scPPIN()
    analyzer2.load_network(edges)
    analyzer2.set_node_weights(pvalues)
    
    # Add confidence scores to edges
    np.random.seed(42)
    weights = {}
    for u, v in edges:
        # Higher confidence for edges between significant genes
        if pvalues.get(u, 1.0) < 0.01 and pvalues.get(v, 1.0) < 0.01:
            weights[(u, v)] = np.random.uniform(0.8, 1.0)
        else:
            weights[(u, v)] = np.random.uniform(0.3, 0.7)
    
    analyzer2.set_edge_weights(weights=weights)
    
    try:
        analyzer2.detect_module(fdr=fdr, edge_weight_scale=0.5)
        
        print(f"\nModule with edge weights detected!")
        print(f"  Nodes: {analyzer2.module.number_of_nodes()}")
        print(f"  Edges: {analyzer2.module.number_of_edges()}")
        print(f"  Genes in module: {list(analyzer2.module.nodes())}")
        
    except Exception as e:
        print(f"Error with edge weights: {e}")
    
    print("\n" + "=" * 60)
    print("Tutorial complete!")
    print("=" * 60)
    print("\nKey features:")
    print("  ✓ Class-based API with method chaining")
    print("  ✓ Automatic network filtering")
    print("  ✓ Edge weight support")
    print("  ✓ Easy visualization")


if __name__ == '__main__':
    main()
