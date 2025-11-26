"""
Example demonstrating the new class-based scPPIN API

This script shows the basic workflow with the new scPPIN class.
"""

import sys
sys.path.insert(0, '..')

import numpy as np
from scppin import scPPIN


def main():
    """Run example with new class-based API."""
    print("=" * 70)
    print("scPPIN Class-Based API Example")
    print("=" * 70)
    
    # Create analyzer instance
    analyzer = scPPIN()
    
    # Step 1: Load network from edge list
    print("\n1. Loading network from edge list...")
    edges = [
        ('TP53', 'MDM2'),
        ('TP53', 'CDKN1A'),
        ('TP53', 'PSEN1'),
        ('MDM2', 'CDKN1A'),
        ('APOE', 'CLU'),
        ('APOE', 'ABCA1'),
    ]
    
    analyzer.load_network(edges)
    print(f"   Network: {analyzer.network.number_of_nodes()} nodes, "
          f"{analyzer.network.number_of_edges()} edges")
    
    # Step 2: Set node weights (p-values)
    print("\n2. Setting node weights (p-values)...")
    pvalues = {
        'TP53': 0.0001,
        'MDM2': 0.0005,
        'CDKN1A': 0.001,
        'PSEN1': 0.005,
        'APOE': 0.01,
        'CLU': 0.02,
        'ABCA1': 0.05,
    }
    analyzer.set_node_weights(pvalues)
    print(f"   Node weights for {len(analyzer.node_weights)} genes")
    
    # Step 3: Detect module
    print("\n3. Detecting functional module...")
    fdr = 0.01
    
    try:
        module = analyzer.detect_module(fdr=fdr)
        
        print(f"   Module detected!")
        print(f"   Nodes: {module.number_of_nodes()}")
        print(f"   Edges: {module.number_of_edges()}")
        print(f"   Genes: {sorted(list(module.nodes()))}")
        
        # Print node scores
        print("\n   Node scores:")
        for node in sorted(module.nodes()):
            score = module.nodes[node].get('score', 'N/A')
            if isinstance(score, float):
                print(f"     {node}: {score:.4f}")
            else:
                print(f"     {node}: {score}")
        
    except ImportError as e:
        print(f"   Error: {e}")
        print("   Note: This example requires pcst_fast to be installed.")
        print("   Install with: pip install pcst-fast")
        return
    
    # Step 4: Example with edge weights
    print("\n" + "=" * 70)
    print("Example with Edge Weights")
    print("=" * 70)
    
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
    
    print(f"\nNetwork with edge weights:")
    for (u, v), w in weights.items():
        print(f"  {u} - {v}: {w:.3f}")
    
    try:
        analyzer2.detect_module(
            fdr=fdr,
            edge_weight_attr='weight',
            c0=0.1
        )
        
        print(f"\nModule with edge weights detected!")
        print(f"  Nodes: {analyzer2.module.number_of_nodes()}")
        print(f"  Edges: {analyzer2.module.number_of_edges()}")
        print(f"  Genes: {sorted(list(analyzer2.module.nodes()))}")
        
    except Exception as e:
        print(f"Error with edge weights: {e}")
    
    # Step 5: Method chaining example
    print("\n" + "=" * 70)
    print("Method Chaining Example")
    print("=" * 70)
    
    try:
        analyzer3 = (scPPIN()
                    .load_network(edges)
                    .set_node_weights(pvalues))
        module3 = analyzer3.detect_module(fdr=0.01)
        
        print(f"\nMethod chaining worked!")
        print(f"  Module nodes: {module3.number_of_nodes()}")
    except ImportError:
        pass
    
    print("\n" + "=" * 70)
    print("Example complete!")
    print("=" * 70)
    print("\nKey features of class-based API:")
    print("  ✓ Object-oriented design with state management")
    print("  ✓ Method chaining for concise workflows")
    print("  ✓ Automatic network filtering to genes with p-values")
    print("  ✓ Edge weight support with normalization")
    print("  ✓ Easy to understand and use")


if __name__ == '__main__':
    main()
