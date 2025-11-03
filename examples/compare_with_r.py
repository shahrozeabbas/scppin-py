"""
Compare Python implementation with R implementation

This script demonstrates that the Python implementation produces
similar results to the original R implementation using the class-based API.
"""

import sys
sys.path.insert(0, '..')

import pandas as pd
from scppin import scPPIN
from scppin.core.network_utils import load_ppin
import matplotlib.pyplot as plt


def load_example_data():
    """Load example data from R package."""
    print("Loading example data...")
    
    # Load network
    network_path = '../scppin/data/networks/biogridHomoSapiens3.5.166.graphml'
    try:
        network = load_ppin(network_path)
        print(f"Network loaded: {network.number_of_nodes()} nodes, "
              f"{network.number_of_edges()} edges")
    except FileNotFoundError:
        print("Network file not found. Please copy from R package:")
        print("  cp ../R/inst/extdata/biogridHomoSapiens3.5.166.graphml "
              "../scppin/data/networks/")
        return None, None
    
    # Load p-values
    try:
        pvalues_df = pd.read_csv('examplePvalues.csv')
        pvalues = dict(zip(pvalues_df['gene'], pvalues_df['pVal']))
        print(f"P-values loaded: {len(pvalues)} genes")
    except FileNotFoundError:
        print("P-values file not found. Please copy from R package:")
        print("  cp ../R/inst/extdata/examplePvalues.csv .")
        return network, None
    
    return network, pvalues


def main():
    """Run comparison with R implementation."""
    print("=" * 70)
    print("scPPIN Python vs R Comparison")
    print("=" * 70)
    
    # Load data
    network, pvalues = load_example_data()
    
    if network is None or pvalues is None:
        print("\nCannot proceed without data files.")
        return
    
    # Create analyzer and detect module (same parameters as R tutorial)
    print("\n" + "=" * 70)
    print("Detecting functional module (FDR = 0.01)...")
    print("=" * 70)
    
    fdr = 0.01
    
    try:
        analyzer = scPPIN()
        analyzer.load_network(network)
        analyzer.set_node_weights(pvalues)
        analyzer.detect_module(fdr=fdr)
        
        module = analyzer.module
        
        print(f"\nPython implementation results:")
        print(f"  Module size: {module.number_of_nodes()} nodes")
        print(f"  Module edges: {module.number_of_edges()} edges")
        print(f"  Genes in module: {sorted(list(module.nodes()))}")
        
        print("\nExpected R implementation results (from tutorial):")
        print("  Module size: 3 nodes")
        print("  Genes: APP, ALDOB, SCD")
        
        # Check if we got the expected genes
        expected_genes = {'APP', 'ALDOB', 'SCD'}
        detected_genes = set(module.nodes())
        
        if expected_genes.issubset(detected_genes):
            print("\n✓ Python implementation includes all expected genes!")
        else:
            missing = expected_genes - detected_genes
            extra = detected_genes - expected_genes
            print(f"\n⚠ Differences detected:")
            if missing:
                print(f"  Missing genes: {missing}")
            if extra:
                print(f"  Extra genes: {extra}")
            print("\nNote: Small differences are expected due to:")
            print("  - Different PCST solver implementations")
            print("  - Numerical precision differences")
            print("  - Random initialization in optimization")
        
        # Visualize
        print("\n" + "=" * 70)
        print("Creating visualization...")
        print("=" * 70)
        
        fig = plt.figure(figsize=(10, 8))
        analyzer.plot_module(fdr=fdr, title='Functional Module (Python scPPIN)')
        plt.savefig('python_module_result.png', dpi=150, bbox_inches='tight')
        print("Saved visualization to 'python_module_result.png'")
        
    except ImportError as e:
        print(f"\nError: {e}")
        print("\nPlease install required packages:")
        print("  pip install pcst-fast")
    except Exception as e:
        print(f"\nError: {e}")
        import traceback
        traceback.print_exc()
    
    print("\n" + "=" * 70)
    print("Comparison complete!")
    print("=" * 70)


if __name__ == '__main__':
    main()
