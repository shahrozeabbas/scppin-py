#!/usr/bin/env python
"""
Verification script for scPPIN installation.

Run this script to verify that scPPIN is properly installed and working.
"""

import sys


def check_imports():
    """Check that all required modules can be imported."""
    print("=" * 70)
    print("Checking imports...")
    print("=" * 70)
    
    required_modules = [
        ('numpy', 'NumPy'),
        ('scipy', 'SciPy'),
        ('networkx', 'NetworkX'),
        ('matplotlib', 'Matplotlib'),
    ]
    
    optional_modules = [
        ('pcst_fast', 'pcst-fast'),
        ('scanpy', 'Scanpy'),
        ('anndata', 'AnnData'),
    ]
    
    all_good = True
    
    # Check required
    print("\nRequired modules:")
    for module, name in required_modules:
        try:
            __import__(module)
            print(f"  ✓ {name}")
        except ImportError:
            print(f"  ✗ {name} - MISSING (required)")
            all_good = False
    
    # Check optional
    print("\nOptional modules:")
    for module, name in optional_modules:
        try:
            __import__(module)
            print(f"  ✓ {name}")
        except ImportError:
            print(f"  - {name} - not installed (optional)")
    
    return all_good


def check_scppin():
    """Check that scPPIN can be imported."""
    print("\n" + "=" * 70)
    print("Checking scPPIN package...")
    print("=" * 70)
    
    try:
        import scppin
        print(f"\n✓ scPPIN imported successfully")
        print(f"  Version: {scppin.__version__}")
        
        # Check main API
        if hasattr(scppin, 'scPPIN'):
            print("\nMain API:")
            print("  ✓ scPPIN class available")
        else:
            print("  ✗ scPPIN class missing")
            return False
        
        return True
        
    except ImportError as e:
        print(f"\n✗ Failed to import scPPIN: {e}")
        print("\nTry installing with: pip install -e .")
        return False


def run_basic_test():
    """Run a basic functionality test."""
    print("\n" + "=" * 70)
    print("Running basic functionality test...")
    print("=" * 70)
    
    try:
        import scppin
        from scppin import scPPIN
        from scppin.core import fit_bum
        import networkx as nx
        import numpy as np
        
        # Create simple test network
        network = nx.Graph()
        network.add_edges_from([
            ('GENE1', 'GENE2'),
            ('GENE2', 'GENE3'),
            ('GENE3', 'GENE4'),
        ])
        
        # Create test p-values
        pvalues = {
            'GENE1': 0.001,
            'GENE2': 0.005,
            'GENE3': 0.01,
            'GENE4': 0.05,
        }
        
        print("\nTest network:")
        print(f"  Nodes: {network.number_of_nodes()}")
        print(f"  Edges: {network.number_of_edges()}")
        print(f"  P-values: {len(pvalues)} genes")
        
        # Test BUM fitting
        print("\nTesting BUM model fitting...")
        pvals_array = np.array(list(pvalues.values()))
        lambda_param, alpha, success = fit_bum(pvals_array)
        
        if success:
            print(f"  ✓ BUM fitting successful")
            print(f"    λ = {lambda_param:.4f}")
            print(f"    α = {alpha:.4f}")
        else:
            print(f"  ⚠ BUM fitting did not converge (may be OK for small sample)")
        
        # Test module detection (requires pcst_fast)
        print("\nTesting module detection...")
        try:
            analyzer = scPPIN()
            analyzer.load_network(network)
            analyzer.set_node_weights(pvalues)
            module = analyzer.detect_module(fdr=0.01)
            print(f"  ✓ Module detection successful")
            print(f"    Module nodes: {module.number_of_nodes()}")
            print(f"    Module edges: {module.number_of_edges()}")
            print(f"    Genes: {sorted(list(module.nodes()))}")
            
        except ImportError as e:
            print(f"  - Module detection skipped: {e}")
            print(f"    Install pcst-fast with: pip install pcst-fast")
            return True  # Not a failure, just missing optional dep
        
        return True
        
    except Exception as e:
        print(f"\n✗ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


def main():
    """Run all verification checks."""
    print("\n" + "=" * 70)
    print("scPPIN Installation Verification")
    print("=" * 70)
    
    # Check imports
    imports_ok = check_imports()
    
    if not imports_ok:
        print("\n" + "=" * 70)
        print("RESULT: FAILED - Missing required dependencies")
        print("=" * 70)
        print("\nInstall missing dependencies with:")
        print("  pip install numpy scipy networkx matplotlib")
        return 1
    
    # Check scPPIN
    scppin_ok = check_scppin()
    
    if not scppin_ok:
        print("\n" + "=" * 70)
        print("RESULT: FAILED - scPPIN not properly installed")
        print("=" * 70)
        print("\nInstall scPPIN with:")
        print("  cd scppin-py")
        print("  pip install -e .")
        return 1
    
    # Run basic test
    test_ok = run_basic_test()
    
    # Final result
    print("\n" + "=" * 70)
    if test_ok:
        print("RESULT: ✓ ALL CHECKS PASSED")
        print("=" * 70)
        print("\nscPPIN is properly installed and working!")
        print("\nNext steps:")
        print("  - Run examples: python examples/basic_tutorial.py")
        print("  - Run tests: pytest tests/ -v")
        print("  - Read documentation: README.md")
        return 0
    else:
        print("RESULT: ✗ SOME TESTS FAILED")
        print("=" * 70)
        print("\nSee error messages above for details.")
        return 1


if __name__ == '__main__':
    sys.exit(main())
