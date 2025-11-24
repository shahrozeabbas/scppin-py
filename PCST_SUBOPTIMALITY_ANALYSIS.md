# PCST Suboptimality Analysis

## Executive Summary

The PCST solver is returning suboptimal modules with a net benefit of **2.18** when greedy algorithms find solutions with net benefit of **108.37** (4873% better). The problem appears to be fundamental to how PCST is being parameterized or called, not a tuning issue with costs or FDR.

## Critical Findings

### 1. PCST is Stuck at Local Optimum

**Current PCST Solution (FDR=0.1, num_clusters=1, base_cost=median):**
- Nodes: 3 (KIF22, LDLRAD4, GADD45B)
- Prize: 2.908
- Cost: 0.729 (1 edge)
- Net Benefit: **2.179**

**Greedy Algorithm Solution (starting from top node CEACAM19):**
- Nodes: 21 
- Prize: 122.945
- Cost: 14.579 (20 edges)
- Net Benefit: **108.366**

**Improvement: 4873% better with greedy!**

### 2. ALL Top Nodes Are Profitable to Add

Every single node in the top 20 (by prize) would be profitable to connect to the current solution:

| Node | Prize | Distance | Cost | Net Benefit | Status |
|------|-------|----------|------|-------------|--------|
| CEACAM19 | 26.69 | 4 edges | 2.92 | +23.77 | ✓ Highly profitable |
| CR1 | 22.41 | 4 edges | 2.92 | +19.49 | ✓ Highly profitable |
| RELB | 21.10 | 3 edges | 2.19 | +18.92 | ✓ Highly profitable |
| ... | ... | ... | ... | ... | All profitable |

**Result:** 20/20 top nodes have positive net benefit, yet PCST ignores all of them.

### 3. Network Connectivity is NOT the Issue

- All high-prize nodes ARE connected to the current solution (3-5 edges away)
- Shortest paths exist in the network
- Connectivity analysis shows the network is fully connected

### 4. Cost Parameters Are NOT the Issue

Tested various base_cost values with NO change in module size:

| Base Cost | Factor | Nodes | Edges | Status |
|-----------|--------|-------|-------|--------|
| 0.729 | Current (median) | 3 | 1 | Still 3 nodes |
| 0.144 | Q10 (10th percentile) | 3 | 1 | Still 3 nodes |
| 0.073 | 0.1 × median | 3 | 1 | Still 3 nodes |
| 0.036 | 0.05 × median | 3 | 1 | Still 3 nodes |
| **0.007** | **0.01 × median (1% of original)** | **3** | **1** | **Still 3 nodes!** |
| 0.001 | Even more aggressive | 3 | 1 | Still 3 nodes |

**Reduction to 1% of original base_cost has NO effect.** This rules out cost sensitivity.

### 5. FDR is NOT the Issue

**FDR=0.1 (filtered):** 3 nodes selected
**FDR=1.0 (all nodes significant):** Still 3 nodes selected

With FDR=1.0, CEACAM19 (prize 26.69) is fully available but PCST doesn't include it.

### 6. PCST Parameters Make No Difference

Tested different PCST configurations:
- `num_clusters=1` → 3 nodes
- `num_clusters=2` → 3 nodes  
- `num_clusters=5` → 3 nodes
- `pruning='gw'` → 3 nodes
- `pruning='strong'` → 2 nodes (worse)

## Root Cause Analysis

The problem is **NOT:**
- ✗ Cost structure (tested at 1% of original)
- ✗ FDR threshold (tested at FDR=1.0)
- ✗ Network connectivity (all high-value nodes connected)
- ✗ Parameter tuning (tested num_clusters, pruning)
- ✗ Node availability (FDR=1.0 gives access to all nodes)

The problem appears to be **fundamental to PCST itself:**
- PCST is a Prize-Collecting Steiner Tree approximation algorithm
- It uses GW (Goemans-Williamson) rounding by default
- This is an NP-hard problem; solvers use heuristics
- The algorithm appears to be converging to a local optimum

## Key Observation

When we reduced `base_cost` to 0.001 (1/729th of original), PCST *still* selected only 3 nodes. This suggests:

1. **The algorithm is not exploring beyond a local solution**
2. **It's finding a"stable" configuration and not escaping it**
3. **The issue is algorithmic, not parametric**

## Recommendations

1. **Investigate pcst_fast library behavior**
   - Check if there's a convergence criterion we're missing
   - Try running with verbosity to see what the solver is doing
   - Consider alternative PCST solvers

2. **Test alternative algorithms**
   - Compare with other Steiner tree solvers
   - Consider simpler greedy approaches
   - Explore spectral clustering alternatives

3. **Reconsider the problem formulation**
   - Maybe PCST penalty terms need adjustment
   - Consider adding regularization to encourage larger modules
   - Review if edge costs are conceptually correct

4. **Validate the greedy baseline**
   - The greedy algorithm finds 21-node modules with much better quality
   - Consider if greedy with backtracking might be a better practical choice
   - Test greedy against known biological modules

## Files Involved

- `test_cost_to_connect_high_prize.py` - Shows profitable unselected nodes
- `test_greedy_vs_pcst.py` - Comparison of greedy vs PCST
- `test_base_cost_sensitivity.py` - Cost parameter testing
- `test_pcst_with_all_nodes.py` - FDR and parameter testing
- `test_pcst_debug.py` - Direct pcst_fast analysis
- `test_fdr_1_0.py` - FDR threshold testing

## Conclusion

**PCST is solving the mathematical problem correctly but appears to be trapped at a poor local optimum.** The problem is not with how we've parameterized it, but rather a fundamental limitation of the approximation algorithm or how it's being applied to this problem structure.

