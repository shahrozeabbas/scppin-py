Algorithm Overview
==================

scPPIN-py detects functional modules in protein-protein interaction networks by 
integrating single-cell RNA sequencing (scRNA-seq) data using a Prize-Collecting 
Steiner Tree (PCST) approach.

Pipeline Overview
------------------

The algorithm consists of five main steps:

.. list-table::
   :header-rows: 1
   :widths: 10 90

   * - Step
     - Description
   * - **1. Fit BUM Model**
     - Distinguish true signals from noise in p-value distributions using a 
       Beta-Uniform Mixture (BUM) model
   * - **2. Compute Node Scores**
     - Transform p-values to node weights (prizes) using the fitted BUM model
   * - **3. Prepare Edge Costs**
     - Compute edge costs from network topology, optionally incorporating edge 
       weights (confidence scores)
   * - **4. Solve PCST**
     - Find optimal subgraph that maximizes prize (node scores) minus cost 
       (edge costs)
   * - **5. Return Module**
     - Extract and return the connected functional module

Step 1: Fit BUM Model
----------------------

The Beta-Uniform Mixture (BUM) model is used to separate true signals from noise 
in p-value distributions from differential expression analysis. The model 
assumes p-values come from a mixture of:

* **Uniform component** (noise): P-values uniformly distributed in [0, 1]
* **Beta component** (signal): P-values following a Beta distribution

The model is parameterized by:

* :math:`\lambda`: Mixing parameter (proportion of uniform/noise component)
* :math:`\alpha`: Beta distribution shape parameter (:math:`0 < \alpha < 1`)

Smaller :math:`\alpha` values indicate stronger signal enrichment in small p-values.

Step 2: Compute Node Scores
----------------------------

Node scores (prizes) are computed from p-values using the fitted BUM model:

.. math::
   \text{score}(p) = -\log_{10}(p)

Where :math:`p` is the p-value. Nodes with smaller p-values receive higher 
scores (prizes), making them more likely to be included in the detected module.

The network is automatically filtered to only include genes with available 
p-values before module detection.

Step 3: Prepare Edge Costs
---------------------------

Edge costs determine the penalty for including an edge in the module. The default 
cost is based on the network topology, but can be adjusted using edge weights 
(confidence scores):

.. math::
   \text{cost}(e) = \text{base\_cost} \cdot (1 - w(e) \cdot \text{scale}) + c_0

Where:

* :math:`w(e)`: Edge weight (in [0, 1]), higher values indicate more confidence
* :math:`\text{scale}`: Scaling factor for weight influence (default: 1.0)
* :math:`c_0`: Minimum cost to prevent zero-cost edges (default: 0.1 · base_cost)

This formula ensures that:

* Higher edge weights → Lower costs → More likely to include edge
* Minimum cost :math:`c_0` prevents numerical instability from zero-cost edges

If no edge weights are provided, uniform costs are used.

Step 4: Solve PCST
------------------

The Prize-Collecting Steiner Tree (PCST) problem finds a connected subgraph that 
maximizes:

.. math::
   \text{objective} = \sum_{v \in V'} \text{prize}(v) - \sum_{e \in E'} \text{cost}(e)

Where:

* :math:`V'` is the set of nodes in the module
* :math:`E'` is the set of edges in the module
* :math:`\text{prize}(v)` is the node score from Step 2
* :math:`\text{cost}(e)` is the edge cost from Step 3

The solution is guaranteed to be a connected subgraph (tree) that balances 
including high-scoring nodes while minimizing connectivity costs.

Step 5: Return Module
---------------------

The detected module is returned as an igraph Graph containing:

* **Nodes**: Genes significantly contributing to the module
* **Edges**: Protein-protein interactions connecting module genes
* **Node attributes**: Original p-values, computed scores, BUM parameters

False Discovery Rate (FDR) Control
------------------------------------

The FDR threshold controls which nodes are considered "significant" for module 
detection:

* **Lower FDR** (e.g., 0.001): Only genes with very small p-values are included 
  → Smaller, more stringent modules
* **Higher FDR** (e.g., 0.05): Genes with moderate p-values can be included 
  → Larger, less stringent modules

Note: FDR guarantees are approximate when using edge weights, as the cost formula 
adjusts edge costs dynamically.

Original Method
---------------

This implementation follows the method described in:

**Klimm et al. (2020)**: Functional module detection through integration of 
single-cell RNA sequencing data with protein–protein interaction networks. 
BMC Genomics 21, Article number: 756.

`DOI: 10.1186/s12864-020-07144-2 <https://doi.org/10.1186/s12864-020-07144-2>`_

Key Differences from R Implementation
---------------------------------------

* **Class-based API**: Object-oriented design with setup method chaining
* **Edge weight formula**: Uses the author's recommended formula from 
  `GitHub Issue #10 <https://github.com/floklimm/scPPIN/issues/10>`_
* **Automatic filtering**: Network automatically filtered to genes with p-values
* **Performance**: Vectorized NumPy operations provide ~5-10x speedup

Next Steps
----------

* See :doc:`../tutorials/basic_usage` for step-by-step usage guide
* Check :doc:`../api/index` for implementation details
* Explore :doc:`../examples/basic` for complete working examples
