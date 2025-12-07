Parameters Reference
====================

This page documents all parameters for scPPIN-py methods with guidance on when and how to use them.

detect_module() Parameters
---------------------------

The ``detect_module()`` method has several parameters that control module detection:

.. list-table:: detect_module() Parameters
   :header-rows: 1
   :widths: 20 10 70

   * - Parameter
     - Default
     - Description
   * - ``fdr``
     - 0.01
     - False discovery rate threshold. Lower values are more stringent and produce smaller modules. Increase if module is too small, decrease if too large.
   * - ``edge_weight_attr``
     - None
     - Edge attribute name containing weights (e.g., 'weight', 'confidence'). If None, uses uniform edge costs matching the R implementation. Set this when you have confidence scores for interactions.
   * - ``c0``
     - 0.01
     - Minimum edge cost to prevent numerical instability from zero-cost edges. Usually doesn't need to be changed.
   * - ``normalization``
     - 'minmax'
     - Normalization method for edge weights: 'minmax' (scale to [0,1]), 'log1p' (log transform), 'power' (raise to power 6), or None (use weights directly). Only used when ``edge_weight_attr`` is provided.
   * - ``simplify``
     - True
     - Remove self-loops and parallel edges from network. Usually keep as True unless you need to preserve all edges.
   * - ``validate``
     - True
     - Validate network structure (check for empty network, disconnected components, etc.). Keep as True unless you're certain your network is valid.
   * - ``use_max_prize_root``
     - False
     - If True, use the node with highest prize as the PCST root. This makes results more deterministic but may miss optimal solutions. Use False for standard behavior.

When to Adjust Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~

**Module too small?**
* Increase ``fdr`` (e.g., 0.05) to include more genes
* Set ``edge_weight_attr`` if you have confidence scores
* Try ``use_max_prize_root=True`` for different root selection

**Module too large?**
* Decrease ``fdr`` (e.g., 0.001) for stricter filtering
* Use edge weights to prioritize high-confidence interactions

**Have edge confidence scores?**
* Set ``edge_weight_attr='weight'`` (or your attribute name)
* Choose appropriate ``normalization`` method:
  * ``'minmax'``: Default, scales weights to [0, 1]
  * ``'log1p'``: Good for skewed distributions
  * ``'power'``: Emphasizes strong edges (raises to power 6)
  * ``None``: Use weights directly (must already be in [0, 1])

load_network() Parameters
-------------------------

.. list-table:: load_network() Parameters
   :header-rows: 1
   :widths: 20 10 70

   * - Parameter
     - Default
     - Description
   * - ``source``
     - -
     - Network source: file path (str), list of edge tuples, pandas DataFrame, or igraph Graph object
   * - ``weight_column``
     - None
     - Column name in CSV/DataFrame to use as edge weights. If provided, loads weights from that column and sets as 'weight' attribute on edges.
   * - ``fmt``
     - 'auto'
     - File format hint: 'auto' (detect from extension), 'csv', 'graphml', or 'gml'. Usually 'auto' works fine.

Examples:

.. code-block:: python

   # Load from CSV with edge weights
   model.load_network('edges.csv', weight_column='confidence')

   # Load from list
   edges = [('GENE1', 'GENE2'), ('GENE2', 'GENE3')]
   model.load_network(edges)

   # Load GraphML file
   model.load_network('network.graphml', fmt='graphml')

set_node_weights() Parameters
-----------------------------

.. list-table:: set_node_weights() Parameters
   :header-rows: 1
   :widths: 20 10 70

   * - Parameter
     - Default
     - Description
   * - ``weights``
     - -
     - Dictionary mapping gene names to p-values, or AnnData object with rank_genes_groups results
   * - ``groupby``
     - None
     - Key in adata.obs for grouping labels (required when weights is AnnData). Example: 'leiden', 'louvain'
   * - ``group``
     - None
     - Specific group/cluster to extract p-values for (required when weights is AnnData). Example: '0', '1', 'cluster_A'

Examples:

.. code-block:: python

   # From dictionary
   pvalues = {'TP53': 0.0001, 'MDM2': 0.001}
   model.set_node_weights(pvalues)

   # From AnnData (requires scanpy)
   sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')
   model.set_node_weights(adata, groupby='leiden', group='0')

set_edge_weights() Parameters
----------------------------

.. list-table:: set_edge_weights() Parameters
   :header-rows: 1
   :widths: 20 10 70

   * - Parameter
     - Default
     - Description
   * - ``weights``
     - -
     - Dictionary mapping (u, v) tuples to weight values. Only edges present in the network are set (automatically filtered).
   * - ``attr_name``
     - 'weight'
     - Edge attribute name to store weights. Use this name in ``detect_module(edge_weight_attr='weight')``.

Example:

.. code-block:: python

   weights = {
       ('TP53', 'MDM2'): 0.9,
       ('TP53', 'CDKN1A'): 0.8,
   }
   model.set_edge_weights(weights=weights, attr_name='weight')
   module = model.detect_module(fdr=0.01, edge_weight_attr='weight')

plot_module() Parameters
------------------------

.. list-table:: plot_module() Parameters
   :header-rows: 1
   :widths: 20 10 70

   * - Parameter
     - Default
     - Description
   * - ``fdr``
     - 0.01
     - FDR threshold for visualization (should match the FDR used in ``detect_module()``)
   * - ``**kwargs``
     - -
     - Additional plotting arguments passed to matplotlib

Example:

.. code-block:: python

   model.plot_module(fdr=0.01, figsize=(10, 8))

See Also
--------

* :doc:`tutorials/basic_usage` for step-by-step usage examples
* :doc:`api/index` for complete API reference
* :doc:`algorithm/overview` for algorithm details
