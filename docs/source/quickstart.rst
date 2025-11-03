Quick Start
===========

This guide will get you up and running with scPPIN-py in 5 minutes.

Basic Workflow
--------------

The typical workflow consists of four steps:

#. Load a protein-protein interaction network
#. Set node weights (p-values from differential expression)
#. Detect functional module
#. Visualize results

Here's a complete example:

.. code-block:: python

   from scppin import scPPIN

   # Create analyzer instance
   analyzer = scPPIN()

   # 1. Load network from edge list
   analyzer.load_network('edges.csv')

   # 2. Set node weights (p-values)
   pvalues = {'TP53': 0.0001, 'MDM2': 0.001, 'CDKN1A': 0.005}
   analyzer.set_node_weights(pvalues)

   # 3. Detect functional module
   analyzer.detect_module(fdr=0.01)

   # 4. Visualize
   analyzer.plot_module(fdr=0.01)

Method Chaining
---------------

All methods return ``self``, so you can chain operations:

.. code-block:: python

   from scppin import scPPIN

   analyzer = (scPPIN()
              .load_network('edges.csv')
              .set_node_weights({'TP53': 0.0001, 'MDM2': 0.001})
              .detect_module(fdr=0.01))

   # Access results
   print(f"Module has {analyzer.module.number_of_nodes()} nodes")

What to Expect
--------------

After running ``detect_module()``, you'll have:

* ``analyzer.module``: A NetworkX graph containing the detected functional module
* ``analyzer.bum_params``: BUM model parameters (lambda, alpha)
* ``analyzer.node_scores``: Computed node scores

The module is a connected subgraph that maximizes the sum of node prizes (from p-values) 
minus the sum of edge costs (from network topology and optional edge weights).

Understanding the Parameters
-----------------------------

* **FDR threshold** (``fdr=0.01``): False discovery rate threshold. Lower values 
  are more stringent and produce smaller modules.
  
* **Edge weights**: Optional confidence scores for interactions. Higher weights 
  make edges more likely to be included in the module.

Loading Networks
----------------

You can load networks from various sources:

.. code-block:: python

   # From CSV file
   analyzer.load_network('edges.csv')

   # From list of tuples
   edges = [('TP53', 'MDM2'), ('TP53', 'CDKN1A')]
   analyzer.load_network(edges)

   # From DataFrame
   import pandas as pd
   df = pd.DataFrame({'source': ['TP53'], 'target': ['MDM2']})
   analyzer.load_network(df)

   # With edge weights from CSV column
   analyzer.load_network('edges.csv', weight_column='confidence')

Setting Node Weights
--------------------

Node weights are p-values from differential expression analysis:

.. code-block:: python

   # From dictionary
   pvalues = {
       'TP53': 0.0001,  # Very small p-value indicates statistical significance
       'MDM2': 0.001,
       'CDKN1A': 0.005,
   }
   analyzer.set_node_weights(pvalues)

   # Note: Very small p-values indicate highly significant genes
   # These will have higher "prizes" in the module detection algorithm

Next Steps
----------

* Read the :doc:`tutorials/basic_usage` tutorial for a detailed walkthrough
* Explore :doc:`examples/basic` for a complete working example
* Check out the :doc:`api/index` for full API documentation

