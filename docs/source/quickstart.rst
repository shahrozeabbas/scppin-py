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

   # Create model instance
   model = scPPIN()

   # 1. Load network from edge list
   model.load_network('edges.csv')

   # 2. Set node weights (p-values)
   pvalues = {'TP53': 0.0001, 'MDM2': 0.001, 'CDKN1A': 0.005}
   model.set_node_weights(pvalues)

   # 3. Detect functional module
   module = model.detect_module(fdr=0.01)

   # 4. Visualize
   model.plot_module(fdr=0.01)

Method Chaining
---------------

Setup methods return ``self``, so you can chain operations:

.. code-block:: python

   from scppin import scPPIN

   model = (scPPIN()
              .load_network('edges.csv')
              .set_node_weights({'TP53': 0.0001, 'MDM2': 0.001}))
   module = model.detect_module(fdr=0.01)

   # Access results
   print(f"Module has {module.vcount()} nodes")

Setup methods (``load_network()``, ``set_node_weights()``, ``set_edge_weights()``) 
return ``self`` so you can chain them. ``detect_module()`` returns the resulting 
igraph Graph while also storing it on ``model.module``.

What to Expect
--------------

After running ``detect_module()``, you'll have:

* ``model.module``: An igraph Graph containing the detected functional module
* ``model.bum_params``: BUM model parameters (lambda, alpha)
* ``model.node_scores``: Computed node scores

The module is a connected subgraph that maximizes the sum of node prizes (from p-values) 
minus the sum of edge costs (from network topology and optional edge weights).

Next Steps
----------

* Read the :doc:`tutorials/basic_usage` tutorial for a detailed walkthrough
* Explore :doc:`examples/basic` for a complete working example
* Check out the :doc:`api/index` for full API documentation
