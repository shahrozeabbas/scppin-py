API Reference
==============

This section contains the complete API reference for scPPIN-py.

Main API
--------

.. toctree::
   :maxdepth: 2

   scppin

The main API consists of the ``scPPIN`` class, which provides methods for:

* Loading protein-protein interaction networks
* Setting node weights (p-values) and edge weights
* Detecting functional modules
* Visualizing results

Setup methods (``load_network()``, ``set_node_weights()``, ``set_edge_weights()``) 
support method chaining for convenient workflows. ``detect_module()`` returns 
the resulting igraph Graph while also storing it on ``model.module``.
