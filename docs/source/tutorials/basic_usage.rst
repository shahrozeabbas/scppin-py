Basic Usage Tutorial
=====================

This tutorial demonstrates how to use scPPIN-py to detect functional modules 
in protein-protein interaction networks.

Overview
--------

In this tutorial, we'll:

#. Load a protein-protein interaction network
#. Set node weights (p-values from differential expression)
#. Detect a functional module
#. Visualize the results
#. (Optional) Use edge weights for improved detection

Loading Networks
----------------

Networks can be loaded from various sources. Here, we'll create a small example 
network using Alzheimer's disease-related genes:

.. code-block:: python

   from scppin import scPPIN

   # Create analyzer instance
   analyzer = scPPIN()

   # Load network from edge list
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
   print(f"Network: {analyzer.network.number_of_nodes()} nodes, "
         f"{analyzer.network.number_of_edges()} edges")

After loading, the network is automatically normalized (gene names converted to 
uppercase) and filtered to only include genes for which you provide p-values 
(see next step).

Setting Node Weights
--------------------

Node weights are p-values from differential expression analysis. Smaller p-values 
indicate more significant differential expression:

.. code-block:: python

   pvalues = {
       # Highly significant genes (cluster 1)
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
   print(f"Node weights for {len(analyzer.node_weights)} genes")

When you call ``set_node_weights()``, the network is automatically filtered to 
only include genes that have p-values. This ensures that module detection only 
considers genes with expression data.

Understanding P-values
~~~~~~~~~~~~~~~~~~~~~~~

* **Very small p-values (< 0.01)**: Highly significant genes with strong 
  differential expression
* **Small p-values (0.01-0.05)**: Moderately significant genes
* **Large p-values (> 0.05)**: Genes without significant differential expression

Genes with smaller p-values will have higher "prizes" in the module detection 
algorithm, making them more likely to be included in the detected module.

Detecting Modules
-----------------

Use ``detect_module()`` to find the functional module:

.. code-block:: python

   fdr = 0.01  # False discovery rate threshold
   analyzer.detect_module(fdr=fdr)

   print(f"\nModule detected!")
   print(f"  Nodes: {analyzer.module.number_of_nodes()}")
   print(f"  Edges: {analyzer.module.number_of_edges()}")
   print(f"  Genes in module: {list(analyzer.module.nodes())}")

The FDR parameter controls the stringency:
* **Lower FDR** (e.g., 0.001): More stringent, smaller modules
* **Higher FDR** (e.g., 0.05): Less stringent, larger modules

After detection, you can access:

* ``analyzer.module``: NetworkX graph containing the detected module
* ``analyzer.bum_params``: BUM model parameters (lambda, alpha)  
* ``analyzer.node_scores``: Computed node scores

Node Scores
~~~~~~~~~~~

Each node in the module has a score attribute:

.. code-block:: python

   print("\nNode scores:")
   for node in analyzer.module.nodes():
       score = analyzer.module.nodes[node].get('score', 'N/A')
       print(f"  {node}: {score:.4f}" if isinstance(score, float) else f"  {node}: {score}")

Higher scores indicate genes with stronger signal (smaller p-values).

Visualization
-------------

Visualize the detected module:

.. code-block:: python

   import matplotlib.pyplot as plt

   fig, axes = plt.subplots(1, 2, figsize=(14, 6))
   
   # Plot original network
   import networkx as nx
   pos = nx.spring_layout(analyzer.network, seed=42)
   nx.draw(analyzer.network, pos, with_labels=True, node_color='lightblue',
           node_size=500, font_size=8, ax=axes[0])
   axes[0].set_title('Original Network')
   axes[0].axis('off')
   
   # Plot detected module
   analyzer.plot_module(fdr=fdr, ax=axes[1])
   
   plt.tight_layout()
   plt.savefig('module_result.png', dpi=150, bbox_inches='tight')

The ``plot_module()`` method colors nodes by their p-values, with darker colors 
indicating more significant genes.

Using Edge Weights (Optional)
------------------------------

Edge weights represent confidence in protein-protein interactions. Higher weights 
make edges more likely to be included in the module:

.. code-block:: python

   import numpy as np

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
   analyzer2.detect_module(fdr=fdr, edge_weight_scale=0.5)
   
   print(f"\nModule with edge weights:")
   print(f"  Nodes: {analyzer2.module.number_of_nodes()}")
   print(f"  Edges: {analyzer2.module.number_of_edges()}")
   print(f"  Genes: {list(analyzer2.module.nodes())}")

The ``edge_weight_scale`` parameter controls how much edge weights influence the 
cost calculation. See :doc:`../algorithm/overview` for details.

Key Takeaways
-------------

* Networks are automatically normalized and filtered to genes with p-values
* Smaller p-values → Higher node prizes → More likely to be in module
* FDR threshold controls module size (lower = smaller, more stringent)
* Edge weights can improve module detection by prioritizing high-confidence 
  interactions
* All methods support method chaining for concise code

Next Steps
----------

* See :doc:`../examples/basic` for a complete runnable example
* Explore :doc:`../api/index` for full API documentation
* Learn about the algorithm in :doc:`../algorithm/overview`

