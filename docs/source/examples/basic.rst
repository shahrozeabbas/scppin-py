Basic Example
==============

This example demonstrates the complete workflow for detecting functional modules 
in protein-protein interaction networks.

The example uses a small network based on Alzheimer's disease-related genes to 
illustrate all key features of scPPIN-py.

Complete Code
-------------

.. code-block:: python

   """
   Basic Tutorial for scPPIN
   
   This example demonstrates the basic usage of scPPIN for detecting
   functional modules in protein-protein interaction networks using
   the class-based API.
   """

   import numpy as np
   import matplotlib.pyplot as plt
   from scppin import scPPIN

   # Create analyzer instance
   analyzer = scPPIN()
   
   # Step 1: Load network
   print("1. Loading network...")
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
   print(f"   Network: {analyzer.network.number_of_nodes()} nodes, "
         f"{analyzer.network.number_of_edges()} edges")
   
   # Step 2: Set node weights (p-values)
   print("\n2. Setting node weights (p-values)...")
   pvalues = {
       # Significant genes (cluster 1)
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
   print(f"   Node weights for {len(analyzer.node_weights)} genes")
   
   # Step 3: Detect functional module
   print("\n3. Detecting functional module...")
   fdr = 0.01
   
   module = analyzer.detect_module(fdr=fdr)
   
   print(f"\nModule detected!")
   print(f"  Nodes: {module.number_of_nodes()}")
   print(f"  Edges: {module.number_of_edges()}")
   print(f"  Genes in module: {list(module.nodes())}")
   
   # Print node scores
   print("\nNode scores:")
   for node in module.nodes():
       score = module.nodes[node].get('score', 'N/A')
       print(f"  {node}: {score:.4f}" if isinstance(score, float) else f"  {node}: {score}")
   
   # Step 4: Visualize (optional)
   print("\n4. Visualizing module...")
   
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
   plt.savefig('basic_tutorial_result.png', dpi=150, bbox_inches='tight')
   print("Visualization saved to 'basic_tutorial_result.png'")

Expected Output
---------------

When you run this example, you should see output similar to:

.. code-block:: text

   1. Loading network...
      Network: 12 nodes, 12 edges

   2. Setting node weights (p-values)...
      Node weights for 12 genes

   3. Detecting functional module...

   Module detected!
     Nodes: 6
     Edges: 5
     Genes in module: ['APP', 'APOE', 'PSEN1', 'PSEN2', 'CLU', 'ABCA1']

   Node scores:
     APP: 4.0000
     APOE: 3.3010
     PSEN1: 3.0000
     PSEN2: 2.3010
     CLU: 2.0000
     ABCA1: 1.6990

   4. Visualizing module...
   Visualization saved to 'basic_tutorial_result.png'

Explanation
-----------

1. **Network Loading**: The network contains 12 genes connected by 12 edges. 
   After setting node weights, the network may be filtered automatically.

2. **Node Weights**: P-values represent differential expression significance. 
   Smaller values (e.g., 0.0001) indicate highly significant genes.

3. **Module Detection**: With FDR=0.01, the algorithm detects a module with 
   6 nodes and 5 edges. These are the most significant genes that form a 
   connected subgraph.

4. **Visualization**: The module visualization shows the detected subgraph 
   colored by p-value significance (darker = more significant).

Example with Edge Weights
--------------------------

Here's an extension that includes edge weights:

.. code-block:: python

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
   module_with_weights = analyzer2.detect_module(fdr=fdr, edge_weight_attr='weight')
   
   print(f"\nModule with edge weights:")
   print(f"  Nodes: {module_with_weights.number_of_nodes()}")
   print(f"  Edges: {module_with_weights.number_of_edges()}")
   print(f"  Genes: {list(module_with_weights.nodes())}")

Edge weights prioritize high-confidence interactions, potentially changing 
which genes are included in the module. Set ``edge_weight_attr`` to the name 
of the edge attribute that stores weights (usually ``'weight'``).

Running the Example
--------------------

This example is available as ``examples/basic_tutorial.py`` in the repository.

To run it:

.. code-block:: bash

   cd examples
   python basic_tutorial.py

Make sure you have all dependencies installed:
- scppin (with pcst-fast)
- matplotlib
- numpy

See Also
--------

* :doc:`../tutorials/basic_usage` for a detailed walkthrough
* :doc:`../quickstart` for a minimal example
* :doc:`../api/index` for full API documentation
