Basic Example
==============

Complete runnable example - copy and modify for your data.

.. code-block:: python

   """
   Basic Example for scPPIN
   
   Complete workflow for detecting functional modules in protein-protein
   interaction networks.
   """

   import numpy as np
   import matplotlib.pyplot as plt
   from scppin import scPPIN

   # Create model instance
   model = scPPIN()
   
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
   
   model.load_network(edges)
   print(f"   Network: {model.network.vcount()} nodes, "
         f"{model.network.ecount()} edges")
   
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
   
   model.set_node_weights(pvalues)
   print(f"   Node weights for {len(model.node_weights)} genes")
   
   # Step 3: Detect functional module
   print("\n3. Detecting functional module...")
   fdr = 0.01
   
   module = model.detect_module(fdr=fdr)
   
   print(f"\nModule detected!")
   print(f"  Nodes: {module.vcount()}")
   print(f"  Edges: {module.ecount()}")
   print(f"  Genes in module: {list(module.vs['name'])}")
   
   # Print node scores
   print("\nNode scores:")
   for v in module.vs:
       node_name = v['name']
       score = v.get('score', 'N/A')
       print(f"  {node_name}: {score:.4f}" if isinstance(score, float) else f"  {node_name}: {score}")
   
   # Step 4: Visualize (optional)
   print("\n4. Visualizing module...")
   
   # Plot detected module
   model.plot_module(fdr=fdr)
   
   plt.tight_layout()
   plt.savefig('basic_tutorial_result.png', dpi=150, bbox_inches='tight')
   print("Visualization saved to 'basic_tutorial_result.png'")

See Also
--------

* :doc:`../tutorials/basic_usage` for a detailed walkthrough
* :doc:`../quickstart` for a minimal example
* :doc:`../api/index` for full API documentation
