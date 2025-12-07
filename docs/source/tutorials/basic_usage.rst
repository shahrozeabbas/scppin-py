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

   # Create model instance
   model = scPPIN()

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

   model.load_network(edges)
   print(f"Network: {model.network.vcount()} nodes, "
         f"{model.network.ecount()} edges")

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

   model.set_node_weights(pvalues)
   print(f"Node weights for {len(model.node_weights)} genes")

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
   module = model.detect_module(fdr=fdr)

   print(f"\nModule detected!")
   print(f"  Nodes: {module.vcount()}")
   print(f"  Edges: {module.ecount()}")
   print(f"  Genes in module: {list(module.vs['name'])}")

The FDR parameter controls the stringency:
* **Lower FDR** (e.g., 0.001): More stringent, smaller modules
* **Higher FDR** (e.g., 0.05): Less stringent, larger modules

After detection, you can access:

* ``model.module``: igraph Graph containing the detected module
* ``model.bum_params``: BUM model parameters (lambda, alpha)  
* ``model.node_scores``: Computed node scores

Node Scores
~~~~~~~~~~~

Each node in the module has a score attribute:

.. code-block:: python

   print("\nNode scores:")
   for v in module.vs:
       node_name = v['name']
       score = v.get('score', 'N/A')
       print(f"  {node_name}: {score:.4f}" if isinstance(score, float) else f"  {node_name}: {score}")

Higher scores indicate genes with stronger signal (smaller p-values).

Visualization
-------------

Visualize the detected module:

.. code-block:: python

   import matplotlib.pyplot as plt

   # Plot detected module
   model.plot_module(fdr=fdr)
   
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

   # Create new model for edge weights example
   model2 = scPPIN()
   model2.load_network(edges)
   model2.set_node_weights(pvalues)
   
   # Add confidence scores to edges
   np.random.seed(42)
   weights = {}
   for u, v in edges:
       # Higher confidence for edges between significant genes
       if pvalues.get(u, 1.0) < 0.01 and pvalues.get(v, 1.0) < 0.01:
           weights[(u, v)] = np.random.uniform(0.8, 1.0)
       else:
           weights[(u, v)] = np.random.uniform(0.3, 0.7)
   
   model2.set_edge_weights(weights=weights)
   module_with_weights = model2.detect_module(fdr=fdr, edge_weight_attr='weight')
   
   print(f"\nModule with edge weights:")
   print(f"  Nodes: {module_with_weights.vcount()}")
   print(f"  Edges: {module_with_weights.ecount()}")
   print(f"  Genes: {list(module_with_weights.vs['name'])}")

Set ``edge_weight_attr`` to the edge attribute that stores weight values (``'weight'`` by default).

Using Scanpy/AnnData
--------------------

scPPIN-py integrates seamlessly with scanpy and AnnData objects. You can extract 
p-values directly from differential expression results:

.. code-block:: python

   import scanpy as sc
   from scppin import scPPIN

   # Load your single-cell data
   adata = sc.read_h5ad('your_data.h5ad')

   # Preprocess (if not already done)
   sc.pp.normalize_total(adata, target_sum=1e4)
   sc.pp.log1p(adata)

   # Cluster cells
   sc.pp.neighbors(adata)
   sc.tl.leiden(adata, key_added='leiden')

   # Run differential expression to get marker genes per cluster
   sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon')

   # Create scPPIN model
   model = scPPIN()
   model.load_network('ppi_network.csv')

   # Extract p-values for a specific cluster (e.g., cluster '0')
   model.set_node_weights(adata, groupby='leiden', group='0')

   # Detect functional module
   module = model.detect_module(fdr=0.01)

   print(f'Module: {module.vcount()} genes, {module.ecount()} edges')
   print(f'Genes: {module.vs["name"]}')

The ``set_node_weights()`` method automatically extracts p-values from 
``adata.uns['rank_genes_groups']`` for the specified group. Make sure you've run 
``sc.tl.rank_genes_groups()`` first.

Analyzing Multiple Clusters
~~~~~~~~~~~~~~~~~~~~~~~~~~~

You can analyze multiple clusters in a loop:

.. code-block:: python

   # Analyze each cluster
   for cluster in adata.obs['leiden'].unique():
       model = scPPIN()
       model.load_network('ppi_network.csv')
       model.set_node_weights(adata, groupby='leiden', group=cluster)
       
       module = model.detect_module(fdr=0.01)
       print(f'Cluster {cluster}: {module.vcount()} genes in module')

Key Takeaways
-------------

* Networks are automatically normalized and filtered to genes with p-values
* Smaller p-values → Higher node prizes → More likely to be in module
* FDR threshold controls module size (lower = smaller, more stringent)
* Edge weights can improve module detection by prioritizing high-confidence 
  interactions
* Setup methods (load/set operations) support method chaining for concise code

Next Steps
----------

* See :doc:`../examples/basic` for a complete runnable example
* Explore :doc:`../api/index` for full API documentation
* Learn about the algorithm in :doc:`../algorithm/overview`
