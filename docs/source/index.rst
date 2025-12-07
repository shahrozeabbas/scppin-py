Welcome to scPPIN-py
======================

.. raw:: html

   <div style="font-size: 1.2em; margin-bottom: 2em;">
      Pure Python implementation of scPPIN with class-based API
   </div>

**scPPIN-py** detects functional modules in protein-protein interaction networks by integrating single-cell RNA sequencing data. This is a reimplementation of the original R package with a clean, object-oriented API.

.. note::
   **Original Method**: The method is described in `Klimm et al. (2020) <https://doi.org/10.1186/s12864-020-07144-2>`_ (BMC Genomics 21, Article number: 756).

Key Features
------------

* **Class-Based API**: Object-oriented design with setup method chaining
* **Pure Python**: No external C++ binaries (uses ``pcst_fast``)
* **Edge Weights**: Supports confidence scores or computed correlations
* **Author's Formula**: Implements recommended edge cost formula from `GitHub Issue #10 <https://github.com/floklimm/scPPIN/issues/10>`_
* **Scanpy Integration**: Works seamlessly with AnnData objects
* **Fast**: Vectorized NumPy operations, ~5-10x faster than R
* **Easy Installation**: Single ``pip install`` command
* **Automatic Filtering**: Network automatically filtered to genes with p-values

Quick Start
-----------

Install scPPIN-py from GitHub:

.. code-block:: bash

   pip install git+https://github.com/shahrozeabbas/scppin-py.git

Quick Example
~~~~~~~~~~~~~

.. code-block:: python

   from scppin import scPPIN

   # Create model instance
   model = scPPIN()

   # Load network from edge list
   model.load_network('edges.csv')

   # Set node weights (p-values)
   pvalues = {'TP53': 0.0001, 'MDM2': 0.001, 'CDKN1A': 0.005}
   model.set_node_weights(pvalues)

   # Detect functional module
   module = model.detect_module(fdr=0.01)

   # Visualize
   model.plot_module(fdr=0.01)

See the :doc:`quickstart` guide for more details.

Documentation Contents
-----------------------

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   installation
   quickstart

.. toctree::
   :maxdepth: 2
   :caption: User Guide

   tutorials/index
   parameters

.. toctree::
   :maxdepth: 2
   :caption: API Reference

   api/index

.. toctree::
   :maxdepth: 2
   :caption: Algorithm

   algorithm/index

.. toctree::
   :maxdepth: 2
   :caption: Examples

   examples/index

.. toctree::
   :maxdepth: 1
   :caption: Additional Resources

   troubleshooting

Citation
--------

If you use this package, please cite the original paper:

.. code-block:: bibtex

   @article{klimm2020functional,
     title={Functional module detection through integration of single-cell RNA 
            sequencing data with protein--protein interaction networks},
     author={Klimm, Florian and Toledo, Enrique M and Monfeuga, Thomas and 
             Zhang, Fang and Deane, Charlotte M and Reinert, Gesine},
     journal={BMC genomics},
     volume={21},
     number={1},
     pages={1--10},
     year={2020},
     publisher={BioMed Central}
   }

Links
-----

* `Original R Package <https://github.com/floklimm/scPPIN>`_
* `GitHub Repository <https://github.com/shahrozeabbas/scppin-py>`_
* `Original Paper <https://doi.org/10.1186/s12864-020-07144-2>`_

License
-------

GPL-3.0 (same as original R package)
