Installation
============

Quick Install
-------------

Install scPPIN-py directly from GitHub:

.. code-block:: bash

   pip install git+https://github.com/shahrozeabbas/scppin-py.git

.. note::
   The package is not yet available on PyPI. Install from GitHub using the command above.

Development Install
-------------------

If you want to contribute or build documentation:

1. **Clone the repository:**

.. code-block:: bash

   git clone https://github.com/shahrozeabbas/scppin-py.git
   cd scppin-py

2. **Create a virtual environment (recommended):**

.. code-block:: bash

   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate

   .. note::
      On Windows, use ``venv\Scripts\activate`` instead.

3. **Install in development mode:**

.. code-block:: bash

   pip install -e ".[all]"

This installs the package with all optional dependencies including documentation tools.

Dependencies
------------

Core Dependencies (required)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These are automatically installed with scPPIN-py:

* numpy >= 1.20.0,<2.0
* scipy >= 1.7.0
* igraph >= 0.11.0
* matplotlib >= 3.3.0
* pandas >= 1.3.0
* pcst-fast >= 1.0.0

Optional Dependencies
~~~~~~~~~~~~~~~~~~~~~~

Install these for additional functionality:

* **scanpy**: For scanpy integration and AnnData support
  .. code-block:: bash

     pip install "scppin[scanpy]"

  This installs:
  
  * scanpy >= 1.8.0
  * anndata >= 0.8.0

* **dev**: For development and testing
  .. code-block:: bash

     pip install "scppin[dev]"

  This installs:
  
  * pytest >= 7.0.0
  * pytest-cov >= 3.0.0
  * black >= 22.0.0
  * flake8 >= 4.0.0
  * mypy >= 0.950

* **docs**: For building documentation
  .. code-block:: bash

     pip install "scppin[docs]"

Installing pcst-fast
--------------------

The ``pcst-fast`` package is required for solving the Prize-Collecting Steiner Tree problem. 
It should install automatically, but if you encounter issues:

.. code-block:: bash

   pip install pcst-fast

If compilation fails, you may need to install build tools:

On Ubuntu/Debian
~~~~~~~~~~~~~~~~

.. code-block:: bash

   sudo apt-get install build-essential python3-dev

On macOS
~~~~~~~~

.. code-block:: bash

   xcode-select --install

On Windows
~~~~~~~~~~

Install Visual Studio Build Tools from https://visualstudio.microsoft.com/downloads/

Verifying Installation
----------------------

Run this simple test to verify your installation:

.. code-block:: python

   from scppin import scPPIN

   # Run a simple test
   model = scPPIN()
   model.load_network([('A', 'B'), ('B', 'C')])
   model.set_node_weights({'A': 0.001, 'B': 0.005, 'C': 0.01})
   module = model.detect_module(fdr=0.01)

   print(f"Module detected: {module.vcount()} nodes")

If this runs without errors, installation is successful!

Platform-Specific Notes
-----------------------

macOS (Apple Silicon M1/M2)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

* May need to install Rosetta for some dependencies
* Use native ARM Python if possible

Windows
~~~~~~~

* May need Visual Studio Build Tools
* Use Anaconda for easier dependency management

Linux
~~~~~

* Should work out of the box on most distributions
* May need ``python3-dev`` package

.. code-block:: bash

   sudo apt-get install python3-dev  # Ubuntu/Debian
   sudo yum install python3-devel     # CentOS/RHEL

Next Steps
----------

Now that you have scPPIN-py installed, check out the :doc:`quickstart` guide to get started!
