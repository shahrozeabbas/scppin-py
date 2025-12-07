Troubleshooting
================

Common issues and their solutions.

pcst-fast Installation Fails
------------------------------

**Problem**: Installation fails with compilation errors when installing ``pcst-fast``.

**Solutions**:

1. **Install build tools** (required for compiling C extensions):

   **On Ubuntu/Debian:**
   .. code-block:: bash

      sudo apt-get install build-essential python3-dev

   **On macOS:**
   .. code-block:: bash

      xcode-select --install

   **On Windows:**
   Install Visual Studio Build Tools from https://visualstudio.microsoft.com/downloads/

2. **Try installing from conda-forge** (if using conda):
   .. code-block:: bash

      conda install -c conda-forge pcst-fast

3. **Install pcst-fast separately** before installing scppin:
   .. code-block:: bash

      pip install pcst-fast
      pip install scppin

No Module Named 'scppin'
-------------------------

**Problem**: Import error ``ModuleNotFoundError: No module named 'scppin'``.

**Solutions**:

1. **Verify you're in the correct virtual environment**:
   .. code-block:: bash

      which python  # Should point to your venv
      python -m pip list | grep scppin  # Should show scppin installed

2. **Reinstall the package in development mode**:
   .. code-block:: bash

      cd /path/to/scppin-py
      pip install -e .

3. **Check Python version** (requires Python 3.10+):
   .. code-block:: bash

      python --version  # Should be 3.10 or higher

Network File Not Found
-----------------------

**Problem**: FileNotFoundError when loading network files.

**Solutions**:

1. **Check file path**: Ensure the file path is correct relative to your current 
   working directory or use an absolute path.

2. **Copy network files from R package** (if using included networks):
   .. code-block:: bash

      cp ../R/inst/extdata/*.graphml scppin/data/networks/

3. **Verify file format**: Supported formats include CSV, TXT, GraphML, GML, 
   or a list of edge tuples.

Import Errors with Scanpy
--------------------------

**Problem**: ``ImportError`` when trying to use scanpy integration.

**Solution**: Install optional scanpy dependencies:
.. code-block:: bash

   pip install "scppin[scanpy]"

This installs both ``scanpy`` and ``anndata`` packages.

pcst-fast Not Found
--------------------

**Problem**: Runtime error indicating ``pcst-fast`` is not available.

**Solution**: Install pcst-fast separately:
.. code-block:: bash

   pip install pcst-fast

This should install automatically with scppin, but sometimes needs to be 
installed separately if the build fails.

General Import Errors
---------------------

**Problem**: Various ImportError messages.

**Solutions**:

1. **Reinstall the package**:
   .. code-block:: bash

      pip install -e .

2. **Check all dependencies are installed**:
   .. code-block:: bash

      pip install numpy scipy igraph matplotlib pandas pcst-fast

3. **Verify Python version** (requires 3.10+):
   .. code-block:: bash

      python --version

Empty Module Detected
---------------------

**Problem**: ``detect_module()`` runs but returns an empty module.

**Possible causes**:

1. **No genes match p-values**: The network may not contain genes matching 
   the gene names in your p-value dictionary. Check that gene names are 
   correctly normalized (typically uppercase).

2. **FDR too strict**: Try increasing the FDR threshold:
   .. code-block:: python

      model.detect_module(fdr=0.05)  # Less stringent

3. **Network disconnected**: If your network has disconnected components, 
   the algorithm may not find a connecting path.

4. **All p-values too large**: Ensure some genes have small p-values 
   (< 0.05) for module detection.

Memory Issues with Large Networks
----------------------------------

**Problem**: Out of memory errors when working with very large networks.

**Solutions**:

1. **Filter network before loading**: Only include interactions you're 
   interested in.

2. **Use edge filtering**: Many network databases allow filtering by 
   confidence or interaction type.

3. **Work with subnetworks**: Consider analyzing subnetworks of interest 
   rather than the entire PPI network.

4. **Increase system memory**: If possible, use a machine with more RAM.

Platform-Specific Issues
------------------------

macOS (Apple Silicon M1/M2)
~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Some dependencies may need Rosetta. Try:
  .. code-block:: bash

     arch -x86_64 pip install scppin

* Use native ARM Python if possible for better performance.

Windows
~~~~~~~

* Visual Studio Build Tools are required for compiling C extensions.
* Consider using Anaconda for easier dependency management:
  .. code-block:: bash

     conda create -n scppin python=3.10
     conda activate scppin
     pip install scppin

Linux
~~~~~

* Ensure ``python3-dev`` is installed:
  .. code-block:: bash

     sudo apt-get install python3-dev  # Ubuntu/Debian

Getting Help
------------

If you encounter issues not covered here:

* **Check GitHub Issues**: Search existing issues or create a new one at 
  https://github.com/shahrozeabbas/scppin-py/issues
* **Original R Package**: The original R implementation documentation may 
  have additional insights: https://github.com/floklimm/scPPIN
* **Check Dependencies**: Ensure all dependencies are properly installed 
  and up to date

