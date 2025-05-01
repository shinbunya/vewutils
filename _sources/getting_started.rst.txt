Getting Started
===============

Installation
------------

To use these tools, clone the repository and install it by executing ``install.sh``. Install ``conda`` before using the installation script. The installation script creates a new conda environment, the name of which is specified as a command line argument. It installs all dependencies, including a fork of adcircpy. It then installs ``vewutils`` using ``pip`` command.

Run the following commands for installation:

.. code-block:: bash

    git clone https://github.com/shinbunya/vewutils.git vewutils
    cd vewutils
    ./install.sh vewutils # <-- You can specify your own conda environment name instead of vewutils.

Manual Installation
-------------------

If you prefer to install manually, you can use pip:

.. code-block:: bash

    git clone https://github.com/shinbunya/vewutils.git vewutils
    cd vewutils
    pip install .

Or for development:

.. code-block:: bash

    git clone https://github.com/shinbunya/vewutils.git vewutils
    cd vewutils
    pip install -e .

Dependencies
------------

VEW Utils requires and is tested withPython 3.11 and has the following dependencies:

* adcircpy (use a fork from https://github.com/shinbunya/adcircpy)
* pandas
* numpy
* shapely
* plotly
* geopandas
* scipy
* rasterio
* geowombat
* pyproj>=2.6
* yaml
* matplotlib
* xarray
* requests
* pytz
* jupyter
* ipywidgets
* netCDF4
* haversine
* paramiko
* pooch
* psutil
* searvey
* typepigeon<2
* utm
* appdirs

Basic Usage
-----------

After installation, you can import the package in Python:

.. code-block:: python

    import vewutils

Each tool can be accessed as a subpackage:

.. code-block:: python

    from vewutils import channelpaving
    from vewutils import mesh
    from vewutils import vewprocessing
    # etc.

Running Command-Line Tools
--------------------------

Part of the functionality is also available as command-line tools. 
See the :doc:`usage_guides` section for detailed examples of how to use each tool via command line.
