Installation
============

ADCIRC Utils can be installed using pip or conda.

Using pip
---------

.. code-block:: bash

    pip install adcircutils

Using conda
-----------

.. code-block:: bash

    conda install -c conda-forge adcircutils

Development Installation
-----------------------

To install from source for development:

.. code-block:: bash

    git clone https://github.com/sbunya/adcircutils.git
    cd adcircutils
    pip install -e ".[dev]"

Dependencies
------------

ADCIRC Utils requires Python 3.11 or later and has the following dependencies:

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