Tutorials
=========

This section provides step-by-step tutorials for using ADCIRC Utils.

Getting Started
--------------

This tutorial will guide you through the basic setup and usage of ADCIRC Utils.

1. Installation
~~~~~~~~~~~~~~

First, install ADCIRC Utils:

.. code-block:: bash

    pip install adcircutils

2. Basic Usage
~~~~~~~~~~~~~

Import the package and load a grid:

.. code-block:: python

    from adcircutils import Grid
    grid = Grid.from_file('fort.14')

3. Visualization
~~~~~~~~~~~~~~

Plot the grid:

.. code-block:: python

    import matplotlib.pyplot as plt
    grid.plot()
    plt.show()

Advanced Topics
--------------

This tutorial covers more advanced features of ADCIRC Utils.

1. Custom Grid Creation
~~~~~~~~~~~~~~~~~~~~~

Create a custom grid:

.. code-block:: python

    from adcircutils import Grid
    import numpy as np
    
    x = np.linspace(0, 100, 11)
    y = np.linspace(0, 100, 11)
    xx, yy = np.meshgrid(x, y)
    
    grid = Grid.from_arrays(xx, yy)
    grid.write('custom_grid.14')

2. Output Processing
~~~~~~~~~~~~~~~~~~

Process output files:

.. code-block:: python

    from adcircutils import Output
    output = Output.from_file('fort.63.nc')
    water_levels = output.get_water_levels() 