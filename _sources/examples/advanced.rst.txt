Advanced Examples
===============

This section provides advanced examples of how to use ADCIRC Utils.

Creating a Custom Grid
--------------------

.. code-block:: python

    from adcircutils import Grid
    import numpy as np
    
    # Create a simple rectangular grid
    x = np.linspace(0, 100, 11)
    y = np.linspace(0, 100, 11)
    xx, yy = np.meshgrid(x, y)
    
    # Create grid object
    grid = Grid.from_arrays(xx, yy)
    grid.write('custom_grid.14')

Processing Multiple Output Files
-----------------------------

.. code-block:: python

    from adcircutils import Output
    import glob
    
    # Process all output files in a directory
    for file in glob.glob('output/*.63.nc'):
        output = Output.from_file(file)
        max_water_level = output.get_water_levels().max()
        print(f'Max water level in {file}: {max_water_level}')

Custom Boundary Conditions
------------------------

.. code-block:: python

    from adcircutils import Boundary
    import numpy as np
    
    # Create a custom boundary condition
    time = np.arange(0, 24, 0.5)
    water_level = np.sin(time * 2 * np.pi / 12)  # 12-hour tidal cycle
    
    boundary = Boundary.from_arrays(time, water_level)
    boundary.write('tidal_boundary.19') 