Basic Examples
=============

This section provides basic examples of how to use ADCIRC Utils.

Loading a Grid
-------------

.. code-block:: python

    from adcircutils import Grid
    grid = Grid.from_file('fort.14')

Plotting a Grid
--------------

.. code-block:: python

    import matplotlib.pyplot as plt
    from adcircutils import Grid
    
    grid = Grid.from_file('fort.14')
    grid.plot()
    plt.show()

Reading Output Files
------------------

.. code-block:: python

    from adcircutils import Output
    output = Output.from_file('fort.63.nc')
    water_levels = output.get_water_levels() 