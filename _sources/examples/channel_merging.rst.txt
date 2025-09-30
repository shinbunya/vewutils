Channel Mesh Merging
=====================

This example demonstrates the process of merging multiple meshes (channel, non-channel, and background) while handling Vertical Element Wall (VEW). The workflow includes:

* Combining channel and non-channel meshes with VEWs
* Subtracting channel+non-channel coverage from background mesh
* Merging all components into a final mesh
* Adjusting VEW channel elevations and barrier heights
* Transferring nodal attributes and updating Manning's n values

You can find the complete example in the `examples/channelmerging/example.ipynb <https://github.com/shinbunya/vewutils/tree/main/examples/channelmerging/example.ipynb>`_ notebook.

.. image:: /_static/channel_merging.png
   :alt: Channel Merging Example
   :align: center
