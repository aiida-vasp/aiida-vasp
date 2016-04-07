General Purpose Plugins
-----------------------

These are calculation plugins which aim to expose the full flexibility
of the associated programs and will retrieve all output files, making
no assumption about which output will be needed.
This comes at the cost that potentially big files that won't be used
in analysis will be retrieved and permanently stored in your data base.
Therefore it is recommended to use these only where no special purpose
calculation is available for the task.
 

.. toctree::
   :maxdepth: 4
   :caption: Plugins

   vasp-5.X.X - Universal VASP 5 <vasp5>
   wannier90.x <wannier>
