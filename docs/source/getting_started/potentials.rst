=============================
4. Add the potential datasets
=============================

To run a VASP calculation, potentials (the POTCAR files) have to be uploaded to
the database. For more details regarding the handling of the potentials, please see :ref:`potentials`.

Assuming you already have a valid license you can download the `VASP`_ potentials from their portal and save
it to a convenient location. For the example here, we use ``$HOME/myaiida/potpar_PBE.54.tar``. In `AiiDA-VASP`_
we refer to a set of potentials as a potential family. Execute the following command to upload the whole potential family to the database::

  % verdi data vasp-potcar uploadfamily --path=$HOME/myaiida/potpaw_PBE.54.tar --name=PBE.54 --description="PBE potentials version 54"
  POTCAR files found: 327. New files uploaded: 327, Added to Family: 327

We use the POTCAR parser from `pymatgen` to get the metadata and sometimes this issues a warning if it detects unknown metadata flags in the potentials. You can usually ignore these warnings.

The ``name`` and ``description`` are not optional and have to be
specified. The ``path`` could be either an archive, or one could use
a folder name.  It is also possible, not to specify path, where the plugin will
traverse all folders from the folder in which the command above is executed from and try to upload all
potentials it finds to the specified family.

.. _VASP: https://www.vasp.at
.. _AiiDA-VASP: https://github.com/aiida-vasp/aiida-vasp
