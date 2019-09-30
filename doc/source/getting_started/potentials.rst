Upload potential datasets to the AiiDA database
===============================================

To run a VASP calculation, potentials (the POTCAR file content) have to be uploaded to
the database. For more details regarding the handling of the potentials, please see :ref:`potentials`.

If you have a valid license you can download the VASP potentials at their download service
and follow :ref:`production_potentials`. If not, you are probably following a tutorial,
so please follow :ref:`tutorial_potentials`.

.. _tutorial_potentials:
     
Tutorial potentials
-------------------

In case you are following the tutorials, we only have to upload the tutorial potentials, which
are already available on the public VASP wiki. We have conveniently collected them in an archive
for you do download. Let us go trough the procedure of uploading this set of potentials to the database.

1. First, download the file to some convenient location::

  wget https://github.com/aiida-vasp/aiida-vasp/raw/develop/tutorials/vasp_potentials.tar.gz

2. Untar it::
     
   tar xvzf vasp_potentials.tar.gz
   
3. Go into the ``vasp_potentials`` directory. We will now upload the potentials
using the custom commands created for the POTCAR data types::
  
  verdi data vasp-potcar uploadfamily --name pbe --description "A few tutorial PBE potentials"

4. This will upload the potentials into the database and hash them. E.g. you will not be able to have multiple entries of one potential in the database. Also, in the calculations, only the hash is used such that the POTCAR data (which is covered by license) is not revealed. This should complete with a message that three potentials was found and uploaded.

.. _production_potentials:
   
Production ready potentials
---------------------------
   
Usually, most VASP users have a potential set (or several) they like to use. Here we assume this is
the modified PBE potential that are GW ready and was supplied with VASP.5.4.4. First, make sure
you have the archive, ``potpaw_PBE.54.tar`` located somewhere. Here we assume it is placed in
here ``$HOME/myaiida/potpaw_PBE.54.tar``. Execute the following command to upload the whole
family to the database::
  
  verdi data vasp-potcar uploadfamily --path=$HOME/myaiida/potpaw_PBE.54.tar --name=PBE.54 --description="PBE potentials for version 5.4.4"
  skipping file /home/username/potpaw_PBE.54/H_AE/POTCAR - uploading raised <type 'exceptions.IndexError'>list index out of range
  POTCAR files found: 327. New files uploaded: 326, Added to Family: 326

You can ignore the error, which is caused by one of the potentials not complying to the usual standard.
  
The ``name`` and ``description`` are not optional and have to be
specified. The ``path`` could be either an archive, or one could use
a folder name.  It is also possible, not to specify path, but then you the plugin will
traverse all folders from the one the above command is executed from and try to upload all
potentials it finds to the given family.
