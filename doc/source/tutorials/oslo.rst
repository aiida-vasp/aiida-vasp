.. _oslo:

======================
Workshop Oslo initials
======================

We will utilize the AWS resources that has already been distributed for the `AiiDA part of the workshop`_.

Please follow these steps to get set up AiiDA-VASP.

1. You need to _`reset your password` that is set randomly at the computing cluster. This is different from the Jupyter instance you have been using so far. Your ``username`` in this session should be the one you entered as a username. Then you will receive an SMS with a new password that you should use in the next steps. Hopefully you entered a valid telephone number during submission of the application.

2. When the password is reset you need to take your private key `aiida_tutorial_NN` and `aiida_tutorial_NN.pub` and upload the files to the `/home/max.ssh/` folder at the AWS resources and make sure you rename them to `aiida_tutorial` and `aiida_tutorial.pub`

3. Then you need to log in to the Saga computing cluster in Norway using your ``username`` and the password you received on SMS.::

     ssh username@saga.sigma2.no

4. Answer `yes` if it ask you a question. And put the content of `aiida_tutorial_NN.pub` into the ``authorized_keys`` file that you find in the `~/.ssh` folder on the Saga cluster.

5. Log out from the cluster and log in again and make sure you do not need to enter the password.

6. Log into the AWS service and try the same with the command::

     ssh -i ~/.ssh/aiida_tutorial username@saga.sigma2.no

7. Install AiiDA-VASP with the command::

     pip install aiida-vasp

8. Now we need to upload the potentials. For this tutorial we only need the VASP tutorial POTCAR files::

     wget https://github.com/aiida-vasp/aiida-vasp/blob/develop/tutorials/vasp_potentials.tar.gz

9. Untar it::

     tar xvzf vasp_potentials.tar.gz

10. Go into the ``vasp_potentials`` directory. We will now upload the potentials using the custom commands created for the POTCAR data types.::

      verdi data vasp-potcar uploadfamily --name pbe --description "A few tutorial PBE potentials"

11. This will upload the potentials into the database and hash them. E.g. you will not be able to have multiple entries of one potential in the database. Also, in the calculations, only the hash is used such that the POTCAR data (which is covered by license) is not revealed. This should complete with a message that 3 potentials was found and uploaded. You are now ready to try a few calculations in AiiDA-VASP.
   
.. _reset your password: https://www.metacenter.no/user/reset/
.. _AiiDA part of the workshop: https://aiida-tutorials.readthedocs.io/en/latest/pages/2019_SINTEF/index.html 
