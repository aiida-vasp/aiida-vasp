.. _profile:

===================
7. Create a profile
===================

Please consult the documentation of AiiDA on how to create an `AiiDA profile`_. If you followed our instructions when creating the
database, please use the following information when requested by ``verdi setup``::

  Database engine: postgresql_psycopg2
  PostgreSQL host: localhost
  PostgreSQL port: 5432
  AiiDA Database name: aiidadb
  AiiDA Database user: aiida
  AiiDA Database password: <password>

In the rest of this getting started guide we will assume the user registered with the email
address ``mymail@address.com``, which is later used as an identifier (check e.g. next step when
we test the computer connection).

After setting up your profile, now AiiDA daemon can be started by

::

   $ verdi daemon start
   Starting the daemon... RUNNING

.. _AiiDA profile: https://aiida-core.readthedocs.io/en/latest/install/installation.html#setup-instructions
