.. _database:

================================
Install and configure a database
================================

`AiiDA`_ relies on a `Postgresql`_ database. You can run this locally on your computer, or use an external service. This is up to you.
If you are sharing or performing co-calculation with colleagues, we recommend that you perform day-day calculations on a local database
and export the results you want to share when it is convenient. This can then easily be imported in a central database for everyone to access.

What you need to configure the database in `AiiDA`_ is the hostname or IP address of the `Postgresql`_ database server,
a username, a password, a database name and a port number. Please notice that the password can be intercepted, so do not chose something
that you typically use for high security services, even though it is only used on your local computer.

`AiiDA`_ uses the `psycopg2`_ adapter to access the `Postgresql`_ database server from Python.

On Fedora 30, please follow these steps to configure a usable database for `AiiDA`_ for the first time.

#. Let us install `Postgresql`_::

     sudo dnf install postgresql
     sudo dnf install postgresql-server

#. Then we initialize the service::

     sudo initdb

#. And then enable (so that it restarts if you computer restarts) and start it::

     sudo systemctl enable postgresql
     sudo systemctl start postgresql
     
#. The installation created the user ``postgres`` which you can use to administer the databases, users etc.
   Let us now add the database ``aiidadb`` and the user ``aiida`` and give them the necessary privileges and passwords.
   First, let us get into the `Postgresql`_ interactive terminal::

     sudo su - postgres
     psql

#. Now we need to do something important. We need to change the ``ident`` entries in the ``pg_hba.conf`` file to ``md5``.
   In order to locate this file, let us ask `Postgresql`_::
     
     SHOW hba_file;

#. Say the location is ``/var/lib/pgsql/data/pg_hba.conf``, which is pretty typical on Fedora 30. Exit the
   `Postgresql`_ interactive terminal with ``\q``. Also, make sure you exit from the ``postgres`` user. When you are back
   as the regular user, replace all entries
   of ``ident`` with ``md5`` by running the command::

     sudo sed -i 's/ident/md5/g' /var/lib/pgsql/data/pg_hba.conf

   The ``md5`` obscures your password so that it is not sent over the network in clear text. However, it can still be
   obtained from your `AiiDA`_ config file. Also, we would like to note that ``md5`` is outdated and not anymore secure.
   How secure you need this to be is of course something you should consider. The more modern way is ``scram-sha-256``, however
   some adapters does not yet support this since it is a new option that requires `Postgresql`_ 10 or later. This goes for `AiiDA`_
   as well. Currently only ``md5`` is supported. However, according to :ref:`this post <https://www.postgresql.org/docs/11/auth-password.html>`
   it seems as long as ``md5`` is present in ``pg_hba.conf`` the client side determines the encryption used.

#. Now, let us go back into the Postgreqsl interactive terminal again to create the database and user::

     sudo su - postgres
     psql
     CREATE USER aiida WITH PASSWORD '<password>';
     CREATE DATABASE aiidadb OWNER aiida ENCODING 'UTF8' LC_COLLATE='en_US.UTF-8' LC_CTYPE='en_US.UTF-8' TEMPLATE=template0;
     GRANT ALL PRIVILEGES ON DATABASE aiidadb to aiida;

   You supply the password of your choice as ``<password>``. Please notice you need the semicolon after the commands. We have now created
   the necessary database ``aiidadb`` and the user ``aiida`` and given all permissions to the user ``aiida`` on the database ``aiidadb``.
   Also note that we have forced a specific encoding. Please do not forget to do this.

We can now continue to create the profile.

.. _Postgresql: https://www.postgresql.org/
.. _psycopg2: https://github.com/psycopg/psycopg2
.. _AiiDA: https://www.aiida.net
