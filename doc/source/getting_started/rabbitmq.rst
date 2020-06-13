.. _messaging:

===================
3. Install RabbitMQ
===================

`AiiDA`_ relies on `RabbitMQ`_. The installation of this requires super user privileges. On Fedora 30 or later you typically install this with::

  $ sudo dnf install rabbitmq-server
  $ sudo systemctl enable rabbitmq-server
  $ sudo systemctl start rabbitmq-server

On Debian (e.g. Ubuntu) systems you would only have to execute::

  $ sudo apt-get install rabbitmq-server

We can check if `RabbitMQ`_ is running by issuing::

  $ sudo rabbitmqctl status

If it is not running on your Debian system, try a restart. Now, `RabbitMQ`_ is enabled and
running. It should start every time you restart your computer. In case you run into problems,
consult the detailed installation instructions for `RPM based systems`_ (CentOS, Fedora, etc.)
and `Debian based systems`_.

.. _RabbitMQ: https://www.rabbitmq.com/
.. _AiiDA: https://www.aiida.net
.. _RPM based systems: https://www.rabbitmq.com/install-rpm.html
.. _Debian based systems: https://www.rabbitmq.com/install-debian.html
