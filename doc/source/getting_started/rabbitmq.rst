.. _rabbitmq:

================
Install RabbitMQ
================

`AiiDA`_ relies on `RabbitMQ`_. The installation of this requires super user privileges. On Fedora 30 or later you typically install this with::

  sudo dnf install rabbitmq-server
  sudo systemctl enable rabbitmq-server
  sudo systemctl start rabbitmq-server

Now, `RabbitMQ`_ is enabled and running. It should start every time you restart your computer.

.. _RabbitMQ: https://www.rabbitmq.com/
.. _AiiDA: https://www.aiida.net
