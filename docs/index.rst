Documentation for backscatter
=============================

Overview
--------
backscatter contains various packages to process and analyze
SuperDARN data. This package aims to make it easier for users
to use, modify, or understand the software needed to process data
compared to the aging RST C code package.

While people are free to contribute new analysis code, the main
goals of backscatter are to offer the following:

- A way to read and write all SuperDARN DMAP formatted files
- Processing of .rawacf files into .fitacf files
- Processing of .fitacf files into .mapfiles needed to produce convection map plots


Installation
------------
backscatter can be installed or used in a couple of different ways.

The package can be cloned from https://github.com/SuperDARNCanada/backscatter.
In order to work with the package one can simply enter the directory and work
from there, or add it to the Python path.

backscatter can also be installed as a package by running
``python setup.py install`` from within the directory OR
via pip VCS install by running ``pip install git+git://github.com/SuperDARNCanada/backscatter.git`` without having to clone anything!

backscatter makes use of configuration files at import time, and installation
creates a system wide configuration file in /etc/backscatter as well as the
user's home directory.

The setup scripts have dependency handling in them, but in case that fails
or if you choose not to install, the following dependencies are needed:

- setuptools
- ConfigParser
- Numpy >= v1.8

*Note that depending on your system, installation may require root privileges.*

Usage
-----

When importing backscatter, the package will attempt to locate a configuration
file in the current directory firstly, the user's home directory secondly, and /etc/backscatter thirdly. This allows a system wide configuration for 
standard operation while allowing each user to override any options they choose.

backscatter also requires a folder containing the hdw.dat files for each radar
somewhere on the system. The location can be set in the configuration file with the default being /usr/local/hdw.dat. Please clone https://github.com/vtsuperdarn/hdw.dat to /usr/local, or to another directory and reflect the change in your configuration files. Backscatter will fail to import if it can't find the hardware files.

Once successfully imported, you can start using package contents in your scripts.

If installed, backscatter can also be used in a command line way that can be used for batch processing, or as a replacement for RST commands in current scripts. For example, to use fitacf as a utility from command line we can run ``python -m backscatter.fitacf.fitacf in_file out_file`` from anywhere.



Contents:
---------

.. toctree::
   :maxdepth: 2

   backscatter


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

