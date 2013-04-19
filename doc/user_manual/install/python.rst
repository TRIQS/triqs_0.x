.. _python_install:

Python and the scientific python tools
====================================================

On some machines, not all the necessary scientific python tools are present
(e.g. old distribution, computation centers, etc....), or not in the right version.

There is a way to efficiently solve this issue (other than recompiling your own python, and 
all the libs, which is painful) and specially **avoiding any upgrade on your system** (which 
is crucial on cluster machines e.g.).

.. _virtualenv:
.. highlight:: bash

Using virtualenv to install/upgrade some python packages locally
----------------------------------------------------------------------------

In the frequent cases where only one or two packages are missing, there is a simple solution
called virtualenv.

Principle
^^^^^^^^^^^^

virtualenv is a python package that allow you to install in your home a simple copy (or mirror)
or python, and install and upgrade some packages for you, without touching the systems.

.. warning:: This works as long as the other dependencies of the packages on e.g. C++ libraries are satisfied.

Usage 
^^^^^^^^

* First install virtualenv, or ask the system manager to do it.
  On ubuntu ::

    sudo apt-get install python-virtualenv
 
* To create a new python in your home, do e.g. ::

    virtualenv ~/.my_python --system-site-packages
 
  This create in ~/.my_python a clone of the python of the system, with its packages.

  You can then use it with ::
   
    ~/.my_python/bin/python

* To make this your default python, just add in your ~/.bash_aliases (or bashrc) ::

    source $HOME/.my_python/bin/activate

  Then (reopen your terminal ...) ::

    which python 
  
  should return typically::
    
    /home/ME/.my_python/bin/python

* You can install any python package for this local python installation, e.g. ::

   easy_install --upgrade ipython

  and you have the latest ipython notebook....


