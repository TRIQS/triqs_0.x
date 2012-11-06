.. _EPD:
.. highlight:: bash

Enthought Python Distribution (EPD)
---------------------------------------------

When the python itself is tool old, or not present, we recommend
the `Enthought Python Distribution <http://www.enthought.com/>`_ 
which provides a complete package of python and scientific tools.
It is available for most platforms.

Enthought's distribution allows you to have a recent python, ipython with all the necessary libraries
(including the ipython notebook and its dependencies) without upgrading your machine's system (in particular python).

* Take the `EPD Free version <http://www.enthought.com/products/epd_free.php>`_  (the free and light version of Enthought).

* Install it ::

   bash epd_free.sh -b -p EPD_Free_Install_Dir

* Add the following packages (the last two are only to compile the documentation) ::

   EPD_Free_Install_Dir/bin/easy_install h5py
   EPD_Free_Install_Dir/bin/easy_install Sphinx
   EPD_Free_Install_Dir/bin/easy_install pyparsing
   EPD_Free_Install_Dir/bin/easy_install cython

* When compiling TRIQS, you will simply pass the option ::

    cmake ..... -DPYTHON_INTERPRETER=EPD_Free_Install_Dir/bin/python

 and that is all.

.. warning ::
 
 A priori, you could also use the regular, full Enthought distribution (it is free for academics), but we do not recommend it
 at present. Indeed, in this distribution the HDF5 has been compiled without C++ support, which TRIQS requires.
 So in most machines it will work, with h5py compiled with the Enthought HDF5 version, while the C++ code will silently include
 your system's HDF5 headers, if they exist ...



