
The interface
=============


The basic function of the interface to the Wien2k program package is
to take the output of the program that constructs the projected local
orbitals (:program:`dmftproj`), and to store all the necessary information into
an hdf5 file. This latter file is then used to do the DMFT calculation. The
reason for this structure is that this enables the user to have everything
that is necessary to reproduce the calculation in one single hdf5 arxive.

.. index:: Interface to Wien2k

.. _interfacetowien:

The interface to Wien2k
-----------------------

As explained above, this interface produces an hdf5 arxive out of the files that
were written by the band structure package :program:`Wien2k/dmftproj`. 
For this purpose we
use the python module :class:`SumK_LDA_Wien2k_input`. It is initialised as::

  from pytriqs.dft.converters.wien2k_converter import *
  Converter = SumK_LDA_Wien2k_input(Filename = material_of_interest)

The only necessary parameter to this construction is the parameter `Filename`.
It has to be the root of the files produces by dmftproj. For example, if you did a 
calculation for TiO, the :program:`Wien2k` naming convention is that all files are called 
:file:`TiO.*`, so you would give `Filename = TiO`. The constructor opens
an hdf5 arxive, named :file:`material_of_interest.h5`, where all the data is stored.

There are three optional parameters to the Constructor:

  * `LDASubGrp`: We store all data in sub groups of the hdf5 arxive. For the main data
    that is needed for the DMFT loop, we use the sub group specified by this optional parameter.
    If it is not given, the standard value `SumK_LDA` is used as sub group name.
  * `SymmSubGrp`: In this sub group we store all the data for applying the symmetry 
    operations in the DMFT loop. Standard value is `SymmCorr`.
  * `repacking`: If true, and the hdf5 file already exists, the system command :program:`h5repack` 
    is invoked. This command ensures a minimal file size of the hdf5
    file. Standard value is `False`. If you want to use this, be sure
    that :program:`h5repack` is in your path variable!

After initialising the interface module, we can now convert the input text files into the
hdf5 arxive by::

  Converter.convert_DMFT_input()

This reads all the data, and stores it in the sub group `LDASubGrp`, as discussed above. 
In this step, the files :file:`material_of_interest.ctqmcout` and :file:`material_of_interest.symqmc`
have to be present in the working directory.

After this step, all the necessary information for the DMFT loop is stored in the hdf5 arxive, where
the string variable `Converter.HDFfile` gives the file name of the arxive.
You can now proceed with :ref:`LDADMFTmain`.


Data for post-processing
------------------------

In order to calculate some properties using the DMFT self energy, several other routines are
used in order to convert the necessary input from :program:`Wien2k/dmftproj`. For instance, for 
calculating the partial density of states or partial charges consistent with the definition
of :program:`Wien2k`, you have to use::

  Converter.convert_Parproj_input()

This reads the files :file:`material_of_interest.parproj` and :file:`material_of_interest.sympar`.
Again, there are two optional parameters

  * `ParProjSubGrp`: The sub group, where the data for the partial projectors is stored. Standard
    is `SumK_LDA_ParProj`.
  * `SymmParSubGrp`: Sub group for the symmetry operations, standard value is `SymmPar`.

Another routine of the class allows to read the input for plotting the momentum-resolved
spectral function. It is done by::
  
  Converter.convert_bands_input()

The optional parameter, which tells the routine where to store the data is here `BandsSubGrp`, 
and its standard value is `SumK_LDA_Bands`.

After having converted this input, you can further proceed with the :ref:`analysis`.

MPI issues
----------

The interface package is written such that all the operations are done only on the master node.
The broadcasting to the nodes has to be done by hand. The :class:`SumK_LDA`, described in the
following section, takes care of this automatically.

Interfaces to other packages
----------------------------

Because of the modular structure, it is straight forward to extend the TRIQS package 
in order to work with other band-structure codes. The only necessary requirement is that 
the interface module produces an hdf5 arxive, that stores all the data in the specified
form. For the details of what data is stored in detail, see the reference manual.
