
TRIQS as a library
==================

Let's say you want to write a new code that uses the libraries of TRIQS
but that is not meant to be a part of TRIQS. Here's how you should proceed
in order to do that.

First of all you will have to install the development version of TRIQS.
Let's be concrete and say that you download the TRIQS sources in
TRIQS/triqs_dev::

   mkdir TRIQS
   cd TRIQS
   git clone git@bitbucket.org:mferrero/triqs_dev.git

Now we will have to compile and install TRIQS, but because we want
to use TRIQS as a library, we have to make sure the
*Install_dev* option is turned on::

  cd TRIQS/triqs_dev
  ccmake .
  (turn the *Install_dev* option to ON)
  (turn the *Install_dev* option to ON)
  make -j 4
  make test
  make install

  mkdir build
  cd build

turn on the *Install_dev* flag of the
