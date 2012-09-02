
TRIQS as a library
==================

.. highlight:: c

Let's say you want to write a new code that uses the libraries of TRIQS
but that is not meant to be a part of TRIQS. Here's how you should proceed
in order to do that.


Setting the stage
-----------------

Just to be clear, let me start by showing the directory structure
I will use:

.. code-block :: bash

   /home/triqs/triqs_dev --> sources of TRIQS
   /home/triqs/build --> the directory where TRIQS will be compiled
   /home/triqs/install --> TRIQS is installed here
   /home/project/src --> the source code of my project
   /home/project/build --> the directory where I will compile my code

Basically, anything that has to do with the sources of TRIQS, its compilation
and its installation will be under :file:`/home/triqs`, while the new project
will be under :file:`/home/project`. Obviously you can choose any other directory
but this structure will be assumed below.

Install development version of TRIQS
------------------------------------

First of all you will have to install the *development* version of TRIQS.
Let's be concrete and say that you want to download the TRIQS sources in
:file:`/home/triqs/triqs_dev`:

.. code-block :: bash

   $ cd /home
   $ mkdir triqs
   $ cd triqs
   $ git clone git@bitbucket.org:mferrero/triqs_dev.git

Now we want to compile and install TRIQS from these sources,
but because we want to use TRIQS as a library, we have to
make sure the ``Install_dev`` option is turned on. Moreover
we decide to install everything in :file:`/home/triqs/install`:

.. code-block :: bash

  $ cd /home/triqs
  $ mkdir build
  $ cd build
  $ cmake -DBOOST_SOURCE_DIR=path_to_boost -DInstall_dev=ON \
          -DCMAKE_INSTALL_PREFIX=/home/triqs/install ../triqs_dev
  $ make -j 4
  $ make test
  $ make install

If everything went fine, the TRIQS libraries and headers are installed
in :file:`/home/triqs/install`.

Start your project
------------------

Now we want to start our project. We will do this in a directory
:file:`/home/project`. We put the sources in :file:`/home/project/src` and
later build (compile) the project in :file:`/home/project/build`.
Let's start by adding some sources:

.. code-block :: bash

  $ cd /home
  $ mkdir project
  $ cd project
  $ mkdir src
  $ cd src

OK, our project will be just one :file:`main.cpp` file::

  #include <triqs/arrays/array.hpp>
  #include <triqs/arrays/proto/array_algebra.hpp>

  namespace tqa = triqs::arrays;

  int main(){

    tqa::array<double,1> A(10), B(10);
    A()=2; B()=3;
    tqa::array<double,1> C = A+B;
    std::cout << "C = "<< C << std::endl;

  }

As you can see, the code includes headers from TRIQS. Along with
:file:`main.cpp` we write a :file:`CMakeLists.txt` file to compile our project.
In order to make this easy, there is a file called :file:`TRIQSConfig.cmake`
file in :file:`/home/triqs/install/share/triqs/cmake`. Including this file in
your :file:`CMakeLists.txt` automatically defines a certain number of useful
variables, especially the include directories related to the TRIQS headers and
the location of the TRIQS libraries. Here is what your simple
:file:`CMakeLists.txt` can be:

.. code-block :: bash

  cmake_minimum_required(VERSION 2.8)

  project(myproj CXX)

  set(CMAKE_BUILD_TYPE Release)

  SET(CMAKE_INSTALL_RPATH "${TRIQS_PATH}/lib")
  SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
  SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)

  include(${TRIQS_PATH}/share/triqs/cmake/TRIQSConfig.cmake)

  add_executable(example main.cpp)

  include_directories(${TRIQS_INCLUDE} ${EXTRA_INCLUDE} ${CBLAS_INCLUDE} ${FFTW_INCLUDE})
  target_link_libraries(example ${TRIQS_LIBRARY} ${EXTRA_LIBRARY})

We're all set! Everything is ready to compile our project. If we want to build
everything in :file:`/home/project/build`, we do as follows:

.. code-block :: bash

  $ cd /home/project
  $ mkdir build
  $ cd build
  $ cmake -DTRIQS_PATH=/home/triqs/install /home/project/src
  $ make
  $ ./example

That's it! You can modify your sources and then recompile with make. Obviously
with bigger projects you :file:`CMakeLists.txt` file will change, but the
principle remains the same.

