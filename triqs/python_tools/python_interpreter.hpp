
/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011 by M. Ferrero, O. Parcollet
 *
 * TRIQS is free software: you can redistribute it and/or modify it under the
 * terms of the GNU General Public License as published by the Free Software
 * Foundation, either version 3 of the License, or (at your option) any later
 * version.
 *
 * TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
 * WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
 * FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
 * details.
 *
 * You should have received a copy of the GNU General Public License along with
 * TRIQS. If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/

#ifndef PYTHON_INTERPRETER_H_2h378f
#define PYTHON_INTERPRETER_H_2h378f

#include "Python.h"
#include "numpy/arrayobject.h"

// from PEP353. <2.5 compatibility ....
#if PY_VERSION_HEX < 0x02050000 && !defined(PY_SSIZE_T_MIN)
typedef int Py_ssize_t;
#define PY_SSIZE_T_MAX INT_MAX
#define PY_SSIZE_T_MIN INT_MIN
#endif

#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <complex>
#include <boost/python.hpp>


namespace triqs { namespace python_tools { 

 /**
   This class helps to write main.cpp to embed the python interpreter.
   Only one object of this type can be created in a program.
   */
 class python_interpreter { 
  int argc_; char **argv_; static bool active;
  boost::python::dict * AuxiliaryPythonFunctions;

  public:
  /**
   * argc, argv are the usual main arguments.
   * Additionnal_Python_Path : any white-space separated list of paths to add the sys.path at initialisation
   */
  python_interpreter(int argc, char **argv,  std::string Additionnal_Python_Path = "");

  ///
  ~python_interpreter();

  /** Change the prompt of the interpreter */
  python_interpreter & set_prompt (std::string const & prompt);

  /// Tag the svn version, etc.. 
  void print_greeting(); 

  /**
   * Execute some code into the interpreter.
   * e.g. *this << "python code"; will execute "python code" as if you type it in the interactive interpreter.
   */
  python_interpreter & operator<< (std::string const & do_it) { PyRun_SimpleString(do_it.c_str()); return *this; }

  /**
   * Execute some code into the interpreter from the file named filename. 
   * Returns true iif file can be executed without error
   */
  bool execute_file (const char * filename);

  /**
   * Execute some code into the interpreter from the file named filename. 
   * Returns true iif file can be executed without error
   */
  bool execute_file (std::string const & filename) { return execute_file(filename.c_str());}

  /** 
   * Gives the control of execution to the Python interpreter through an ipython shell.
   *  ./my_python_interpreter script.py options  : runs as usual, just pure python
   *  ./my_python_interpreter -i script.py options : runs the scripts, and then launch ipython shell at the end of the script.
   */
  void mainloop();

  /** 
   * Gives the control of execution to the Python interpreter.  
   *   It also includes the support for ipython.
   *   ./my_python_interpreter script.py options  : runs as usual, just pure python
   *   ./my_python_interpreter -i script.py options : runs the scripts, and gets normal interactive python at the end of the script.
   */
  void mainloop_ipython();

  /**
   *  Defines a new function in the interpreter in a private dict, not accessible from python level
   *  name : name of the function
   *  args : arguments
   *  code : the code
   *  It is equivalent to say in python : 
   *   def name (args):
   *     code
   */
  void def_auxiliary_function(std::string name, std::string args, std::string code);

  /**
   * Get back the auxiliary function
   */
  boost::python::object get_auxiliary_function(std::string F);

  };

} }

#endif

