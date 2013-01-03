
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

#include "python_interpreter.hpp"

// blitz array
#include <triqs/deprecated/blitzext/pyarray_converters.hpp>

#define AS_STRING(X) AS_STRING2(X)
#define AS_STRING2(X) #X

using namespace  boost::python;

//------------------------------------------------------------ 

namespace triqs { namespace python_tools { 

 void translatorString (std::string const& x) { PyErr_SetString(PyExc_RuntimeError, x.c_str());};

 bool python_interpreter::active = false;

 //------------------------------------------------------------ 

 python_interpreter::python_interpreter(int argc, char **argv, std::string Additionnal_Python_Path):
  argc_(argc), argv_(argv) 
 {
  if (active) throw std::string("Internal error : only one python_interpreter can be created at a time");
  active= true;

  Py_Initialize(); // init the interpreter

  // init the boost python array converters
  triqs::blitzext::init_boost_array_converters();
  // triqs_arrays::BoostPythonConverters::Register_Array_Boost_Converters();

  // register the exception translation
  // register_exception_translator<OP_Tools::GF_Exception>(translatorGF);
  register_exception_translator<std::string>(translatorString);

  // I set up the sys.argv variable first here because it is used in the init of the tables modules used e.g. in HDFArchive
  PySys_SetArgv(argc,argv);

  // init the private dict 
  AuxiliaryPythonFunctions = new boost::python::dict(); 
  // add one function to call by keyword
  def_auxiliary_function("KeywordCall","O_to_be_called, kwargs"," return O_to_be_called (**kwargs)");

  //repair the sys.executable
  //std::stringstream fs; fs<<"import glob,sys,os; sys.executable = os.path.abspath(glob.glob('"<<argv[0]<<"')[0])";
  std::stringstream fs; fs<<"import sys,os; sys.executable = os.path.abspath(os.environ['_'])";
  *this<< fs.str();

  // Add Additionnal_Python_Path to the sys.path
  //std::string PYTHON_PATH_ADD_= "import sys; sys.path =[\""+ std::string( AS_STRING(PYTHONPATH_ADD)) +"\"] + '" + Additionnal_Python_Path +"'.split() +  sys.path"; 
  std::string PYTHON_PATH_ADD_= "import sys; sys.path = '" + Additionnal_Python_Path +"'.split() +  sys.path"; 
  *this << PYTHON_PATH_ADD_;

#ifndef __HOSTNAME__
#define __HOSTNAME__ "UNKOWN"
#endif

#ifndef __COMPILEDBY__
#define __COMPILEDBY__ "UNKOWN"
#endif
 }

 //------------------------------------------------------------ 

 python_interpreter::~python_interpreter() { 
  delete AuxiliaryPythonFunctions;
 }

 //------------------------------------------------------------ 

 void python_interpreter::print_greeting() { 
  std::stringstream fs; 
  fs<< "myprint = lambda x : sys.stderr.write(\"%s\\n\"%x)\n" 
   << "myprint ('Compiled on "<<__DATE__<<" at " 
   << __TIME__<<" on machine "<< AS_STRING(__HOSTNAME__)<< " by "<< AS_STRING(__COMPILEDBY__) " ')\n";
  *this<<fs.str();
 }

 //------------------------------------------------------------ 

 bool python_interpreter::execute_file( const char * filename ) { 
  FILE *fp = fopen(filename,"r");
  if (!fp)  return false;
  PyErr_Clear();
  if (!PyRun_SimpleFile(fp,(char *)filename)) {PyErr_Print(); return false;}
  return true;
 }

 /*
    void python_interpreter::executeScript(bool requireScript) { 
    if (argc_<=1) return;
    FILE *fp = fopen(argv_[1],"r");
    if (!fp)  { 
    if (requireScript)
    {TRIQS_RUNTIME_ERROR<<"File "<< argv_[1]<<" not found ";}
    else 
    {cout<<"File "<< argv_[1]<<" not found "<<endl;}
    }
    else {
    PyErr_Clear();
    if (!PyRun_SimpleFile(fp,(char *)argv_[1])) {PyErr_Print();}
    }
    }
    */

 //------------------------------------------------------------ 

 void python_interpreter::mainloop_ipython() { 

   PyRun_SimpleString( 
    "import sys\n"
    "import IPython\n"
    "if IPython.__version__ >= '0.11' : \n"
    "   from IPython.frontend.terminal.ipapp import launch_new_instance\n"
    "   sys.exit(launch_new_instance())\n"
    "else:\n"
    "  import IPython.Shell\n"
    "  IPython.Shell.start().mainloop()\n"
    );
     //"import IPython\n"
     //"IPython.embed()\n");
     //"IPython.Shell.start(locals()).mainloop()\n");

 }

 //------------------------------------------------------------ 

 void python_interpreter::mainloop() { 
  // now start the interpreter main loop and start executing the commands given to the new interpreter...
  Py_Main(argc_,argv_);
 }

 //------------------------------------------------------------ 

 python_interpreter & python_interpreter::set_prompt (const std::string & prompt) {
  boost::python::object sys_module =  boost::python::import("sys");
  boost::python::setattr(sys_module,"ps1", boost::python::object(prompt.c_str())); 
  return *this;
 }

 //------------------------------------------------------------ 

 boost::python::object python_interpreter::get_auxiliary_function(std::string F){ assert(active); return (*AuxiliaryPythonFunctions)[F]; }

 //------------------------------------------------------------ 

 void python_interpreter::def_auxiliary_function(std::string name, std::string args, std::string code) { 
  std::string s = "def " + name + "(" + args + " ) : \n " + code + "\n";
  boost::python::object result = exec(s.c_str(), *AuxiliaryPythonFunctions,*AuxiliaryPythonFunctions);
 }
} }


