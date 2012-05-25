
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

#include <triqs/gf_local/MeshGF.hpp>
#include <triqs/gf_local/TailGF.hpp>
#include <triqs/gf_local/GF_Bloc_Base.hpp>
#include <triqs/gf_local/GF_Bloc_ImFreq.hpp>
#include <triqs/gf_local/GF_Bloc_ReFreq.hpp>
#include <triqs/gf_local/GF_Bloc_ImTime.hpp>
#include <triqs/gf_local/GF_Bloc_ReTime.hpp>
#include <triqs/gf_local/GF_Bloc_ImLegendre.hpp>
#include <triqs/gf_local/GF_C.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/def.hpp>
#include <boost/shared_ptr.hpp>

using namespace python;

#define  decl_Base(DataType,NAME)\
  class_<GF_Bloc_Base<DataType> >(NAME, init<object,object,PyObject *,boost::shared_ptr<MeshGF>, boost::shared_ptr<TailGF> >()) \
    .def_readonly("_data_c_array", &GF_Bloc_Base<DataType>::data_as_numpy) \
    .def_readonly("_tail",  &GF_Bloc_Base<DataType>::tail_ptr)		\
    .def_readonly("mesh",  &GF_Bloc_Base<DataType>::mesh_ptr)		\
    .def_readonly("N1",  &GF_Bloc_Base<DataType>::N1)			\
    .def_readonly("N2",  &GF_Bloc_Base<DataType>::N2)			\
    .def_readonly("_IndicesL",  &GF_Bloc_Base<DataType>::IndicesL)	\
    .def_readonly("_IndicesR",  &GF_Bloc_Base<DataType>::IndicesR)	\
    .def_readonly("Statistic",  &GF_Bloc_Base<DataType>::Statistic)	\
    .def_readonly("Beta",  &GF_Bloc_Base<DataType>::Beta)		\
    .def("copyFrom", &GF_Bloc_Base<DataType>::operator=)			\
    .def("save",&GF_Bloc_Base<DataType>::save,save_overlo("Save the Green's function into text files.")) \
    .def("load",&GF_Bloc_Base<DataType>::load,"Load the Green function from text files on all nodes. Inverse of save") \
    .def("zero",&GF_Bloc_Base<DataType>::zero,"Puts the GF to 0")	\
 ; 

BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(save_overlo, save, 1,2);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(F_overloads, setFromFourierOf, 1, 2);
BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(Finv_overloads, setFromInverseFourierOf, 1, 2);

void translatorString(std::string const& x) { PyErr_SetString(PyExc_RuntimeError, x.c_str()); };

#include <triqs/deprecated/blitzext/pyarray_converters.hpp>

BOOST_PYTHON_MODULE(_pytriqs_GF) {

  register_exception_translator<std::string>(translatorString);
  triqs::blitzext::init_boost_array_converters();

  docstring_options doc_options;
  doc_options.disable_cpp_signatures();
  doc_options.disable_py_signatures();

  enum_<Type_GF>("GF_Type")
    .value("Real_Time", Real_Time)
    .value("Imaginary_Time", Imaginary_Time)
    .value("Imaginary_Frequency",Imaginary_Frequency)
    .value("Real_Frequency",Real_Frequency)
    .value("Imaginary_Legendre",Imaginary_Legendre)
    ;

  enum_<Statistic_GF>("GF_Statistic")
    .value("Boson", Boson)
    .value("Fermion",Fermion)
    ;

  // **********  MeshGF  ******************

  class_<MeshGF, boost::shared_ptr<MeshGF> >("MeshGF", init<Type_GF,Statistic_GF,double,PyObject *>()) 
    .def_readwrite("Statistic", &MeshGF::Statistic)
    .def_readonly("TypeGF",  &MeshGF::typeGF)
    .def_readonly("Beta",  &MeshGF::Beta)
    .def("__len__",&MeshGF::len)
    .def("__iter__",&MeshGF::__iter__)
    .def("next",&MeshGF::next)
    .def("Bose_Fermi",&MeshGF::Bose_Fermi)
    .def("__reduce__",&MeshGF::__reduce__)
    .def("__reduce_to_dict__",&MeshGF::__reduce_to_dict__)
    .def("__factory_from_dict__",&MeshGF::__factory_from_dict__)
    .staticmethod("__factory_from_dict__")
    .def(self == other<MeshGF>())
    ;

  // **********  TailGF  ******************

  class_<TailGF, boost::shared_ptr<TailGF> >("TailGF", init<int,int,boost::python::list,boost::python::list>() )
    .def (init<const TailGF &, boost::python::object, boost::python::object>())
    .def ("copyFrom",&TailGF::copyFrom)
    .def ("invert",&TailGF::invert,"Replace itself by the expansion of the inverse of the function.")
    .def ("__call__",&TailGF::__call__,"Sets the expansion to 0")
    .def ("zero",&TailGF::zero,"Sets the expansion to 0")
    .def("changeOrderMax",&TailGF::changeOrderMax, "Changes the orderMax")
    .def("Coefs",&TailGF::AllCoefs)
    .def("__repr__",&TailGF::__repr__)
    .def("__getitem__",&TailGF::__getitem__)
    .def("__setitem__",&TailGF::__setitem__)
    .add_property("OrderMin",&TailGF::OrderMin)
    .add_property("OrderMax",&TailGF::OrderMax_py)
    .def("__reduce__",&TailGF::__reduce__)
    .def("__reduce_to_dict__",&TailGF::__reduce_to_dict__)
    .def("__factory_from_dict__",&TailGF::__factory_from_dict__)
    .staticmethod("__factory_from_dict__")
    .def(self == other<TailGF>())
    .def("from_L_T_R",&TailGF::from_L_T_R,"Multiply Left & right with Matrices :  G <- L* G2 * R")
    .def("transpose",&TailGF::transpose,"Transpose the array : new view as in numpy")
    .def("conjugate",&TailGF::conjugate ,"Transpose the array : new view as in numpy")
    .def(self += other<TailGF>())
    .def(self -= other<TailGF>())
    .def(self *= other<TailGF>())
    ;
 
  // **********   Base ******************

  decl_Base(COMPLEX,"_GF_Bloc_Base_C_C");
  decl_Base(double,"_GF_Bloc_Base_R_C");

  // **********   ImFreq ******************

  char imfreq_doc_density[] = "Computes the density :math:`G(\\tau = 0^-)`";
  char imfreq_doc_fourier[] = "\
Sets to the Fourier transform of Gt.\n\
\n\
:**Parameters**: - Gt - an imaginary-time Green's function\n\
:**Return type**: a new Green's function\n";
  char imfreq_doc_legendre[] = "Transforms from Legendre and sets it's tail";

  class_<GF_Bloc_ImFreq, bases<GF_Bloc_Base<COMPLEX> > >("GFBloc_ImFreq", init<object,object,PyObject *,boost::shared_ptr<MeshGF>, boost::shared_ptr<TailGF> >()) 
    .def("density",&GF_Bloc_ImFreq::density, imfreq_doc_density)
    .def("setFromFourierOf",&GF_Bloc_ImFreq::setFromFourierOf,F_overloads(imfreq_doc_fourier))
    .def("setFromLegendre",&GF_Bloc_ImFreq::setFromLegendre,imfreq_doc_legendre)
    ;

  // **********   ReFreq ******************

  char refreq_doc_pade[] = "Sets to the analytic continuation of Gw using Pade approximants.";

  class_<GF_Bloc_ReFreq, bases<GF_Bloc_Base<COMPLEX> > >("GFBloc_ReFreq", init<object,object,PyObject *,boost::shared_ptr<MeshGF>, boost::shared_ptr<TailGF> >()) 
    .def("density",&GF_Bloc_ReFreq::density, "Computes the density :math:`G(\\tau = 0^-)`")
    .def("setFromFourierOf",&GF_Bloc_ReFreq::setFromFourierOf,"Sets to the Fourier transform of Gt")
    .def("setFromPadeOf",&GF_Bloc_ReFreq::setFromPadeOf,
         (python::arg("Gw"), python::arg("N_Matsubara_Frequencies") = 100, python::arg("Freq_Offset") = .0),refreq_doc_pade)
    ;

  // **********   ImTime ******************

  class_<GF_Bloc_ImTime, bases<GF_Bloc_Base<double> > >("GFBloc_ImTime", init<object,object,PyObject *,boost::shared_ptr<MeshGF>, boost::shared_ptr<TailGF> >()) 
    .def("setFromLegendre",&GF_Bloc_ImTime::setFromLegendre,"Transforms from Legendre and set it's tail")
    .def("setFromInverseFourierOf",&GF_Bloc_ImTime::setFromInverseFourierOf,Finv_overloads("Sets to the inverse Fourier transform of Gw"))
    ;

 // **********   ReTime ******************

  class_<GF_Bloc_ReTime, bases<GF_Bloc_Base<COMPLEX> > >("GFBloc_ReTime", init<object,object,PyObject *,boost::shared_ptr<MeshGF>, boost::shared_ptr<TailGF> >()) 
    .def("density",&GF_Bloc_ReTime::density, "Computes the density :math:`G(\\tau = 0^-)`")
    .def("__call__",&GF_Bloc_ReTime::operator(), "Evaluate the function by linear interpolation $ ")
    .def("setFromInverseFourierOf",&GF_Bloc_ReTime::setFromInverseFourierOf,"Sets to the inverse Fourier transform of Gw")
    ;

  // **********   ImLegendre ******************

  char imleg_doc_enforce[] = "\
Enforces the discontinuity :math:`G(0^+)+G(\\beta^-)` by projecting the Legendre\
coefficients on the corresponding subspace. Then recompute the tail.\n\
\n\
Parameters\n\
----------\n\
A : numpy array\n\
  a matrix giving the discontinuities\n";

  class_<GF_Bloc_ImLegendre, bases<GF_Bloc_Base<double> > >("GFBloc_ImLegendre", init<object,object,PyObject *,boost::shared_ptr<MeshGF>, boost::shared_ptr<TailGF> >()) 
    .def("density",&GF_Bloc_ImLegendre::density,"Computes the density :math:`G(\\tau = 0^-)`")
    .def("determine_tail",&GF_Bloc_ImLegendre::determine_tail,"Set the tail from the Legendre coefficients")
    .def("enforce_discontinuity",&GF_Bloc_ImLegendre::enforce_discontinuity_py,imleg_doc_enforce)
    .def("setFromMatsubara",&GF_Bloc_ImLegendre::setFromImTime,"Sets from a Matsubara Green's function")
    .def("setFromMatsubara",&GF_Bloc_ImLegendre::setFromImFreq,"Sets from a Matsubara Green's function")
    ;

  //*************  GF_C converter ****************
 
   GF_C_details::register_converter<GF_Bloc_ImFreq>();
   GF_C_details::register_converter<GF_Bloc_ReFreq>();
   GF_C_details::register_converter<GF_Bloc_ImTime>();
   GF_C_details::register_converter<GF_Bloc_ReTime>();
   GF_C_details::register_converter<GF_Bloc_ImLegendre>();
 
};
