/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2012 by M. Ferrero, O. Parcollet
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
#ifndef TRIQS_MVF_MESH_H
#define TRIQS_MVF_MESH_H
#include <triqs/lazy/core.hpp>
#include <triqs/arrays/array.hpp>
#include <triqs/arrays/matrix.hpp>
#include <triqs/utility/proto/tools.hpp>
#include <triqs/arrays/h5/simple_read_write.hpp>
#include <triqs/arrays/expressions/matrix_algebra.hpp>
#include <triqs/arrays/expressions/array_algebra.hpp>
#include "./meshes.hpp"
#include "triqs/utility/complex_ops.hpp"

namespace triqs { namespace gf { 

 namespace tqa= triqs::arrays; namespace tql= triqs::lazy; namespace mpl= boost::mpl; 
 namespace proto=boost::proto; namespace bf = boost::fusion; namespace tup = triqs::utility::proto; 
 using tqa::range;

 /*------------------------------------------------------------------------------
  *-----------------------------------------------------------------------------*/

 template<typename MeshType, bool IsView, typename TailType> class mv_fnt_on_mesh;

 template<typename MeshType, bool IsView> class mv_fnt_on_mesh<MeshType,IsView,void>  { 
  friend class mv_fnt_on_mesh<MeshType,!IsView, void>;
  public:
   typedef MeshType mesh_type;
   typedef typename mesh_type::domain_type domain_type;
   
   
   typedef typename domain_type::gf_result_type element_type;


   typedef tqa::array <element_type,3,arrays::Option::Fortran>             data_non_view_type;
   typedef typename data_non_view_type::view_type                                data_view_type; 
   typedef typename mpl::if_c<IsView, data_view_type, data_non_view_type>::type  data_type;

   domain_type const & domain() const { return mymesh.domain();}
   mesh_type const & mesh() const { return mymesh;}
   data_view_type data_view()             { return data;} 
   const data_view_type data_view() const { return data;}

  protected:
   mesh_type mymesh;
   data_type data;
   
   mv_fnt_on_mesh() {} // all arrays of zero size (empty)
   mv_fnt_on_mesh (size_t N1, size_t N2, mesh_type const & m) : mymesh(m), data(N1,N2,m.size()){ data()=0;}
   mv_fnt_on_mesh (mesh_type const & m, data_view_type const & dat) : mymesh(m), data(dat){}
   mv_fnt_on_mesh(mv_fnt_on_mesh const & x): mymesh(x.mymesh), data(x.data){}// for clarity, would be synthetized anyway
   mv_fnt_on_mesh(mv_fnt_on_mesh<MeshType,!IsView,void> const & x): mymesh(x.mymesh), data(x.data){}
   // int is there to avoid any confusion with copy construction 
   template<typename GfType> mv_fnt_on_mesh(GfType const & x, int): mymesh(x.mesh()), data(x(arg_range_type())) {}
   // slice constructor. Only for view.
   template<typename A, typename B,bool V> mv_fnt_on_mesh (mv_fnt_on_mesh<MeshType, V,void> const & m, A a, B b) : 
     mymesh(m.mesh()),data(m.data(range(a),range(b), range())){static_assert(IsView, "Internal Error");}

   void operator = (mv_fnt_on_mesh const & rhs) { mymesh = rhs.mymesh; data = rhs.data; }// for clarity, would be synthetized anyway
   void operator = (mv_fnt_on_mesh<MeshType,!IsView,void> const & rhs) { mymesh = rhs.mymesh; data = rhs.data; }
   template<typename RHS> void operator = (RHS const & rhs) { // the general version  
    mymesh = rhs.mesh(); 
    BOOST_AUTO( sha , rhs(0).shape());
    tqa::resize_or_check_if_view( data, tqa::make_shape(sha[0],sha[1],mymesh.size()) );
    for (size_t u=0; u<mesh().size(); ++u) data(range(),range(),u) = rhs(mesh()[u]); 
   }

  public:
   TRIQS_NVL_HAS_AUTO_ASSIGN(); 
   template<typename F> friend void triqs_nvl_auto_assign (mv_fnt_on_mesh & x, F f) { // mesh is invariant in this case... 
    for (size_t u=0; u<x.mesh().size(); ++u) x.data(range(),range(),u) = f(x.mesh()[u]); 
   }

   typedef typename mesh_type::domain_type::point_type    arg0_type;
   typedef meshes::mesh_pt<mesh_type>                     mesh_pt_type;
   typedef typename mesh_type::range_type                 arg_range_type;

   typedef tqa::matrix_view<element_type,       arrays::Option::Fortran>  mv_type;
   typedef tqa::matrix_view<const element_type, arrays::Option::Fortran>  const_mv_type;

   mv_type       operator() (mesh_pt_type const & x)       { return data(range(),range(),x.i);}
   const_mv_type operator() (mesh_pt_type const & x) const { return data(range(),range(),x.i);}

   mv_type       operator() (arg0_type const & i)       { return mesh().interpolate(*this,i);}
   const_mv_type operator() (arg0_type const & i) const { return mesh().interpolate(*this,i);}

   data_view_type       operator() (arg_range_type const & R)       { return data(range(),range(),R);}
   const data_view_type operator() (arg_range_type const & R) const { return data(range(),range(),R);}

   /// Save the Green function in i omega_n (as 2 columns).
   void save(std::string file,  bool accumulate=false) const {}

   /// Load the GF
   void load(std::string file){}

   /// HDF5 saving .... // shall we use boost::file_system for paths ???
   // NOT OK on gcc since CommonFG is abstract... Grrr....
   //friend void h5_write( H5::CommonFG file, std::string path, mv_fnt_on_mesh const & g) {} 
   //friend void h5_read( H5::CommonFG file, std::string path, mv_fnt_on_mesh & g) {} // exceptions will propagate from the array read/write... 

 };
 /*------------------------------------------------------------------------------
  *-----------------------------------------------------------------------------*/
 template<typename MeshType, bool IsView, typename TailType> class mv_fnt_on_mesh : 
  public  mv_fnt_on_mesh<MeshType,IsView, void>  { 
   typedef mv_fnt_on_mesh<MeshType,IsView,void> B;
   friend class mv_fnt_on_mesh<MeshType,!IsView, TailType>;
   public:
   typedef typename TailType::non_view_type tail_non_view_type; 
   typedef typename TailType::view_type     tail_view_type; 
   typedef typename mpl::if_c<IsView, tail_view_type, tail_non_view_type>::type   tail_type;

   protected:
   tail_type tail;

   mv_fnt_on_mesh():B() {} // all arrays of zero size (empty)
   mv_fnt_on_mesh (size_t N1, size_t N2, typename B::mesh_type const & m) : B(N1,N2,m), tail(N1,N2,m.mesh_tail){}
   mv_fnt_on_mesh (typename B::mesh_type const & d, typename B::data_view_type const & dat, tail_view_type const & t) : B(d,dat), tail(t){}
   mv_fnt_on_mesh(mv_fnt_on_mesh const & x):                             B(x), tail(x.tail){}
   mv_fnt_on_mesh(mv_fnt_on_mesh<MeshType,!IsView, TailType> const & x): B(x), tail(x.tail){} // crucial 
   template<typename GfType> mv_fnt_on_mesh(GfType const & x, int ) : B(x,0), tail(x(domains::infty())){}
   template<typename A,typename C,bool V> // slice constructor
    mv_fnt_on_mesh (mv_fnt_on_mesh<MeshType,V,TailType> const & m,A a,C c):B(m,a,c), tail(m.tail.slice(a,c)){} 

   template<typename RHS> void operator = (RHS const & rhs) { B::operator=(rhs); tail = rhs (domains::infty()); }
   void operator = (mv_fnt_on_mesh const & rhs)             { B::operator=(rhs); tail = rhs.tail;} 

   public:
   TRIQS_NVL_HAS_AUTO_ASSIGN(); 
   template<typename F> friend void triqs_nvl_auto_assign (mv_fnt_on_mesh & x, F f) { // mesh is invariant in this case... 
    triqs_nvl_auto_assign( static_cast<B &>(x), f); // calling as a the base class first ....
    x.tail =  f (tail_non_view_type(x.data.shape()[0], x.data.shape()[1],x.mesh().mesh_tail ));
    // if f is an expression, replace the placeholder with a simple tail. If f is a function callable on infty, 
    // it uses the fact that tail_non_view_type can be caster into domains::infty (See below).
   }

   using B::operator();
   typename B::mv_type       operator() (typename B::arg0_type const & i)       {return this->mesh().interpolate(*this,i);}
   typename B::const_mv_type operator() (typename B::arg0_type const & i) const {return this->mesh().interpolate(*this,i);}

   tail_view_type       operator() ( domains::infty const & x)       {return tail;}
   const tail_view_type operator() ( domains::infty const & x) const {return tail;}

   /// Save the Green function in i omega_n (as 2 columns).
   void save(std::string file,  bool accumulate=false) const {}

   /// Load the GF
   void load(std::string file){}

   /// HDF5 saving .... // shall we use boost::file_system for paths ???
   // NOT OK on gcc since CommonFG is abstract... Grrr....
   //friend void h5_write( H5::CommonFG file, std::string path, mv_fnt_on_mesh const & g) {} 
   //friend void h5_read( H5::CommonFG file, std::string path, mv_fnt_on_mesh & g) {} // exceptions will propagate from the array read/write... 

  };
}}
#endif
