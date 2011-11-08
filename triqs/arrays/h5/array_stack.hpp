/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011 by O. Parcollet
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
#ifndef TRIQS_ARRAYS_H5_STACK_H
#define TRIQS_ARRAYS_H5_STACK_H
#include "../array.hpp"
#include "./common.hpp"

namespace triqs { namespace arrays { namespace h5 { 
 using namespace H5;

 namespace details { // to be replaced by ellipsis 
  template<class T, size_t N, class Opt> array_view<T,N-1,Opt> slice0( array<T,N,Opt> const & A, size_t ind);
#define AUX(z,p,unused) BOOST_PP_COMMA_IF(p) range()
#define IMPL(z, NN, unused) \
  template<class T, class Opt> array_view<T,BOOST_PP_INC(NN),Opt> slice0( array<T,BOOST_PP_INC(NN)+1,Opt> const & A, size_t ind) {\
   return A(ind,BOOST_PP_REPEAT(BOOST_PP_INC(NN),AUX,nil));}
  BOOST_PP_REPEAT(ARRAY_NRANK_MAX , IMPL, nil);
#undef IMPL
#undef AUX
 }

 /**
  *  Hdf5 array stack 
  */
 template<class T, size_t dim>
  class array_stack {
   size_t bufsize, step, _size; 
   static const bool T_is_complex = boost::is_complex<T>::value;
   static const unsigned int RANK = dim + 1 + (T_is_complex ? 1 : 0);
   mini_vector<hsize_t,RANK> dims, offset, maxdims, dim_chunk, buffer_dim, zero;
   DataSet dataset;
   array<T,dim+1> buffer;

   public :

   size_t size() const { return _size;}

   template <typename FileGroupType >
    array_stack( FileGroupType file_or_group, std::string const & name, mini_vector<size_t,dim> const & a_dims, size_t bufsize_)  {
     mini_vector<hsize_t,RANK> dim_chunk;
     bufsize = bufsize_; step = 0; _size =0; 
     for (size_t i =1; i<=dim; ++i) { dims[i] = a_dims[i-1];}
     if (T_is_complex) { dims[RANK-1] =2; }
     maxdims = dims; buffer_dim = dims; dim_chunk = dims;
     dims[0] = 0; maxdims[0] = H5S_UNLIMITED; dim_chunk[0]=1; buffer_dim[0] = bufsize;
     mini_vector<size_t,dim+1> s; for (size_t i =0; i<=dim; ++i) {s[i] =  buffer_dim[i];} 
     buffer.resize(s);
     DataSpace mspace1( RANK, dims.ptr(), maxdims.ptr());
     DSetCreatPropList cparms; cparms.setChunk( RANK, dim_chunk.ptr() ); // Modify dataset creation properties, i.e. enable chunking.
     try { 
      if (h5::exists(file_or_group, name.c_str())) file_or_group.unlink( name.c_str());  
      dataset = file_or_group.createDataSet( name.c_str(), native_type_from_C(typename remove_complex<T>::type()), mspace1, cparms );
      if (boost::is_complex<T>::value)  write_attribute(dataset,"__complex__","1");
     }
     TRIQS_ARRAYS_H5_CATCH_EXCEPTION;
    }

   ~array_stack() {flush();} 
   array_view<T,dim> operator() () { return details::slice0(buffer, step); } 
   void operator++() { ++step; ++_size; if (step==bufsize) flush();  } 
   void flush() { save_buffer(); step=0;}
   template<class AType> void operator << ( AType const & A) { (*this)() = A; ++(*this);}

   protected:
   void save_buffer () {
    if (step==0) return;
    dims[0] += step;
    buffer_dim[0] = step; 
    dataset.extend(dims.ptr());
    DataSpace fspace1 = dataset.getSpace (), mspace = data_space(buffer); 
    fspace1.selectHyperslab( H5S_SELECT_SET, buffer_dim.ptr(), offset.ptr() );
    mspace.selectHyperslab(  H5S_SELECT_SET, buffer_dim.ptr(), zero.ptr() );
    try { dataset.write( data(buffer), data_type_mem(buffer), mspace, fspace1 ); }
    TRIQS_ARRAYS_H5_CATCH_EXCEPTION;
    offset [0] += step;
   }
  };
}}} // namespace
#endif

