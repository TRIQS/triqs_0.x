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
#ifndef TRIQS_ARRAYS_HDF5_ARRAY_PROXY_H
#define TRIQS_ARRAYS_HDF5_ARRAY_PROXY_H
//#define TRIQS_ARRAYS_DEBUG_H5_SLICE

#include "../indexmaps/cuboid/cuboid_domain.hpp"
#include "./common.hpp"
#include "../impl/sliceable_object.hpp"
#include "../impl/tuple_tools.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>

namespace triqs { namespace arrays { 

 namespace h5 { template<int Rf> struct array_proxy_option {}; } 

 namespace Option { 
  template<class IM2, int Rf> struct Im_Opt_2_Opt<IM2, h5::array_proxy_option<Rf> > { typedef h5::array_proxy_option<Rf> type;}; 
 }

 namespace h5 {

  template<int Rank, int Rank_full>
   class index_system { 

    public :
     static const unsigned int rank= Rank, rank_full = Rank_full;
     typedef mini_vector<size_t, Rank>      v_type;
     typedef mini_vector<size_t, Rank_full> v_type_full;
     typedef indexmaps::cuboid_domain<Rank> domain_type;
     domain_type const & domain() const { return mydomain;}

     index_system (v_type_full const & total_lengths_) {
      total_lens_ = total_lengths_;  lens_= total_lengths_;
      for (size_t i =0; i<rank_full; ++i)  stri_[i] = 1; 
      for (size_t i =0; i<rank; ++i)  dims[i] = lens_[i]; 
      mydomain = domain_type (dims);
     }

     index_system (v_type const & dims_, v_type_full const & total_lengths_,  v_type_full const & lengths_, 
       v_type_full const & strides_, v_type_full const & offset_ ) {
      dims = dims_;total_lens_ = total_lengths_;  lens_= lengths_; stri_ = strides_; off_=offset_;
      mydomain = domain_type (dims);
     }

     index_system (DataSpace const & ds, bool is_complex) { 
      int rf = ds.getSimpleExtentNdims();
      if ( rf != rank_full  + (is_complex ? 1 : 0) ) TRIQS_RUNTIME_ERROR <<  "H5 : dimension error";
      int ndims = ds.getSimpleExtentDims( &lens_[0], NULL);
      for (size_t i =0; i<rank; ++i) { dims[i] = lens_[i]; stri_[i] = 1; off_[i]= 0; }
      total_lens_=dims;
      mydomain = domain_type (dims);
     }

     template<bool C> DataSpace dataspace() const { return h5::dataspace_from_LS<rank_full,C> (total_lens_, lens_,stri_,off_); }
     size_t size() const { size_t _size = 1; for (size_t i =0; i<rank; ++i) _size *= lens_[i]; return _size;}
     mini_vector<size_t, rank_full> total_lengths() const { return this->total_lens_;}
     mini_vector<size_t, rank_full> lengths() const { return this->lens_;}
     mini_vector<size_t, rank_full> strides() const { return this->stri_;}

    private: 
     mini_vector<hsize_t,rank_full> lens_, off_, stri_; 
     v_type dims;
     v_type_full total_lens_;
     domain_type mydomain;
   };

  namespace slicer_impl { 

   template <bool BC> inline void _check_BC ( int N, int ind, size_t B) { }
   template <> inline void _check_BC<true> (int N, int ind, size_t B) { 
    bool cond = (ind >= 0) && (ind < B);
    if (!cond) TRIQS_ARRAYS_KEY_ERROR << " index "<<N<<" is out of domain: \n    " << ind <<" is not within [0,"<< B <<"[\n";
   }

   template<int Rank_in, int Rank_out, int N,int P, int c, bool BoundCheck> struct slice_calc { 

    typedef mini_vector<size_t,Rank_in> const &  i_type;
    typedef mini_vector<size_t,Rank_in > &       o_type;
    typedef mini_vector<size_t,Rank_out > &      os_type; 

    template< typename ArgsTuple>
     static void invoke(i_type li, i_type si, os_type lo_c, o_type lo,  o_type so, o_type offset, ArgsTuple const & args ) {
      const int dP = boost::is_base_of<typename boost::tuples::element<N,ArgsTuple>::type, range >::type::value ;
      one_step(li[N], si[N],lo[N],so[N], offset[N] ,boost::tuples::get<N>(args));
      lo_c[P] = lo[N];
      slice_calc<Rank_in,Rank_out,N+1,P+dP,c-1, BoundCheck>::invoke(li,si,lo_c,lo,so, offset, args);
     }

    static void one_step(size_t li, size_t si, size_t & lo, size_t & so, size_t & offset, size_t R){
     _check_BC <BoundCheck> (N, R, li);
     offset = R;  lo =1; so= 1;
    }

    static void one_step(size_t li, size_t si, size_t & lo, size_t & so, size_t & offset, range R){
     _check_BC <BoundCheck> (N, R.first(),li);
     _check_BC <BoundCheck> (N, (R.last()==-1 ? li : R.last()) -1 ,li);
     lo  = ((R.last()==-1 ? li : R.last()) - R.first() + R.step()-1 )/R.step() ; //  python behaviour
     so  = R.step();  offset = R.first()  ;
    }
   };
   // stop the recursion
   template<int Ri, int Ro, int N, int P, bool BC> struct slice_calc <Ri,Ro,N,P,0,BC> : slice_calc<Ri,Ro,N,P,1,BC> {
    template<class T1,class T2,class T3,class T4,class T5,class T6, class T7> static void invoke(T1,T2,T3,T4,T5,T6,T7 ) {}
   };

  }//namespace slicer_impl 
 } // h5

 namespace indexmaps { 

  template<int R, int Rf, typename ArgsTuple> 
   struct slicer < h5::index_system<R,Rf>,  ArgsTuple>  { 

    static const unsigned int len = boost::tuples::length<ArgsTuple>::value;
    static_assert(len>=R, "Too few arguments in slice"); 
    static_assert(len<=R, "Too many arguments in slice");
    static_assert( (R==Rf), "Can not slice array_proxy twice (not implemented)");

    static const unsigned int R2 = R - TupleTools::CountHowManyInt<ArgsTuple>::value;
    typedef h5::index_system< R2, Rf> return_type; 

    static return_type invoke ( h5::index_system<R,R> const & X, ArgsTuple args) { 
     mini_vector<size_t, R2> newdims;
     mini_vector<size_t,R> newoffset, newlengths, newstrides;
     h5::slicer_impl::slice_calc<R,R2,0,0,R,true>::invoke(X.lengths(),X.strides(),newdims, newlengths,newstrides, newoffset, args);
#ifdef TRIQS_ARRAYS_DEBUG_H5_SLICE
     std::cerr<<"-----------------------------------------------"<<std::endl;
     std::cerr<<"H5 Slicing "<< X.lengths().to_string()<<X.strides().to_string()<<newlengths.to_string()<<newstrides.to_string()<< newoffset.to_string() << args<<std::endl;
#endif
     return return_type(newdims,X.total_lengths(), newlengths,newstrides,newoffset);
    };
   }; 

  template<int R, int Rf> struct slicer < h5::index_system<R,Rf>, boost::tuple<> >  { typedef h5::index_system< R,Rf>  return_type;}; 
 }

 namespace h5 { // -----------------------------------------------------------------------------

  /// The array proxy
  template<typename ValueType, int Rank, int Rank_f = Rank >
   class array_proxy :
    Tag::has_special_assign, //Tag::has_immutable_array_interface,
    public sliceable_object <
    ValueType,
    h5::index_system<Rank,Rank_f>, 
    array_proxy_option<Rank_f>,
    Tag::h5_array_proxy, 
    ViewFactory,
    indexmaps::slicer,
    array_proxy<ValueType,Rank, Rank_f>
    >
  {
   public : 
    typedef ValueType value_type;
    typedef std::pair< boost::shared_ptr<H5::CommonFG> ,std::string> storage_type; 
    static const bool T_is_complex = boost::is_complex<ValueType>::value;
    typedef index_system< Rank, Rank_f> indexmap_type;
    static const unsigned int rank = Rank;

    /// Opens a proxy on an existing array. The dataset must exists
    template< class FileGroupType >
     array_proxy ( FileGroupType file_group, std::string const & name) :
      indexmap_ ( indexmap_type(file_group.openDataSet( name.c_str() ).getSpace(),T_is_complex) ) { 
       if (!h5::exists(file_group, name)) TRIQS_RUNTIME_ERROR<< " h5 : no dataset"<< name << " in file "; 
       storage_ = std::make_pair( boost::make_shared<FileGroupType>(file_group),name);
       DataSet dataset = file_group.openDataSet( name.c_str() );
       try { if (T_is_complex) write_attribute(dataset,"__complex__","1"); } 
       catch (...) {} // catch if the attribute already exists...
      }

    /// Constructs a proxy on a new data set of the dimension of the domain D.  The data must not exist. 
    template< class FileGroupType, class LengthType >
     array_proxy ( FileGroupType file_group, std::string const & name_, LengthType L, bool overwrite = false) :
      indexmap_ ( indexmap_type (L) ) 
   { 
    if (h5::exists(file_group, name_)) {
     if (overwrite) file_group.unlink(name_.c_str());  
     else TRIQS_RUNTIME_ERROR<< " h5 : dataset"<< name_ << " already exists in the file "; 
    }
    storage_ = std::make_pair( boost::make_shared<FileGroupType>(file_group),name_);
    DataSpace ds  = indexmap_.template dataspace<T_is_complex>(); //(indexmap_type::rank_full, &indexmap_.lengths()[0], &indexmap_.strides()[0]  );
    DataSet dataset = file_group.createDataSet( name_.c_str(), data_type_file(ValueType()), ds);
    if (T_is_complex) write_attribute(dataset,"__complex__","1");
   }

    /// Shallow copy
    array_proxy(const array_proxy & X):indexmap_(X.indexmap()),storage_(X.storage_){}

    /// for slice construction
    array_proxy (const indexmap_type & IM, const storage_type & ST): indexmap_(IM),storage_(ST){ }

   public:
    indexmap_type const & indexmap() const {return indexmap_;}
    storage_type const & storage() const {return storage_;}
    storage_type & storage() {return storage_;}
    const H5::CommonFG * file_group() const { return storage_.first.get();}
    std::string const & name() const { return storage_.second;}

    typedef typename indexmap_type::domain_type domain_type; 
    domain_type const & domain() const { return indexmap_.domain();}

    typedef typename domain_type::index_value_type shape_type;
    shape_type const & shape() const { return domain().lengths();}

    size_t num_elements() const { return domain().number_of_elements();}
    bool is_empty() const { return this->num_elements()==0;}

    template< typename ISP>// put in the file...
     void operator=(ISP const & X) { 
      try { 
       BOOST_AUTO(C,  make_const_cache_C_order(X));
       //typename result_of::cache<false,Tag::C, ISP >::type C(X);
       DataSet dataset = file_group()->openDataSet( name().c_str() );
       dataset.write( h5::data(C.view()), h5::data_type_mem(C.view()),h5::data_space(C.view()),indexmap().template dataspace<T_is_complex>()); 
      } 
      TRIQS_ARRAYS_H5_CATCH_EXCEPTION;
     }

    template<typename LHS> // from the file to the array or the array_view... 
     void assign_invoke (LHS & lhs) const { 
      static_assert((is_value_or_view_class<LHS>::value), "LHS is not a value or a view class "); 
      h5::resize_or_check(lhs,  indexmap().domain().lengths());
      try { 
       DataSet dataset = file_group()->openDataSet( name().c_str() );
       BOOST_AUTO(C,  make_cache_C_order(lhs));
       //typename result_of::cache<true,Tag::C, LHS >::type C(lhs);
       dataset.read( h5::data(C.view()), h5::data_type_mem(C.view()), h5::data_space(C.view()),
	 indexmap().template dataspace<T_is_complex>() ); 
      }
      TRIQS_ARRAYS_H5_CATCH_EXCEPTION;
     }
   protected:
    indexmap_type indexmap_;
    storage_type storage_;

  };
 }

 template < class V, int R, int Rf > 
  struct ViewFactory< V, R, h5::array_proxy_option<Rf>, Tag::h5_array_proxy > { typedef h5::array_proxy <V,R,Rf> type; };
}}
#endif

