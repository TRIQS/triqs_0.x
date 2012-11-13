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
#ifndef TRIQS_ARRAY_IMPL_INDEX_STORAGE_PAIR_H
#define TRIQS_ARRAY_IMPL_INDEX_STORAGE_PAIR_H

#include "./common.hpp"
#include "../storages/shared_block.hpp"
#include "./assignment.hpp"
#include "./option.hpp"
#include "./sliceable_object.hpp"
#include "triqs/utility/exceptions.hpp"
#include "triqs/utility/typeid_name.hpp"
#include "triqs/utility/view_tools.hpp"

#include <boost/type_traits/add_const.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/if.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>

#include <type_traits>
#ifdef TRIQS_WITH_PYTHON_SUPPORT
#include "../python/numpy_extractor.hpp"
#endif

namespace triqs { namespace arrays { 

 template <bool Const, typename IndexMapIterator, typename StorageType > class iterator_adapter;

 template < class V, int R, class Opt, class ViewTag > struct ViewFactory;

 template <typename IndexMapType, typename StorageType, typename Opt,  typename ViewTag > 
  class indexmap_storage_pair : Tag::indexmap_storage_pair, TRIQS_MODEL_CONCEPT(ImmutableArray),
  public sliceable_object < typename StorageType::value_type,
  IndexMapType,
  Opt,
  ViewTag, 
  ViewFactory,
  indexmaps::slicer,
  indexmap_storage_pair<IndexMapType,StorageType,Opt,ViewTag> > {    

   public : 
    typedef typename StorageType::value_type value_type;
    //static_assert(std::is_default_constructible<value_type>::value, "array/array_view and const operate only on values");
    typedef StorageType storage_type;
    typedef IndexMapType indexmap_type;
    static const unsigned int rank = IndexMapType::domain_type::rank;

   protected:
    indexmap_type indexmap_;
    storage_type storage_;

    indexmap_storage_pair() {}

    indexmap_storage_pair (const indexmap_type & IM, const storage_type & ST):
     indexmap_(IM),storage_(ST){
#ifdef TRIQS_ARRAYS_CHECK_IM_STORAGE_COMPAT
      if (ST.size() != IM.domain().number_of_elements()) 
       TRIQS_RUNTIME_ERROR<<"index_storage_pair construction : storage and indices are not compatible";
#endif
     }

    /// The storage is allocated from the size of IM.
    indexmap_storage_pair (const indexmap_type & IM): indexmap_(IM),storage_(){
     this->storage_ = StorageType(this->indexmap_.domain().number_of_elements(), typename Opt::InitTag() );
    }

    /// Shallow copy
    indexmap_storage_pair(const indexmap_storage_pair & X):indexmap_(X.indexmap()),storage_(X.storage_){}

#ifdef TRIQS_WITH_PYTHON_SUPPORT
    indexmap_storage_pair (PyObject * X, bool allow_copy, const char * name ) { 
     try { 
      numpy_interface::numpy_extractor<indexmap_type,value_type> E(X, allow_copy);
      this->indexmap_ = E.indexmap(); this->storage_  = E.storage();
     }
     catch(numpy_interface::copy_exception s){// intercept only this one...
      TRIQS_RUNTIME_ERROR<< " construction of a "<< name <<" from a numpy  "
       <<"\n   T = "<< triqs::utility::typeid_name(value_type())
       //  <<"\n   rank = "<< IndexMapType::domain_type::rank//this->rank
       <<"\n   Opt = "<< triqs::utility::typeid_name(Opt())
       <<"\nfrom the python object \n"<< numpy_interface::object_to_string(X) 
       <<"\nThe error was :\n "<<s.what();
     }
    }

#endif

     void swap_me( indexmap_storage_pair & X) {
     std::swap(this->indexmap_,X.indexmap_); std::swap (this->storage_, X.storage_);
    }

    friend void swap( indexmap_storage_pair & A, indexmap_storage_pair & B) { A.swap_me(B);}

   public:

    indexmap_type const & indexmap() const {return indexmap_;}
    storage_type const & storage() const {return storage_;}
    storage_type & storage() {return storage_;}

    /// data_start is the starting point of the data of the object
    /// this it NOT &storage()[0], which is the start of the underlying blokc
    /// they are not equal for a view in general
    value_type const * restrict data_start() const { return &storage_[indexmap_.start_shift()];}
    value_type * restrict data_start() { return &storage_[indexmap_.start_shift()];}

    typedef typename indexmap_type::domain_type domain_type; 
    domain_type const & domain() const { return indexmap_.domain();}

    typedef typename domain_type::index_value_type shape_type;
    shape_type const & shape() const { return domain().lengths();}

    size_t len(size_t i) const { return this->shape()[i]; }

    size_t num_elements() const { return domain().number_of_elements();}

    //bool is_empty() const { return this->num_elements()==0;}
    bool is_empty() const { return this->storage_.empty();}

    //  Evaluation. Slices are made by Sliceable object 
    typedef typename domain_type::index_value_type key_type;
    template<typename KeyType> value_type const & operator[](KeyType const & key) const { return storage_[indexmap_[key]]; } 
    template<typename KeyType> value_type & operator[](KeyType const & key) { return storage_[indexmap_[key]]; } 

    // Iterators
    typedef iterator_adapter<true,typename IndexMapType::iterator, StorageType> const_iterator;
    typedef iterator_adapter<false,typename IndexMapType::iterator, StorageType> iterator;
    const_iterator begin() const {return const_iterator(indexmap(),storage(),false);}
    const_iterator end() const {return const_iterator(indexmap(),storage(),true);}
    iterator begin() {return iterator(indexmap(),storage(),false);}
    iterator end() {return iterator(indexmap(),storage(),true);}

   protected:

    void resize (domain_type const & d) {
     this->indexmap_ = IndexMapType(d);// build a new one with the lengths of IND
     // optimisation. Construct a storage only if the new index is not compatible (size mismatch).
     if (this->storage_.size() != this->indexmap_.domain().number_of_elements())
      this->storage_ = StorageType(this->indexmap_.domain().number_of_elements(), typename Opt::InitTag() );
    }

    template<typename Xtype>
     void resize_and_clone_data( Xtype const & X) { indexmap_ = X.indexmap(); storage_ = X.storage().clone(); }

    //  BOOST Serialization
    friend class boost::serialization::access;
    template<class Archive>
     void serialize(Archive & ar, const unsigned int version) {
      ar & boost::serialization::make_nvp("storage",this->storage_);
      ar & boost::serialization::make_nvp("indexmap",this->indexmap_);
     }
  };// end class

 // pretty print of the array
 template <typename I, typename S ,class Opt, typename V> 
  std::ostream & operator << (std::ostream & out, const triqs::arrays::indexmap_storage_pair<I,S,Opt,V> & A) {
   //std::cerr<< "Lengths = "<<A.indexmap().lengths()<<"Strides = "<<A.indexmap().strides()<< "  ";
   if (A.storage().size()==0) out<<"empty ";
   else pretty_print(out, A.indexmap().domain(),A);
   return out;
  }

}}//namespace triqs::arrays 
#endif

