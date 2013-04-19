/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011-2013 by O. Parcollet
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
#include "./flags.hpp"
#include "../storages/shared_block.hpp"
#include "./assignment.hpp"
#include "../indexmaps/cuboid/foreach.hpp"
#include "triqs/utility/exceptions.hpp"
#include "triqs/utility/typeid_name.hpp"
#include "triqs/utility/view_tools.hpp"
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>
#include <type_traits>
#ifdef TRIQS_WITH_PYTHON_SUPPORT
#include "../python/numpy_extractor.hpp"
#include "../python/array_view_to_python.hpp"
#endif

namespace triqs { namespace arrays {

 template <bool Const, typename IndexMapIterator, typename StorageType > class iterator_adapter;

 template < class V, int R, ull_t OptionFlags, ull_t TraversalOrder,  class ViewTag > struct ViewFactory;

 template <typename IndexMapType, typename StorageType, ull_t OptionFlags, ull_t TraversalOrder,  typename ViewTag >
  class indexmap_storage_pair : Tag::indexmap_storage_pair, TRIQS_MODEL_CONCEPT(MutableArray) {

   public :
    typedef typename StorageType::value_type value_type;
    static_assert(std::is_constructible<value_type>::value, "array/array_view and const operate only on values");
    //static_assert(!std::is_const<value_type>::value,         "no const type");
    typedef StorageType storage_type;
    typedef IndexMapType indexmap_type;
    static constexpr unsigned int rank = IndexMapType::domain_type::rank;
   protected:

    indexmap_type indexmap_;
    storage_type storage_;

    // ------------------------------- constructors --------------------------------------------

    indexmap_storage_pair() {}

    indexmap_storage_pair (const indexmap_type & IM, const storage_type & ST):
     indexmap_(IM),storage_(ST){
#ifdef TRIQS_ARRAYS_CHECK_IM_STORAGE_COMPAT
      if (ST.size() != IM.domain().number_of_elements())
       TRIQS_RUNTIME_ERROR<<"index_storage_pair construction : storage and indices are not compatible";
#endif
     }

    indexmap_storage_pair (indexmap_type && IM, storage_type && ST):
     indexmap_(std::move(IM)),storage_(std::move(ST)){
#ifdef TRIQS_ARRAYS_CHECK_IM_STORAGE_COMPAT
      if (ST.size() != IM.domain().number_of_elements())
       TRIQS_RUNTIME_ERROR<<"index_storage_pair construction : storage and indices are not compatible";
#endif
     }


    /// The storage is allocated from the size of IM.
    indexmap_storage_pair (const indexmap_type & IM): indexmap_(IM),storage_(){
     this->storage_ = StorageType(this->indexmap_.domain().number_of_elements(), typename flags::init_tag<OptionFlags>::type() );
    }

   public:
    /// Shallow copy
    indexmap_storage_pair(const indexmap_storage_pair & X):indexmap_(X.indexmap()),storage_(X.storage_){}
    indexmap_storage_pair(indexmap_storage_pair && X):indexmap_(std::move(X.indexmap())),storage_(std::move(X.storage_)){}
   protected:

#ifdef TRIQS_WITH_PYTHON_SUPPORT
    indexmap_storage_pair (PyObject * X, bool allow_copy, const char * name ) {
     //std::cout << " Enter IPS ref count = "<< X->ob_refcnt << std::endl;
     try {
      numpy_interface::numpy_extractor<indexmap_type,value_type> E(X, allow_copy);
      this->indexmap_ = E.indexmap(); this->storage_  = E.storage();
     }
     catch(numpy_interface::copy_exception s){// intercept only this one...
      TRIQS_RUNTIME_ERROR<< " construction of a "<< name <<" from a numpy  "
       <<"\n   T = "<< triqs::utility::typeid_name(value_type())
       // lead to a nasty link pb ???
       // linker search for IndexMapType::domain_type::rank in triqs.so
       // and can not resolve it ???
       //<<"\n   rank = "<< IndexMapType::domain_type::rank//this->rank
       <<"\n   OptionFlags = "<< OptionFlags
       <<"\nfrom the python object \n"<< numpy_interface::object_to_string(X)
       <<"\nThe error was :\n "<<s.what();
     }
     //std::cout << " Leave IPS ref count = "<< X->ob_refcnt << std::endl;
    }
#endif
    // ------------------------------- swap --------------------------------------------

    void swap_me( indexmap_storage_pair & X) {
     std::swap(this->indexmap_,X.indexmap_); std::swap (this->storage_, X.storage_);
    }

    friend void swap( indexmap_storage_pair & A, indexmap_storage_pair & B) { A.swap_me(B);}

    // ------------------------------- ==  --------------------------------------------
    // at your own risk with floating value, but it is useful for int, string, etc....
    // in particular for tests
    friend bool operator==( indexmap_storage_pair const & A, indexmap_storage_pair const & B) {
     if (A.shape() != B.shape()) return false;
     auto ita = A.begin(); auto itb = B.begin();
     for (;ita != A.end();++ita, ++itb) {if (!(*ita == *itb)) return false;}
     return true;
    }

    friend bool operator!=( indexmap_storage_pair const & A, indexmap_storage_pair const & B) { return (!(A==B));}

   public:
    // ------------------------------- data access --------------------------------------------

    indexmap_type const & indexmap() const {return indexmap_;}
    storage_type const & storage() const {return storage_;}
    storage_type & storage() {return storage_;}

#ifdef TRIQS_WITH_PYTHON_SUPPORT
    PyObject * to_python() const { return numpy_interface::array_view_to_python(*this);}
#endif

    /// data_start is the starting point of the data of the object
    /// this it NOT &storage()[0], which is the start of the underlying blokc
    /// they are not equal for a view in general
    value_type const * restrict data_start() const { return &storage_[indexmap_.start_shift()];}
    value_type * restrict data_start() { return &storage_[indexmap_.start_shift()];}

    typedef typename indexmap_type::domain_type domain_type;
    domain_type const & domain() const { return indexmap_.domain();}

    typedef typename domain_type::index_value_type shape_type;
    shape_type const & shape() const { return domain().lengths();}

    size_t shape(size_t i) const { return domain().lengths()[i];}

    size_t len(size_t i) const { return this->shape()[i]; }

    size_t num_elements() const { return domain().number_of_elements();}

    //bool is_empty() const { return this->num_elements()==0;}
    bool is_empty() const { return this->storage_.empty();}

#ifdef TRIQS_ARRAYS_USE_OLD_KEY_EVAL
    //  Evaluation on keys.
    // deprecated : only used in one example ?
    typedef typename domain_type::index_value_type key_type;
    template<typename KeyType> value_type const & operator[](KeyType const & key) const { return storage_[indexmap_[key]]; }
    template<typename KeyType> value_type & operator[](KeyType const & key) { return storage_[indexmap_[key]]; }
#endif

    // ------------------------------- operator () --------------------------------------------

    // Evaluation and slices
    template<typename... Args>
     typename std::enable_if<
     (!clef::is_any_lazy<Args...>::value) && (indexmaps::slicer<indexmap_type,Args...>::r_type::domain_type::rank==0)
     , value_type &>::type
     operator()(Args const & ... args) {  return storage_[indexmap_(args...)]; }

    template<typename... Args>
     typename std::enable_if<
     (!clef::is_any_lazy<Args...>::value) && (indexmaps::slicer<indexmap_type,Args...>::r_type::domain_type::rank==0)
     , value_type const &>::type
     operator()(Args const & ... args) const { return storage_[indexmap_(args...)]; }

    template<bool is_const, typename ... Args> struct result_of_call_as_view {
     typedef typename indexmaps::slicer<indexmap_type,Args...>::r_type IM2;
     typedef typename std::conditional<is_const, typename std::add_const<value_type>::type, value_type>::type V2;
     static_assert(IM2::domain_type::rank !=0, "Internal error");
     typedef typename ViewFactory<V2,IM2::domain_type::rank, OptionFlags, IM2::traversal_order_in_template,  ViewTag >::type type;
    };

    template<typename... Args>   // non const version
     typename boost::lazy_enable_if_c<
     (!clef::is_any_lazy<Args...>::value) && (indexmaps::slicer<indexmap_type,Args...>::r_type::domain_type::rank!=0)
     , result_of_call_as_view<false,Args...>
     >::type // enable_if
     operator()(Args const & ... args) {
      return typename result_of_call_as_view<false,Args...>::type ( indexmaps::slicer<indexmap_type,Args...>::invoke(indexmap_,args...), storage()); }

    template<typename... Args>  // const version
     typename boost::lazy_enable_if_c<
     (!clef::is_any_lazy<Args...>::value) && (indexmaps::slicer<indexmap_type,Args...>::r_type::domain_type::rank!=0)
     , result_of_call_as_view<true,Args...>
     >::type // enable_if
     operator()(Args const & ... args) const {
      return typename result_of_call_as_view<true,Args...>::type ( indexmaps::slicer<indexmap_type,Args...>::invoke(indexmap_,args...), storage()); }

    typedef typename ViewFactory<value_type,domain_type::rank, OptionFlags, TraversalOrder, ViewTag >::type view_type;

    // Interaction with the CLEF library : calling with any clef expression as argument build a new clef expression
    // NB : this is ok because indexmap_storage_pair has a shallow copy constructor ...
    // so A(i_) if A is an array will NOT copy the data....
    template< typename... Args>
     typename clef::result_of::make_expr_call<indexmap_storage_pair,Args...>::type
     operator()( Args&&... args ) const {
      static_assert(sizeof...(Args) <= indexmap_type::rank, "Incorrect number of variable in call");// not perfect : ellipsis ...
      return make_expr_call(*this,args...);
      //return make_expr_call( view_type(*this),args...);
     }

    template<typename Fnt> friend void triqs_clef_auto_assign (indexmap_storage_pair & x, Fnt f) { assign_foreach(x,f);}

    /// Equivalent to make_view
    typename ViewFactory<typename std::add_const<value_type>::type,domain_type::rank, OptionFlags, TraversalOrder, ViewTag >::type
     operator()() const { return *this; }
    typename ViewFactory<value_type,domain_type::rank, OptionFlags, TraversalOrder, ViewTag >::type
     operator()() { return *this; }

    // ------------------------------- Iterators --------------------------------------------
    typedef iterator_adapter<true,typename IndexMapType::iterator, StorageType> const_iterator;
    typedef iterator_adapter<false,typename IndexMapType::iterator, StorageType> iterator;
    const_iterator begin() const {return const_iterator(indexmap(),storage(),false);}
    const_iterator end() const {return const_iterator(indexmap(),storage(),true);}
    iterator begin() {return iterator(indexmap(),storage(),false);}
    iterator end() {return iterator(indexmap(),storage(),true);}

   protected:

    // ------------------------------- resize --------------------------------------------
    //
    void resize (domain_type const & d) {
     this->indexmap_ = IndexMapType(d);// build a new one with the lengths of IND
     // optimisation. Construct a storage only if the new index is not compatible (size mismatch).
     if (this->storage_.size() != this->indexmap_.domain().number_of_elements())
      this->storage_ = StorageType(this->indexmap_.domain().number_of_elements(), typename flags::init_tag<OptionFlags>::type()  );
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

    // pretty print of the array
    friend std::ostream & operator << (std::ostream & out, const indexmap_storage_pair & A) {
     if (A.storage().size()==0) out<<"empty ";
     else indexmaps::pretty_print(out,A);
     return out;
    }
  };// end class

}}//namespace triqs::arrays
#endif

