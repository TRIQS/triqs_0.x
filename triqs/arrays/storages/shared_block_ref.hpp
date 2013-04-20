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
#ifndef TRIQS_STORAGE_SHARED_POINTER_REF_H
#define TRIQS_STORAGE_SHARED_POINTER_REF_H
#include "./shared_block.hpp"

namespace triqs { namespace arrays { 
 namespace storages {

  template<typename ValueType>
   class shared_block_ref { 

    typedef typename boost::add_const<ValueType>::type const_value_type;
    typedef typename boost::remove_const<ValueType>::type non_const_value_type;

    public:
    typedef ValueType value_type;
    typedef shared_block<value_type> clone_type;

    explicit shared_block_ref() { sbl = nullptr; data_ = nullptr;}

    shared_block_ref(const shared_block<non_const_value_type> & X)    : sbl(&X) { init_data(); }
    shared_block_ref(const shared_block<const_value_type> & X)    : sbl(&X) { init_data(); }
    //shared_block_ref(const shared_block_ref & X): sbl(X.sbl) { init_data(); }
    
    shared_block_ref(shared_block_ref const & X) = default;
    shared_block_ref (shared_block_ref &&) = default;
    shared_block_ref & operator=(const shared_block_ref & X) = default;

    clone_type clone() const { return sbl->clone();}

    value_type & operator[](size_t i) const { return data_[i];}
    bool empty() const {return sbl->empty();}
    size_t size() const {return sbl->size();}
 
    shared_block<value_type> & get () { return *sbl;}
    shared_block<value_type> const & get () const { return *sbl;}

    operator shared_block<value_type> const & () const { return *sbl;}

    private:
    const shared_block<value_type> * sbl;
    value_type * restrict data_; // for optimization on some compilers.
    void init_data(){ data_ = (sbl != nullptr ? &((*sbl)[0]) : nullptr); }
   };
 }

 namespace details { 
  template<bool Const, typename T> struct make_const_type<Const,storages::shared_block_ref<T> > { 
   typedef storages::shared_block_ref<typename make_const_type<Const,T>::type> type;
  };
 }
}}//namespace triqs::arrays 
#endif

