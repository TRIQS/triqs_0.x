/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011-2012 by M. Ferrero, O. Parcollet
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
#ifndef TRIQS_DETMANIP_CTHYB_H
#define TRIQS_DETMANIP_CTHYB_H

#include <triqs/det_manip/det_manip.hpp>

template<typename FunctionTypeArg>
class detManip : public triqs::det_manip::det_manip<FunctionTypeArg> {

  typedef typename boost::unwrap_reference<FunctionTypeArg>::type  FunctionType;
  typedef typename FunctionType::result_type                       value_type;
  typedef typename FunctionType::argument_type                     xy_type;
  typedef triqs::arrays::matrix<value_type>                        matrix_type;

  public:

  detManip(FunctionTypeArg F,size_t init_size): triqs::det_manip::det_manip<FunctionTypeArg>(F,init_size) {}

  //----------------------- ITERATORS ---------------------------------
  /**
     Iterator on the C operators in the time order
  */
  class C_iterator : public std::iterator<std::bidirectional_iterator_tag, xy_type> {
    // public:
    int index; detManip * D;
  public:
    C_iterator():index(0),D(NULL) {}
    C_iterator(detManip * D_, int n=1):D(D_) { assert(D); index = std::min(size_t(n),D->size()+1) ; assert(index<=D->size()+1); }
    C_iterator(const C_iterator& mit) : index(mit.index),D(mit.D) {}
    C_iterator&  operator=(const C_iterator& rhs) {assert (rhs.D); D = rhs.D; index= rhs.index; return *this;}
    C_iterator& operator++() {++index;return *this;}
    C_iterator& operator--() {--index;return *this;}
    bool operator==(const C_iterator& rhs) {assert (D==rhs.D); return index==rhs.index;}
    bool operator!=(const C_iterator& rhs) {assert (D==rhs.D); return index!=rhs.index;}
    xy_type operator*() {assert(D);  return (D->y_values[D->col_num[index-1]]);}
    xy_type operator->() {assert(D); return D->y_values[D->col_num[index-1]];}
  };
  
  inline C_iterator C_begin() { return C_iterator(this);}
  inline C_iterator C_end()   { return C_iterator(this,this->N+1);} 

  /**
     Reverse Iterator on the C operators in the time order
  */
  class C_reverse_iterator : public std::iterator<std::bidirectional_iterator_tag, xy_type> {
    // public:
    int index; detManip * D;
  public:
    C_reverse_iterator():index(0),D(NULL) {}
    C_reverse_iterator(detManip * D_, int n=1):D(D_) { assert(D);
      index = std::min(size_t(n),D->size()+1) ; assert(index<=D->size()+1);
      index = D->size()+1 - index;
    }
    C_reverse_iterator(const C_reverse_iterator& mit) : index(mit.index),D(mit.D) {}
    C_reverse_iterator&  operator=(const C_reverse_iterator& rhs) {assert (rhs.D); D = rhs.D; index= rhs.index; return *this;}
    C_reverse_iterator& operator++() {--index;return *this;}
    C_reverse_iterator& operator--() {++index;return *this;}
    bool operator==(const C_reverse_iterator& rhs) {assert (D==rhs.D); return index==rhs.index;}
    bool operator!=(const C_reverse_iterator& rhs) {assert (D==rhs.D); return index!=rhs.index;}
    xy_type operator*() {assert(D);  return (D->y_values[D->col_num[index-1]]);}
    xy_type operator->() {assert(D); return D->y_values[D->col_num[index-1]];}
  };
  
  inline C_reverse_iterator C_rbegin() { return C_reverse_iterator(this);}
  inline C_reverse_iterator C_rend()   { return C_reverse_iterator(this,this->N+1);} 

  /**
     Iterator on the Cdagger operators in the time order

     \code
     // if DET is a detManip object : 
     for (Cdagger_iterator p= DET->Cdagger_begin(); (p != DET->Cdagger_end()) ; ++p) 
       p is casted into a xy_type
     \endcode

  */
  class Cdagger_iterator : public std::iterator<std::bidirectional_iterator_tag, xy_type> {
    //public:
  int index; detManip * D;
  public:
    Cdagger_iterator():index(0),D(NULL) {}
    Cdagger_iterator(detManip * D_, int n=1):D(D_) { assert(D); index = std::min(size_t(n),D->size()+1) ; assert(index<=D->size()+1); }
    Cdagger_iterator(const Cdagger_iterator& mit) : index(mit.index),D(mit.D) {}
    Cdagger_iterator&  operator=(const Cdagger_iterator& rhs) {assert (rhs.D); D = rhs.D; index= rhs.index; return *this;}
    Cdagger_iterator& operator++() {++index;return *this;}
    Cdagger_iterator& operator--() {--index;return *this;}
    bool operator==(const Cdagger_iterator& rhs) {assert (D==rhs.D); return index==rhs.index;}
    bool operator!=(const Cdagger_iterator& rhs) {assert (D==rhs.D); return index!=rhs.index;}
    xy_type operator*() {assert(D);  return (D->x_values[D->row_num[index-1]]);}
    xy_type operator->() {assert(D); return D->x_values[D->row_num[index-1]];}
  };
  
  inline Cdagger_iterator Cdagger_begin()  { return Cdagger_iterator(this);}
  inline Cdagger_iterator Cdagger_end()    { return Cdagger_iterator(this,this->N+1);} 

  /**
     Reverse Iterator on the Cdagger operators

     \code
     // if DET is a detManip object : 
     for (Cdagger_iterator p= DET->Cdagger_begin(); (p != DET->Cdagger_end()) ; ++p) 
       p is casted into a xy_type
     \endcode

  */
  class Cdagger_reverse_iterator : public std::iterator<std::bidirectional_iterator_tag, xy_type> {
    //public:
  int index; detManip * D;
  public:
    Cdagger_reverse_iterator():index(0),D(NULL) {}
    Cdagger_reverse_iterator(detManip * D_, int n=1):D(D_) { assert(D); 
      index = std::min(size_t(n),D->size()+1) ; assert(index<=D->size()+1);
      index = D->size()+1 - index;
    }
    Cdagger_reverse_iterator(const Cdagger_reverse_iterator& mit) : index(mit.index),D(mit.D) {}
    Cdagger_reverse_iterator&  operator=(const Cdagger_reverse_iterator& rhs) {assert (rhs.D); D = rhs.D; index= rhs.index; return *this;}
    Cdagger_reverse_iterator& operator++() {--index;return *this;}
    Cdagger_reverse_iterator& operator--() {++index;return *this;}
    bool operator==(const Cdagger_reverse_iterator& rhs) {assert (D==rhs.D); return index==rhs.index;}
    bool operator!=(const Cdagger_reverse_iterator& rhs) {assert (D==rhs.D); return index!=rhs.index;}
    xy_type operator*() {assert(D);  return (D->x_values[D->row_num[index-1]]);}
    xy_type operator->() {assert(D); return D->x_values[D->row_num[index-1]];}
  };
  
  inline Cdagger_reverse_iterator Cdagger_rbegin()  { return Cdagger_reverse_iterator(this);}
  inline Cdagger_reverse_iterator Cdagger_rend()    { return Cdagger_reverse_iterator(this,this->N+1);} 

  /**
     Iterator returning the objects C[j], C_dagger[i], M(j,i) for all i,j.
     Usage : 
       for ( C_Cdagger_M_iterator p(DET); !p.atEnd(); ++p) { 
        p.C(), p.Cdagger(), p.M() are C[j], C_dagger[i], M(j,i) resp. }
     NB : The i,j are returning in an ARBITRARY ORDER !
  */
  class C_Cdagger_M_iterator { 
    //public: struct value { xy_type C, Cdagger; VALFNT M;};
  protected:
    int i,j;
    const int N;
    const std::vector<xy_type> & c, & cdagger;
    const matrix_type & Minv;
  public:
    C_Cdagger_M_iterator(const detManip & D, int n=1):
      i(1), j(1), N(D.size()), c(D.y_values), cdagger(D.x_values), Minv(D.mat_inv) {}

    C_Cdagger_M_iterator(const Cdagger_iterator& mit);//not impl

    C_Cdagger_M_iterator&  operator=(const Cdagger_iterator& rhs) {
      assert (rhs.N==N); assert(&rhs.c = &c); assert(&rhs.cdagger = &cdagger);//assert(&rhs.m = &m);  
      i = rhs.i; j= rhs.j;return *this;}

    C_Cdagger_M_iterator& operator++() { ++j; if (j>N) { j=1; ++i; };  return *this;}
    bool atEnd() const {return (i>N);}
    xy_type C() const {return c[j-1];}
    xy_type Cdagger() const {return cdagger[i-1];}
    value_type M() const { return Minv(j-1,i-1); }
    void init() { i=1; j=1;}
  };

  //-----------------------------------------------------------------

  friend class C_Cdagger_M_iterator;
  
  //-----------------------------------------------------------------
  /**
     Returns the i-th Cdagger operator in the order.
     i is between 1 and this->size()
  */
  inline Cdagger_iterator select_Cdagger(int i)  { 
    assert (i>0); assert (i<=this->N);
    Cdagger_iterator Op(this);
    for (int j=1; (j<i); ++Op, ++j) { }
    return Op; 
  }
  
  /**
     Returns the i-th Cdagger operator in the reverse order.
     i is between 1 and this->size()
  */
  inline Cdagger_reverse_iterator select_reverse_Cdagger(int i)  { 
    assert (i>0); assert (i<=this->N);
    Cdagger_reverse_iterator Op(this);
    for (int j=1; (j<i); ++Op, ++j) { }
    return Op; 
  }

  /**
     Returns the i-th C operator in the order.
     i is between 1 and this->size()
  */
  inline  C_iterator select_C(int i)  { 
    assert (i>0); assert (i<=this->N);
    C_iterator Op(this);
    for (int j=1; (j<i); ++Op, ++j) { }
    return Op; 
  }
    
  /**
     Returns the i-th C operator in the reverse order.
     i is between 1 and this->size()
  */
  inline  C_reverse_iterator select_reverse_C(int i)  { 
    assert (i>0); assert (i<=this->N);
    C_reverse_iterator Op(this);
    for (int j=1; (j<i); ++Op, ++j) { }
    return Op; 
  }

  //----------------------- RECOMPUTE ---------------------------------

  void recomputeFrom(std::vector<xy_type> x_array, std::vector<xy_type> y_array) { 

    this->last_try=0; this->sign =1;
    this->n_opts=0; this->n_opts_max_before_check = 100;

    this->N = x_array.size();
    if (this->N==0) {
      this->det = 1;
      this->reserve(30);
      this->row_num.clear(); this->col_num.clear();
      this->x_values.clear(); this->y_values.clear(); 
      this->mat_inv() = 0;
      return;
    } else {
      this->reserve(2*this->N);
      this->row_num.clear(); this->col_num.clear();
      this->x_values.clear(); this->y_values.clear(); 
      this->mat_inv() = 0;
    }

    for (size_t i=0; i<this->N; ++i) { 
      this->row_num.push_back(i);
      this->col_num.push_back(i); 
      this->x_values.push_back(x_array[i]);
      this->y_values.push_back(y_array[i]);
    }

    for (size_t i=0; i<this->N; ++i) { 
      for (size_t j=0; j<this->N; ++j) {
        this->mat_inv(i,j) = this->f(this->x_values[i], this->y_values[j]);
      }
    }

    triqs::arrays::range R(0,this->N);
    this->det = triqs::arrays::determinant(this->mat_inv(R,R));
    this->mat_inv(R,R) = inverse(this->mat_inv(R,R));

  }

};

#endif
