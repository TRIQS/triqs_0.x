
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

#ifndef TRIQS_BASE_TAIL_GF
#define TRIQS_BASE_TAIL_GF

#include "Python.h"
#include <boost/python.hpp>
#include <boost/python/slice.hpp>
#ifndef NO_USE_BLITZ
  #include "triqs/deprecated/blitzext/blitz_op.hpp"
  #include "triqs/deprecated/blitzext/PyArray.hpp"
  using namespace Blitz_OP;
#endif
#include <triqs/utility/mathlib.hpp>
#include <triqs/utility/exceptions.hpp>

using namespace std;
using namespace boost;
using namespace blitz;


/* A small class to take care of the tails of GF.
   It stores the coefficients of the large omega expansion in matrix
   from M1/omega^OrderMin to M2/omega^OrderMax.
*/
class TailGF { 
protected :
  // Indices are the left and right indices of the matrix
  const python::list IndicesL,IndicesR;

public:
  const int N1,N2; // Size of the matrix
  const int OrderMinMIN, OrderMaxMAX; // maxi dimension of the expansion. Fixed at construction.
 
protected:
  PyArray<COMPLEX,3> M; // the coefficients of the tail.
  PyArray<int,2> OrderMaxArray; //  OrderMax of the expansion for each elements
  // since one can take a sliced view of the tail, this may depend on the element.
  // it will be renormalized to min(OrderMaxArray) by calling OrderMax(), 
  // which is done before any matrix operation.

public : 
 
  // Construct a new tail, with the OrderMax, OrderMin and coefs = 0
  TailGF(int OrderMinMIN_, int OrderMaxMAX_, python::list IndicesL_,  python::list IndicesR_);

  // Make a slice VIEW of the tail t : t[slL,slR]
  TailGF(const TailGF & t, python::object slL, python::object slR);

protected:
  // for reconstruction
  TailGF(int OrderMinMIN_, int OrderMaxMAX_, python::object Arr, int OrderMax_, 
	 python::list IndicesL_,  python::list IndicesR_ );

public : 

  // Make a VIEW of the same tail
  TailGF(const TailGF & t);

  // -------   Access to data -----------------

 // Starting coefficient 
  inline int OrderMin() const { int omin=OrderMinMIN;const int omax(OrderMax()); while ((omin<omax) && (max(abs(M(omin,ALL,ALL))) < ZERO)) omin++; return omin;}

  // where is the dev stopping ? o(1/omega^OrderMax+1)
  // for non-const, we normalize the OrderMaxArray.
  inline int OrderMax_py() { return OrderMax();}// boost python can not resolve the const overloading
  inline int OrderMax() { int r = min(OrderMaxArray); OrderMaxArray = r; return r;}
  inline int OrderMax() const { int r = min(OrderMaxArray);  return r;}
  
  // Sets the expansion to 0. Reinit the Tail.
  inline void zero(){ M=0; OrderMaxArray = OrderMaxMAX;}

  // I hope you know what you are doing ...(used when fitting higher order coef of tail)
  inline void changeOrderMax(int newOrderMax) { //if (newOrderMax>OrderMax()) TRIQS_RUNTIME_ERROR<<"I can not INCREASE OrderMax. Sorry...";
    int i = OrderMax(); if (newOrderMax>i) {M(Range(i+1,newOrderMax),ALL,ALL) = 0;} OrderMaxArray=newOrderMax;}

  // whether the coefficient order exists and makes sense
  inline bool has_coef(int order) const { 
    return ( (order <= OrderMax()) && (order >= OrderMinMIN));
  }

  // [i] return a VIEW of the 1/omega^i term if the coef exists.
  inline Array<COMPLEX,2> operator[](int order) const { 
    if (!has_coef(order)) 
      TRIQS_RUNTIME_ERROR<<"TailGF:: operator[] :: order "<<order<<" incorrect. OrderMin="<< OrderMin()<<"  OrderMax = "<<OrderMax();
    return Array<COMPLEX,2> (M(order,ALL,ALL));}

  //Evaluate the tail at frequency omega.
  Array<COMPLEX,2> eval(COMPLEX omega) const;

  // True iif the expansion starts at 1/omega^i with i>0
  inline bool is_decreasing_at_infinity() const { return (OrderMin() >=1);}

  // ------------ for Python interface ------------

  // Returns coefficient i as ArrayViewWithIndexConverter
  python::object __getitem__(int i) const;

  // Sets the coefficient i to val
  void __setitem__(int i,const PyArray<COMPLEX,2> & val);

  // Returns all coefficients as a dict : order -> ArrayViewWithIndexConverter
  python::object AllCoefs() const;

  //
  python::object __repr__() const;

  // Like eval, put returns a ArrayViewWithIndexConverter for python usage.
  python::object __call__(COMPLEX omega) const;

  // for python  = property
  inline void copyFrom(const TailGF & t){ *this = t;}

  // -------------  Various operators -----------------

  // Copy t into this.
  TailGF & operator = (const TailGF & t); 

  // Are the 2 tail equal ?
  bool operator==(const TailGF &other) const;

  TailGF & operator += (const TailGF & t);
  TailGF & operator -= (const TailGF & t);
     
  inline TailGF & operator *= (COMPLEX alpha) {M *=alpha; return *this;}
  inline TailGF & operator /= (COMPLEX alpha) {M /=alpha; return *this;}

  TailGF & operator *= (const TailGF & t); 

  TailGF transpose() const; // return a new view with transpose array, as numpy
  TailGF conjugate(bool is_matsubara_expansion) const; // return a new view with conjugate array , as numpy
  // if is_matsubara_expansion, then add a (-1)^n to take into account the fact that the denominator is (i omega_n)
  
  // ------------ Other operations -------------------
  
  // Give a copy of the tail
  TailGF copy();

  // Replace itself by the expansion of the inverse of the function.
  void invert();

  // Sets the tail to be sum_i (L M_i R)/omega^i (matrix product)
  void from_L_T_R(const PyArray<COMPLEX,2> & L, const TailGF & T2, const PyArray<COMPLEX,2> & R);

  // --------------- IO ---------------------

  // Save to text file 
  void save(string file,  bool accumulate) const;

  // load from text file. inverse of save.
  void load(string file);

  // Give the reducing protocol to the object.
  python::object __reduce_to_dict__() const;
  static python::object __factory_from_dict__(const python::object & dic);

  python::tuple __reduce__() const;
protected: 
  // check compatibility 
  inline void check( const TailGF & t);
  
  template<typename TYPE>
  void check(const Array<TYPE,2> & tab) { 
    if (tab.extent(0)!=N1) TRIQS_RUNTIME_ERROR<<"Tail operator : can not add/substract this array : first dimension incorrect";	
    if (tab.extent(1)!=N2) TRIQS_RUNTIME_ERROR<<"Tail operator : can not add/substract this array : second dimension incorrect"; 
    if (! ( (OrderMinMIN <=0) && (OrderMaxMAX>=0)))
      TRIQS_RUNTIME_ERROR<<"Adding a constant impossible for this tail : it does not fit in OrderMinMIN:OrderMaxMAX window"; 
  if (OrderMax()<0) TRIQS_RUNTIME_ERROR<<"Adding a constant has no effect on this expansion since OrderMax<0. Please check what you are doing !"; 
  }
  
  inline void check_square() const {
    if (N1!=N2) TRIQS_RUNTIME_ERROR<<"This operation can only be done for square matrix";
  }

  // Puts M to 0 *outside* of M(range(wmin,wmax))
  void cleanM_outside(int wmin, int wmax);
  

};


#endif

