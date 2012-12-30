
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

#ifndef DETMANIP_H
#define DETMANIP_H

#include <vector>
#include <iterator>
#include <blitz/array.h>
#include <triqs/deprecated/blitzext/blitz_op.hpp>
#include <triqs/utility/blas_headers.hpp>

/**
   Standard manipulations on determinants.

   This class contains a matrix Minv and its inverse M.
   Minv_ij = Delta(tau_i,tauprime_j)
   where Delta is a callable object, tau, tauprime are some object of TAUTYPE_PTR
   which should be a kind of pointer/iterator to a struct containing a time tau.
   
   It can handle the basic operations :
     - add a row and a col at position (i0,j0) resp. (try_insert)
     - remove a row and a col at position (i0,j0) resp. (try_remove)
     both operations returns the value of the det Minv, should the operation 
     be completed.
     - accept_move : complete the change initiated in try_
     
     in all moves, the matrix Minv, and the inverse M are 
     updated with fast BLAS operations.
     
     M and Minv give accesses the elements of the matrices.

     types : 
      - DELTATYPE : The matrix is Minv_ij = Delta(tau_i,tauprime_j)
      - VALTYPE : COMPLEX or double
      - TAUTYPE : must have a default constructor, operator =

   IMPLEMENTATION : 
     - The storage is done in a compact way, using different permutations
     for row and columns, but that should be transparent for the user
       
 */
template<class DELTATYPE,class VALTYPE,class TAUTYPE_PTR>
class detManip
{
public:

  /// See roll_matrix for usage.
  enum RollDirection {None,Up, Down,Left,Right};

  /** Constructor.
     - Delta is the function : Minv_ij = Delta(tau_i, tau_j')
     - Nmax is the maximum size of the Matrix.
       if a bigger size is needed, the matrix and all work arrays
       will resize themselves to 2*Nmax automatically, but this
       will induces some performance penalty if it happens too often.
  */
  detManip(const DELTATYPE & Delta,int Nmax_,VALTYPE Eta_);
  
  /**
     Copy constructor. Makes a deep copy of the data.
  */
  detManip(const detManip<DELTATYPE,VALTYPE,TAUTYPE_PTR> & X);

private:
  // forbid operator =
  detManip<DELTATYPE,VALTYPE,TAUTYPE_PTR> & operator=(const detManip<DELTATYPE,VALTYPE,TAUTYPE_PTR> & X) {
    assert(0);return *this;}
public:

  /**
     Insertion of a row and a column.

     Consider the insertion of a colum Delta(tau_i - tau_j0') and a row Delta(tau_i0 - tau_j')
     at colum j0 and row i0. i.e. the new colum will be col j0, row i0.
     1<= i0,j0 <= N+1, where N is the current size.
     To be precise, the col. j0 will become j0+1, etc, the new col.
     coming at position j0.
     Inserting at N+1 will simply at the new col at the end therefore.
     Returns the ratio of detMinv_new/ det Minv ( should the operation be completed.)
     This routine does NOT make any modification. It has to be completed with accept_move().
  */
  VALTYPE try_insert(int i0, int j0, TAUTYPE_PTR tau, TAUTYPE_PTR taup);
  
  /**
     Removal of a row and a column.

     Remove the colj0 and row i0 from the matrix.
     Returns the ratio of detMinv_new/ det Minv ( should the operation be completed.)
     This routine does NOT make any modification. It has to be completed with accept_move().
  */
  VALTYPE try_remove(int i0, int j0); 

  /**
     Change a column.

     Change the column j0 and the corresponding Tau
     Returns the ratio of detMinv_new/ det Minv ( should the operation be completed.)
     This routine does NOT make any modification. It has to be completed with accept_move().
  */
  VALTYPE try_change_col(int j0, TAUTYPE_PTR Taup);

  /**
     Change a row.

     Change the row i0 and the corresponding Tau
     Returns the ratio of detMinv_new/ det Minv ( should the operation be completed.)
     This routine does NOT make any modification. It has to be completed with accept_move().
  */
  VALTYPE try_change_row(int i0, TAUTYPE_PTR Tau);

  /**
     Finish the move of the last try_xxx called.
     If not try_xxx has been done, or the last operation was accept_move, 
     this is a FATAL error to the code.
   */
  void accept_move();
  
  /**
     "Cyclic Rolling" of the determinant.

     Right : Move the Nth col to the first col cyclically.
     Left : Move the first col to the Nth, cyclically.
     Up : Move the first row to the Nth, cyclically.
     Down : Move the Nth row to the first row cyclically.
     Return -1 is the roll changes the sign of the det, 1 otherwise
  */
  inline int roll_matrix(RollDirection roll);

   /**
     Move the last inserted/removed column by shift steps
     in the direction given by RollDirection
     Return -1 if the move changes the sign of the det, 1 otherwise
   */
  inline int move_matrix(RollDirection roll, int shift);

  /**
     Give the det Minv of the current state of the matrix.
   */
  inline VALTYPE determinant() const {return sign*det;}

  /**
    Gives the det divided by the max between det and Eta
  */
  inline VALTYPE detovermax() const {return ((abs(det) < Eta) && (abs(det) > 1e-4) ? abs(det)/Eta : 1);}
  
  /**
     Returns M(i,j)
   */
  inline VALTYPE M(int i,int j) const {return _M(col_num(i),row_num(j));}
  // warning : need to invert the 2 permutations.
  
  /**
     Returns Minv(i,j).
   */
  inline VALTYPE Minv(int i,int j) const {return _Minv(row_num(i),col_num(j));}

  /**
     Full "update"
     The matrix and the inverse are recomputed from the tau, tauprime
  */
  void regenerate() const;

  /**
     Recompute the Matrix, from a new set of tau,tauprime.
  */
  void recomputeFrom(std::vector<TAUTYPE_PTR> Tau, std::vector<TAUTYPE_PTR> Taup);

  /**
     Put to size 0
  */
  void reinit();

  // DEBUG
  blitz::Array<VALTYPE,2> M_test();
  blitz::Array<VALTYPE,2> Minv_test();

  // quick access to the matrix, with no reordering.
  // OBSOLETE : use the C_Cdagger_M iterator !
  
  //const Array<VALTYPE,2> & get_M() const {return _M;}
  //const vector< TAUTYPE_PTR> & get_Cdags() const { return tau;}
  //const vector<TAUTYPE_PTR> & get_Cs() const { return tauprime;}
  

  /******************************************************

       Iterator classes
   
  *******************************************************/

  /**
     Iterator on the C operators in the time order
  */
  class C_iterator : public std::iterator<std::bidirectional_iterator_tag, TAUTYPE_PTR> {
    // public:
    int index; detManip * D;
  public:
    C_iterator():index(0),D(NULL) {}
    C_iterator(detManip * D_, int n=1):D(D_) { assert(D); index = std::min(n,D->NumberOfC()+1) ; assert(index<=D->NumberOfC()+1); }
    C_iterator(const C_iterator& mit) : index(mit.index),D(mit.D) {}
    C_iterator&  operator=(const C_iterator& rhs) {assert (rhs.D); D = rhs.D; index= rhs.index; return *this;}
    C_iterator& operator++() {++index;return *this;}
    C_iterator& operator--() {--index;return *this;}
    bool operator==(const C_iterator& rhs) {assert (D==rhs.D); return index==rhs.index;}
    bool operator!=(const C_iterator& rhs) {assert (D==rhs.D); return index!=rhs.index;}
    TAUTYPE_PTR operator*() {assert(D);  return (D->tauprime[D->col_num(index)]);}
    TAUTYPE_PTR operator->() {assert(D); return D->tauprime[D->col_num(index)];}
  };
  
  inline C_iterator C_begin() { return C_iterator(this);}
  inline C_iterator C_end()   { return C_iterator(this,N+1);} 

  /**
     Reverse Iterator on the C operators in the time order
  */
  class C_reverse_iterator : public std::iterator<std::bidirectional_iterator_tag, TAUTYPE_PTR> {
    // public:
    int index; detManip * D;
  public:
    C_reverse_iterator():index(0),D(NULL) {}
    C_reverse_iterator(detManip * D_, int n=1):D(D_) { assert(D);
      index = std::min(n,D->NumberOfC()+1) ; assert(index<=D->NumberOfC()+1);
      index = D->NumberOfC()+1 - index;
    }
    C_reverse_iterator(const C_reverse_iterator& mit) : index(mit.index),D(mit.D) {}
    C_reverse_iterator&  operator=(const C_reverse_iterator& rhs) {assert (rhs.D); D = rhs.D; index= rhs.index; return *this;}
    C_reverse_iterator& operator++() {--index;return *this;}
    C_reverse_iterator& operator--() {++index;return *this;}
    bool operator==(const C_reverse_iterator& rhs) {assert (D==rhs.D); return index==rhs.index;}
    bool operator!=(const C_reverse_iterator& rhs) {assert (D==rhs.D); return index!=rhs.index;}
    TAUTYPE_PTR operator*() {assert(D);  return (D->tauprime[D->col_num(index)]);}
    TAUTYPE_PTR operator->() {assert(D); return D->tauprime[D->col_num(index)];}
  };
  
  inline C_reverse_iterator C_rbegin() { return C_reverse_iterator(this);}
  inline C_reverse_iterator C_rend()   { return C_reverse_iterator(this,N+1);} 

  /**
     Iterator on the Cdagger operators in the time order

     \code
     // if DET is a detManip object : 
     for (Cdagger_iterator p= DET->Cdagger_begin(); (p != DET->Cdagger_end()) ; ++p) 
       p is casted into a TAUTYPE_PTR
     \endcode

  */
  class Cdagger_iterator : public std::iterator<std::bidirectional_iterator_tag, TAUTYPE_PTR> {
    //public:
  int index; detManip * D;
  public:
    Cdagger_iterator():index(0),D(NULL) {}
    Cdagger_iterator(detManip * D_, int n=1):D(D_) { assert(D); index = std::min(n,D->NumberOfC()+1) ; assert(index<=D->NumberOfC()+1); }
    Cdagger_iterator(const Cdagger_iterator& mit) : index(mit.index),D(mit.D) {}
    Cdagger_iterator&  operator=(const Cdagger_iterator& rhs) {assert (rhs.D); D = rhs.D; index= rhs.index; return *this;}
    Cdagger_iterator& operator++() {++index;return *this;}
    Cdagger_iterator& operator--() {--index;return *this;}
    bool operator==(const Cdagger_iterator& rhs) {assert (D==rhs.D); return index==rhs.index;}
    bool operator!=(const Cdagger_iterator& rhs) {assert (D==rhs.D); return index!=rhs.index;}
    TAUTYPE_PTR operator*() {assert(D);  return (D->tau[D->row_num(index)]);}
    TAUTYPE_PTR operator->() {assert(D); return D->tau[D->row_num(index)];}
  };
  
  inline Cdagger_iterator Cdagger_begin()  { return Cdagger_iterator(this);}
  inline Cdagger_iterator Cdagger_end()    { return Cdagger_iterator(this,N+1);} 

  /**
     Reverse Iterator on the Cdagger operators

     \code
     // if DET is a detManip object : 
     for (Cdagger_iterator p= DET->Cdagger_begin(); (p != DET->Cdagger_end()) ; ++p) 
       p is casted into a TAUTYPE_PTR
     \endcode

  */
  class Cdagger_reverse_iterator : public std::iterator<std::bidirectional_iterator_tag, TAUTYPE_PTR> {
    //public:
  int index; detManip * D;
  public:
    Cdagger_reverse_iterator():index(0),D(NULL) {}
    Cdagger_reverse_iterator(detManip * D_, int n=1):D(D_) { assert(D); 
      index = std::min(n,D->NumberOfC()+1) ; assert(index<=D->NumberOfC()+1);
      index = D->NumberOfC()+1 - index;
    }
    Cdagger_reverse_iterator(const Cdagger_reverse_iterator& mit) : index(mit.index),D(mit.D) {}
    Cdagger_reverse_iterator&  operator=(const Cdagger_reverse_iterator& rhs) {assert (rhs.D); D = rhs.D; index= rhs.index; return *this;}
    Cdagger_reverse_iterator& operator++() {--index;return *this;}
    Cdagger_reverse_iterator& operator--() {++index;return *this;}
    bool operator==(const Cdagger_reverse_iterator& rhs) {assert (D==rhs.D); return index==rhs.index;}
    bool operator!=(const Cdagger_reverse_iterator& rhs) {assert (D==rhs.D); return index!=rhs.index;}
    TAUTYPE_PTR operator*() {assert(D);  return (D->tau[D->row_num(index)]);}
    TAUTYPE_PTR operator->() {assert(D); return D->tau[D->row_num(index)];}
  };
  
  inline Cdagger_reverse_iterator Cdagger_rbegin()  { return Cdagger_reverse_iterator(this);}
  inline Cdagger_reverse_iterator Cdagger_rend()    { return Cdagger_reverse_iterator(this,N+1);} 

  /**
     Iterator returning the objects C[j], C_dagger[i], M(j,i) for all i,j.
     Usage : 
       for ( C_Cdagger_M_iterator p(DET); !p.atEnd(); ++p) { 
        p.C(), p.Cdagger(), p.M() are C[j], C_dagger[i], M(j,i) resp.
     NB : The i,j are returning in an ARBITRARY ORDER !
  */
  class C_Cdagger_M_iterator { 
    //public: struct value { TAUTYPE_PTR C, Cdagger; VALFNT M;};
  protected:
    int i,j, index; 
    const int N, N2,Nmax;
    const std::vector<TAUTYPE_PTR> & c, & cdagger;
    const blitz::Array<VALTYPE,2> & m;
    const VALTYPE * restrict M_ptr;
  public:
    C_Cdagger_M_iterator(const detManip & D, int n=1):
      i(1),j(1),index(0),N(D.NumberOfC()),N2(N*N-1),Nmax(D.Nmax), c(D.tauprime), cdagger(D.tau),m(D._M),M_ptr(D._M.dataFirst()) {}

    C_Cdagger_M_iterator(const Cdagger_iterator& mit);//not impl

    C_Cdagger_M_iterator&  operator=(const Cdagger_iterator& rhs) {
      assert (rhs.N==N); assert(&rhs.c = &c); assert(&rhs.cdagger = &cdagger);assert(&rhs.m = &m);  
      i = rhs.i; j= rhs.j;return *this;}

    C_Cdagger_M_iterator& operator++() { ++index; ++j; if (j>N) { j=1; ++i; index += Nmax-N; };  return *this;}
    //C_Cdagger_M_iterator& 
    //    void operator++() { ++j; if (j>N) { j=1; ++i; };} //  return *this;}
    bool atEnd() const {return (i>N);}
    TAUTYPE_PTR C() const {return c[j];}
    TAUTYPE_PTR Cdagger() const {return cdagger[i];}
    //VALTYPE M() const { return m(j,i);}
    VALTYPE M() const { return M_ptr[index];}
    void init() { i=1; j=1;}
  };

  //-----------------------------------------------------------------

  friend class C_Cdagger_M_iterator;
  
  /** Size of the matrix */ // change this name to size
  inline int NumberOfC() const { return N;}

  //-----------------------------------------------------------------
  /**
     Returns the i-th Cdagger operator in the order.
     i is between 1 and this->NumberOfC()
  */
  inline Cdagger_iterator select_Cdagger(int i)  { 
    assert (i>0); assert (i<=N);
    Cdagger_iterator Op(this);
    for (int j=1; (j<i); ++Op, ++j) { }
    return Op; 
  }
  
  /**
     Returns the i-th Cdagger operator in the reverse order.
     i is between 1 and this->NumberOfC()
  */
  inline Cdagger_reverse_iterator select_reverse_Cdagger(int i)  { 
    assert (i>0); assert (i<=N);
    Cdagger_reverse_iterator Op(this);
    for (int j=1; (j<i); ++Op, ++j) { }
    return Op; 
  }

  /**
     Returns the i-th C operator in the order.
     i is between 1 and this->NumberOfC()
  */
  inline  C_iterator select_C(int i)  { 
    assert (i>0); assert (i<=N);
    C_iterator Op(this);
    for (int j=1; (j<i); ++Op, ++j) { }
    return Op; 
  }
    
  /**
     Returns the i-th C operator in the reverse order.
     i is between 1 and this->NumberOfC()
  */
  inline  C_reverse_iterator select_reverse_C(int i)  { 
    assert (i>0); assert (i<=N);
    C_reverse_iterator Op(this);
    for (int j=1; (j<i); ++Op, ++j) { }
    return Op; 
  }

   //-----------------------------------------------------------------
  const VALTYPE Eta;
  const DELTATYPE & delta; 
  int Nmax_current() const {return Nmax;}

  protected:
  VALTYPE det,newdet,ksi;
  TAUTYPE_PTR tmp_tau, tmp_taup;
  int Nmax,N,sign,newsign;
  blitz::Array<VALTYPE,2> _M,_Minv;
  blitz::Array<VALTYPE,1> MB,MC;
  blitz::Array<int,1> row_num,col_num;
  std::vector<TAUTYPE_PTR > tau,tauprime;
  int last_try,i00,j00,ireal,jreal;
  VALTYPE * restrict _Mfirst,* restrict _Minvfirst;
  VALTYPE * restrict MBfirst,* restrict MCfirst;
  
  // pointer computations
  inline VALTYPE * MptrC(int j) {return (_Mfirst + (j-1)*Nmax);}
  inline VALTYPE * MptrL(int i) {return (_Mfirst + (i-1));}
  inline VALTYPE * MinvptrC(int j) {return (_Minvfirst + (j-1)*Nmax);}
  inline VALTYPE * MinvptrL(int i) {return (_Minvfirst + (i-1));}

  inline VALTYPE max_a_b(VALTYPE a, VALTYPE b) const {
    return ((abs(a) < b) && (abs(a) > 1e-4) ? (a<0 ? -b : b) : a);
  }

  // resize the matrix. If 0, it multiplies the size by 2 (and does not call regenerate).
  void resize(int the_size=0); 
  template<class T> inline void for_resize(blitz::Array<T,1> & a){a.resizeAndPreserve(Nmax);assert(a.isStorageContiguous());}
  template<class T> inline void for_resize(blitz::Array<T,2> & a){a.resizeAndPreserve(Nmax,Nmax);assert(a.isStorageContiguous());}
 
};

 
//------------------------------------------------------------------------------------

template<class DELTATYPE,class VALTYPE,class TAUTYPE_PTR>
detManip<DELTATYPE,VALTYPE,TAUTYPE_PTR>::detManip(const DELTATYPE & Delta,int Nmax_,VALTYPE Eta_):
  Eta(Eta_),
  delta(Delta),
  Nmax(Nmax_), N(0),
  _M(Nmax,Nmax,blitz::fortranArray),
  _Minv(Nmax,Nmax,blitz::fortranArray),
  MB(Nmax,blitz::fortranArray),
  MC(Nmax,blitz::fortranArray),
  row_num(Nmax,blitz::fortranArray),
  col_num(Nmax,blitz::fortranArray)
{
  tau.resize(Nmax+1);
  tauprime.resize(Nmax+1);
  MB = 0;MC=0;
  _M=0;_Minv=0;
  last_try=0;
  det = 1;
  sign =1;
  _Mfirst = _M.dataFirst();
  _Minvfirst = _Minv.dataFirst();
  MBfirst = MB.dataFirst();
  MCfirst = MC.dataFirst();
}

//------------------------------------------------------------------------------------

template<class DELTATYPE,class VALTYPE,class TAUTYPE_PTR>
detManip<DELTATYPE,VALTYPE,TAUTYPE_PTR>::detManip(const detManip<DELTATYPE,VALTYPE,TAUTYPE_PTR> & X):
  Eta(X.Eta),
  delta(X.delta),
  Nmax(X.Nmax), N(X.N),
  _M(X._M.copy()),
  _Minv(X._Minv.copy()),
  MB(X.MB.copy()),
  MC(X.MC.copy()),
  row_num(X.row_num.copy()),
  col_num(X.col_num.copy()),
  tau(X.tau),tauprime(X.tauprime),
  last_try(X.last_try),i00(X.i00),j00(X.j00),ireal(X.ireal),jreal(X.jreal)
{
  det=X.det; newdet=X.newdet; ksi=X.ksi;
  sign = X.sign; newsign = X.newsign;
  _Mfirst = _M.dataFirst();
  _Minvfirst = _Minv.dataFirst();
  MBfirst = MB.dataFirst();
  MCfirst = MC.dataFirst();
}

//------------------------------------------------------------------------------------

template<class DELTATYPE,class VALTYPE,class TAUTYPE_PTR>
void detManip<DELTATYPE,VALTYPE,TAUTYPE_PTR>::recomputeFrom(std::vector<TAUTYPE_PTR> Tau, std::vector<TAUTYPE_PTR> Taup)
 {
   assert (Tau.size() == Taup.size());
   // recompute the matrix 
   N =Tau.size(); 
   if  (N>Nmax) resize(2*N); // put some margin..
  
   for (int i=1; i<=N; ++i) { tau[i] = Tau[i-1];tauprime[i] = Taup[i-1];}
   for (int i=1; i<=N; ++i) 
     for (int j=1; j<=N; ++j)
       _Minv(i,j) =delta(tau[i],tauprime[j]);
   _M = _Minv;
   //    cout<<"Minv  "<<_Minv(Range(1,N),Range(1,N))<<" "<<N<<"  "<<Nmax<<endl;
   if (N==0) 
     det =1; 
   else {
     try {det = Blitz_OP::Inverse_Matrix_with_lapack_and_Return_Determinant(_Mfirst,N,Nmax, " In detManip.recomputeFrom ");}
     catch (std::string) {det=0;} // this can happen !
   }
   for (int i =1; i<=N; i++) { row_num(i)= i;col_num(i)= i;}
   sign = 1;// sign of the permutation due to storage....
   last_try = 0;
}
//------------------------------------------------------------------------------------

template<class DELTATYPE,class VALTYPE,class TAUTYPE_PTR>
void detManip<DELTATYPE,VALTYPE,TAUTYPE_PTR>::reinit()
 {
   N = 0;
   sign = 1;det =1;
   last_try = 0;
}

//------------------------------------------------------------------------------------------

template<class DELTATYPE,class VALTYPE,class TAUTYPE_PTR>
VALTYPE detManip<DELTATYPE,VALTYPE,TAUTYPE_PTR>::try_insert(int i0, int j0, TAUTYPE_PTR Tau, TAUTYPE_PTR Taup)
{
  // check input and store it for accept_move
  i00=i0;j00=j0; last_try = 1;
  tmp_tau = Tau;tmp_taup = Taup;
  if (N==Nmax) resize();
  assert(i0<=N+1);  assert(j0<=N+1);
  assert(i0>0); assert(j0>0);

  // treat empty matrix separately (blas on vectors of size 0 ??)
  if (N==0) {
    _Minv(1,1) = delta(tmp_tau , tmp_taup);//delta(Tau , Taup);
    newdet = _Minv(1,1);
    newsign = 1;
    return max_a_b(newdet,Eta);
  }
  else {
    // I add the row and col and the end. If the move is rejected,
    // no effect since N will not be changed : Minv(i,j) for i,j>N 
    // has no meaning.
    
    VALTYPE * restrict C = MinvptrC(N+1), * restrict L = MinvptrL(N+1);
    
    for (int i= 1; i<=N;i++, C++, L+=Nmax)
      {
	//_Minv(i,N+1) = delta(tau[i] , Taup);
	//_Minv(N+1,i) = delta(Tau , tauprime[i]);
        *C  = delta(tau[i] , tmp_taup);
	*L  = delta(tmp_tau , tauprime[i]);
      }

    // MB(i) = sum_j M(i,j) Minv(j,N+1) using BLAS
    const char CN = 'N';
    const double Done = 1, Dzero = 0;
    const int one=1;
    MYFORTRAN(dgemv)(&CN,N,N,Done, _Mfirst,Nmax,MinvptrC(N+1) ,one,Dzero,MBfirst,one);
 
    MB(N+1) = -1;
     
    // ksi = Delta(tau,tauP) - CMB using BLAS
    ksi = delta(Tau , Taup);
    _Minv(N+1,N+1) = ksi;
    ksi -=  MYFORTRAN(ddot)(N,MinvptrL(N+1) ,Nmax,MBfirst,one);

    // compute the newdet
    newdet = det*ksi;
    newsign = ((i0 + j0)%2==0 ? sign : -sign); // since N-i0 + 1 + N-j0 +1 = i0+j0 [2]
    
    VALTYPE r1 = max_a_b(newdet,Eta);
    VALTYPE r2 = max_a_b(det,Eta);

    return (r1/r2)*(newsign*sign);
  }
} 

//------------------------------------------------------------------------------------------

template<class DELTATYPE,class VALTYPE,class TAUTYPE_PTR>
VALTYPE detManip<DELTATYPE,VALTYPE,TAUTYPE_PTR>::try_remove(int i0, int j0)
{
  assert(i0<=N);  assert(j0<=N); 
  assert(i0>0); assert(j0>0);
  i00=i0;j00=j0;last_try = 2;
  jreal = col_num(j00);
  ireal = row_num(i00);
  // compute the newdet
  // first we resolve the ireal,jreal, with the permutation of the Minv, then we pick up what
  // will become the 'corner' coefficient, if the move is accepted, after the exchange of row and col.
  // See accept_move
  VALTYPE ksi = _M(jreal,ireal);
  newdet = det*ksi;
  newsign = ((i0 + j0)%2==0 ? sign : -sign);

  VALTYPE r1 = max_a_b(newdet,Eta);
  VALTYPE r2 = max_a_b(det,Eta);

  return (r1/r2)*(newsign*sign);
}

//------------------------------------------------------------------------------------------

template<class DELTATYPE,class VALTYPE,class TAUTYPE_PTR>
VALTYPE detManip<DELTATYPE,VALTYPE,TAUTYPE_PTR>::try_change_col(int j0, TAUTYPE_PTR Taup)
{
  assert(j0<=N); assert(j0>0);
  j00=j0;last_try = 3;
  jreal = col_num(j0);
  tmp_taup = Taup;

  // Compute the col B.
  for (int i= 1; i<=N;i++)
    MC(i) = delta(tau[i] , tmp_taup) - _Minv(i,jreal);

  // MB(i) = sum_j M(i,j) *MC(j) using BLAS 
  const char CN = 'N';
  const double Done = 1, Dzero = 0;
  const int one=1;
  MYFORTRAN(dgemv)(&CN,N,N,Done, _Mfirst,Nmax,MCfirst ,one,Dzero,MBfirst,one);

  // compute the newdet
  VALTYPE ksi = (1+MB(jreal));
  newdet = det*ksi;
  newsign = sign;
 
  VALTYPE r1 = max_a_b(newdet,Eta);
  VALTYPE r2 = max_a_b(det,Eta);

  return (r1/r2);
}

//------------------------------------------------------------------------------------------

template<class DELTATYPE,class VALTYPE,class TAUTYPE_PTR>
VALTYPE detManip<DELTATYPE,VALTYPE,TAUTYPE_PTR>::try_change_row(int i0, TAUTYPE_PTR Tau)
{
  assert(i0<=N); assert(i0>0);
  i00=i0;last_try = 4;
  ireal = row_num(i0);
  tmp_tau = Tau;

  // Compute the col B.
  for (int i= 1; i<=N;i++)
    MB(i) = delta(tmp_tau, tauprime[i] ) - _Minv(ireal,i);

  // MC(i) = sum_j M(j,i) *MB(j) using BLAS !! transposed !
  const char CT = 'T';
  const double Done = 1, Dzero = 0;
  const int one=1;
  MYFORTRAN(dgemv)(&CT,N,N,Done, _Mfirst,Nmax,MBfirst ,one,Dzero,MCfirst,one);

  // compute the newdet
  VALTYPE ksi = (1+MC(ireal));
  newdet = det*ksi;
  newsign = sign;

  VALTYPE r1 = max_a_b(newdet,Eta);
  VALTYPE r2 = max_a_b(det,Eta);

  return (r1/r2);
}

//------------------------------------------------------------------------------------------

template<class DELTATYPE,class VALTYPE,class TAUTYPE_PTR>
void detManip<DELTATYPE,VALTYPE,TAUTYPE_PTR>::accept_move()
{
  // for BLAS
  const char CT = 'T';
  const double Done = 1, Dzero = 0;
  const int one=1;

  // according to which was the last try_xxx called, I finished the moves.
  switch(last_try){
  case(1): // finish the insert move

    // store the new tau. They are see through the same permutations as 
    // rows and cols resp.
    tau[N+1]= tmp_tau;  tauprime[N+1]= tmp_taup;

    // special empty case again
    if (N==0)
      {
	N=1;
	row_num(1) = 1;
	col_num(1) = 1;
	_M(1,1) = 1/_Minv(1,1);
      }
    else 
      {
	// MC(i) = sum_j M(j,i) Minv(N+1,j) using BLAS
	MYFORTRAN(dgemv)(&CT,N,N,Done, _Mfirst,Nmax,MinvptrL(N+1) ,Nmax,Dzero,MCfirst,one);
	MC(N+1) = -1;
	
	N++;
	
	// keep the real position of the row/col
	// since we insert a col/row, we have first to push the col at the right
	// and then say that col i00 is stored in N, the last col.
	// same for rows.
	for (int i =N-1; i>=i00; i--) row_num(i+1)= row_num(i);
	row_num(i00) = N;
	for (int i =N-1; i>=j00; i--) col_num(i+1)= col_num(i);
	col_num(j00) = N;
		
	// Minv is ok, we need to complete 
	ksi = 1/ksi;

	// compute the change to the inverse
	// M += ksi MB MC with BLAS
	// first put the 0 
	//_M(Range(1,N),N) = 0;_M(N,Range(1,N)) = 0;
	int i=0; for (VALTYPE * restrict p= MptrC(N),* restrict q=MptrL(N); i<N;i++,p++,q+=Nmax) {*p=0; *q=0;}

	MYFORTRAN(dger)(N,N,ksi,MBfirst,one,MCfirst,one,_Mfirst,Nmax);
      }

    break;
  
  case(2):  // finish the remove move

    // special case
    if (N==1) { N=0; row_num(1) = 0; col_num(1) = 0;  break;}

    // repack the matrix _Minv
    // swap the rows ireal and N, jreal and N in _Minv
    // Remember that for M row/col is interchanged by inversion, transposition.
    // TEST the speed of blas 1 versus Blitz for these operations ?
   
    if (jreal !=N){
      //_Minv(Range(1,N),jreal) = _Minv(Range(1,N),N);
      MYFORTRAN(dcopy)(N, MinvptrC(N),one, MinvptrC(jreal),one);
      tauprime[jreal] = tauprime[N]; 
      MYFORTRAN(dswap)(N, MptrL(jreal) ,Nmax, MptrL(N),Nmax);
    }
    
    if (ireal !=N){
      //_Minv(ireal,Range(1,N)) = _Minv(N,Range(1,N));
      MYFORTRAN(dcopy)(N, MinvptrL(N),Nmax, MinvptrL(ireal),Nmax);
      tau[ireal] = tau[N];
      MYFORTRAN(dswap)(N, MptrC(ireal) ,one, MptrC(N) ,one);
    }

    N--;

    // M <- a - d^-1 b c with BLAS
    ksi = - 1/_M(N+1,N+1);
    MYFORTRAN(dger)(N,N,ksi,MptrC(N+1),one,MptrL(N+1),Nmax,_Mfirst,Nmax);

    // modify the permutations
    // see case 1
    for (int i =i00; i<=N; i++) {row_num(i)= row_num(i+1);}
    for (int i =j00; i<=N; i++) {col_num(i)= col_num(i+1);}
    for (int i =1; i<=N; i++) 
      {
	if (col_num(i)==N+1) col_num(i)=jreal;
	if (row_num(i)==N+1) row_num(i)=ireal;
      }
    row_num(N+1) = 0; col_num(N+1) = 0; 

   break;

  case(3): // finish the change_col move
    // modifying the M_inv
    //_Minv(ALL,jreal) += MC(ALL); // MC here is the B colum : see try_change_col
    { //optimised version
      VALTYPE * restrict PP =MinvptrC(jreal);
      for (int u =0; u<N; ++u) PP[u] += MCfirst[u]; 
    }
    
    tauprime[jreal] = tmp_taup;

    // modifying M : Mij += ksi Bi Mnj 
    // using Shermann Morrison formula.
    // implemented in 2 times : first Bn=0 so that Mnj is not modified ! and then change Mnj
    // Cf notes : simply multiply by -ksi
    ksi = - 1/(1+ MB(jreal));
    MB(jreal) = 0;
    MYFORTRAN(dger)(N,N,ksi,MBfirst,one,MptrL(jreal),Nmax,_Mfirst,Nmax);
    ksi *= -1;
    MYFORTRAN(dscal)(N,ksi,MptrL(jreal),Nmax);
    //_M(jreal,ALL) *= -ksi;

    break;
 
  case(4): // finish the change_row move
    // modifying the M_inv
    //_Minv(ireal,ALL) += MB(ALL); // MB here is the C row 
    { // optimised version
      VALTYPE * restrict PP =MinvptrL(ireal);
      for (int u =0,u2=0; u<N; ++u,u2+=Nmax) PP[u2] += MBfirst[u]; 
    }
   
    tau[ireal] = tmp_tau;

    // modifying M : M ij += ksi Min Cj
    // using Shermann Morrison formula.
    // impl. Cf case 3
    ksi = - 1/(1+ MC(ireal));
    MC(ireal) = 0;
    MYFORTRAN(dger)(N,N,ksi,MptrC(ireal),one,MCfirst,one,_Mfirst,Nmax);
    ksi *= -1;
    MYFORTRAN(dscal)(N,ksi,MptrC(ireal),one);
    //_M(ALL,ireal) *= -ksi;

    break;

  case(0):  // make sure we have a move to finish : this is a fatal error if we don't
    TRIQS_RUNTIME_ERROR<<"misuing detManip";
    
  }
  
  det = newdet;
  sign = newsign;
}

//----------------------------------------------------------------------

template<class DELTATYPE,class VALTYPE,class TAUTYPE_PTR>
inline int detManip<DELTATYPE,VALTYPE,TAUTYPE_PTR>::roll_matrix(RollDirection roll)
{
  int tmp;
  switch (roll) {
  case(None) : 
    return 1;
  case(Down) : 
    tmp = row_num(N);
    for (int i =N-1; i>=1; i--) row_num(i+1)= row_num(i);
    row_num(1) = tmp;
    break;
  case(Up) : 
    tmp = row_num(1);
    for (int i =1; i<=N-1; i++) row_num(i)= row_num(i+1);
    row_num(N) = tmp;
    break;
  case(Right) : 
    tmp = col_num(N);
    for (int i =N-1; i>=1; i--) col_num(i+1)= col_num(i);
    col_num(1) = tmp;
    break;
  case(Left): 
    tmp = col_num(1);
    for (int i =1; i<=N-1; i++) col_num(i)= col_num(i+1);
    col_num(N) = tmp;
    break;
  default:
    assert(0);
  }
  // signature of the cycle of order N : (-1)^(N-1)
  if ((N-1)%2==1) { sign *=-1; return -1;}
  return 1;
}

//----------------------------------------------------------------------

template<class DELTATYPE,class VALTYPE,class TAUTYPE_PTR>
inline int detManip<DELTATYPE,VALTYPE,TAUTYPE_PTR>::move_matrix(RollDirection roll, int shift)
{
  int tmp, num;
  switch (roll) {
  case(None) : 
    return 1;
  case(Down) : 
    num = i00;
    tmp = row_num(num);
    for (int i =num; i<num+shift; i++) row_num(i)= row_num(i+1);
    row_num(num+shift) = tmp;
    break;
  case(Up) : 
    num = i00;
    tmp = row_num(num);
    for (int i =num; i>num-shift; i--) row_num(i)= row_num(i-1);
    row_num(num-shift) = tmp;
    break;
  case(Right) : 
    num = j00;
    tmp = col_num(num);
    for (int i =num; i<num+shift; i++) col_num(i)= col_num(i+1);
    col_num(num+shift) = tmp;
    break;
  case(Left): 
    num = j00;
    tmp = col_num(num);
    for (int i =num; i>num-shift; i--) col_num(i)= col_num(i-1);
    col_num(num-shift) = tmp;
    break;
  default:
    assert(0);
  }
  // signature of the cycle of order N : (-1)^(N-1)
  if (shift%2==1) { sign *=-1; return -1;}
  return 1;
}

//----------------------------------------------------------------------

template<class DELTATYPE,class VALTYPE,class TAUTYPE_PTR>
void detManip<DELTATYPE,VALTYPE,TAUTYPE_PTR>::resize(int the_size)
{
  //std::cout <<" Node : " << myMPI_myid() << std::endl;
  boost::mpi::communicator world;
  std::cout <<" Node : " << world.rank() << std::endl;
  std::cout <<" WARNING : RESIZING THE MATRIX : it will be slow if it happens too often !"<<std::endl;
  // double the size of the arrays.
  Nmax = (the_size !=0  ? the_size : 2*Nmax);
  for_resize(_M);  for_resize(_Minv);
  for_resize(MB);  for_resize(MC);
  for_resize(row_num);  for_resize(col_num);
  tau.resize(Nmax+1);  tauprime.resize(Nmax+1);
  _Mfirst = _M.dataFirst();
  _Minvfirst = _Minv.dataFirst();
  MBfirst = MB.dataFirst();
  MCfirst = MC.dataFirst();
  if (the_size==0) regenerate();
  //cout<<" size done"<<endl;
}

//----------------------------------------------------------------------

template<class DELTATYPE,class VALTYPE,class TAUTYPE_PTR>
void detManip<DELTATYPE,VALTYPE,TAUTYPE_PTR>::regenerate() const
{
  if (N!=0) {
    blitz::Array<VALTYPE,2> res(N,N,blitz::fortranArray);

    for (int i=1; i<=N;i++)
      for (int j=1; j<=N;j++)
        res(i,j) =  delta(tau[i],tauprime[j]);

    blitz::Range R(1,N);

    if ((max(abs(res-_Minv(R,R)))/(max(abs(_Minv(R,R)))+1.e-15))>1.e-10) {
      std::cout <<" Node : " << myMPI_myid() << std::endl;
      std::cout<<"Minv SHOULD BE "<<res;
      std::cout<<"AND is "<<_Minv(R,R)<<std::endl; }
    
    _Minv(R,R) = res;

    /*
    cout<<_Minv(R,R)<<_M(R,R)<<endl;
    cout<<res<<endl;
    */

    Blitz_OP::Inverse_Matrix(res);

    //cout<<"APRES INV"<<res<<endl;

    if ((max(abs(res-_M(R,R)))/(max(abs(_M(R,R)))+1.e-15))>1.e-10) {
      std::cout <<" Node : " << myMPI_myid() << std::endl;
      std::cout<<"M SHOULD BE "<<res;
      std::cout<<"AND is "<<_M(R,R)<<std::endl; }
    
    _M(R,R) = res;
  
  }
}


//----------------------------------------------------------------------

template<class DELTATYPE,class VALTYPE,class TAUTYPE_PTR>
blitz::Array<VALTYPE,2> detManip<DELTATYPE,VALTYPE,TAUTYPE_PTR>::Minv_test()
{
  blitz::Array<VALTYPE,2> res(N,N,blitz::fortranArray);
  blitz::Array<VALTYPE,2> res2(N,N,blitz::fortranArray),res3(N,N,blitz::fortranArray);

  res = 0;
  if (N==0) return res;

  for (int i=1; i<=N;i++)
    for (int j=1; j<=N;j++)
      res2(i,j) =  delta(tau[row_num(i)],tauprime[col_num(j)]);

  for (int i=1; i<=N;i++)
    for (int j=1; j<=N;j++)
      res(i,j) = Minv(i,j);

  if (max(abs(res-res2))>1.e-10) {
    std::cout<<"SHOULD BE "<<res2;
    std::cout<<"AND is  "<<res<<std::endl;
    assert(0);}
  
  set_Scalar(res3,1);
  blitz::Range R(1,N);
  res2 = matmul(_Minv(R,R),_M(R,R));
  if (max(abs(res3-res2))>1.e-10) {
    std::cout<<"Inversion not ok"<<res2<<_Minv(R,R)<<_M(R,R);
    //assert(0);
  }

  return res;
}


template<class DELTATYPE,class VALTYPE,class TAUTYPE_PTR>
blitz::Array<VALTYPE,2> detManip<DELTATYPE,VALTYPE,TAUTYPE_PTR>::M_test()
{
  blitz::Array<VALTYPE,2> res(N,N,blitz::fortranArray);
  for (int i=1; i<=N;i++)
    for (int j=1; j<=N;j++)
      res(i,j) = M(i,j);
  return res;
}

#endif
