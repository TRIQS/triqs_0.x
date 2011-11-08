
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

#ifndef BLITZ_OP
#define BLITZ_OP
#define USE_BLITZ_ARRAY
  
#include <vector>
#include <fstream>
#include <blitz/array.h>
#include <triqs/utility/report_stream.hpp>
#include <triqs/deprecated/utility/my_mpi.hpp>

/**
 * \defgroup blitzop Extensions of Blitz
 @{
*/

#define FATAL(s) { std::stringstream fs_for_fatal; fs_for_fatal<<s; throw fs_for_fatal.str();}

namespace Blitz_OP {

///
#define ALL blitz::Range::all()

/**
 Nicer printing of an array
*/
template <class T>
std::string pretty_print(const blitz::Array<T,1> & A)
{std::stringstream fs; for (int u =A.lbound(blitz::firstDim);u<=A.ubound(blitz::firstDim);u++) fs<<"  "<<A(u)<<"  ";
  return(fs.str());}

/**
 Nicer printing of an array
*/
template <class T>
std::string pretty_print(const blitz::Array<T,2> & A)
{ std::stringstream fs;
  for(int i=A.lbound(0);i<=A.ubound(0);i++)
    for(int j=A.lbound(1);j<=A.ubound(1);j++)
      fs<< i<< "  "<<j<<"  "<<A(i,j)<<std::endl;
  return(fs.str());}


/// Is it fortran ordered ?
template <class T, int n>
bool is_fortran_ordering(const blitz::Array<T,n> & A) 
{bool r=true; for (int i=0; i<n;i++) r = (r && (A.ordering(i)==i));return(r);}

/// Is it C ordered ?
template <class T, int n>
bool is_C_ordering(const blitz::Array<T,n> & A) 
{bool r=true; for (int i=0; i<n;i++) r = (r && (A.ordering(i)==n-i-1));return(r);}


inline const char * c_string(const std::string & s){return (s.c_str());}

/** Interface Blitz Fortran (or C)

    If f is a Fortran function : 
        subroutine f(A,n1,n2,...)
	  TYPE, dimension(n1,n2,...):: A
	  integer n1,n2,...
	  ....
    where TYPE can be integer, double, complex<double>
    properly interfaced into C with 
      extern "C" { void f_(TYPE [], int &);}
      
    and if M is a blitz::Array<TYPE,n>, then 
    f_(ref<TYPE,n>(M)) will be ok.
    
    Indeed : at construction, ref will make sure the array is in Fortran ordering (optionally one can skip this, see constructor)
    and contiguous (mandatory test), then the compiler will recast it into a pointer of type TYPE to the data (or a temporary
    copy). After the function has returned, the destructor of ref will copy the data back to the correct place if needed.

    NB : higher dimensionnal arrays are passed the same way as 1d : transformed into a pointer.
    NB : note that on the Fortran side : 
      - the function MUST NOT be in a module : the name convention for function in a module
        is compiler unfortunately dependent
      - automatic array in f90 do not work : one HAS to call f90 function like above, by passing explicitely
      the dimension (if necessary, write a little wrapper).

    For illustration, see the blitz_lapack for interface of some LAPACK routines (the one I used at least once !).
*/
template<class T,int n>
class ref
{
public:
  ref (blitz::Array<T,n> & M, bool check_order_storage=true ):Mcopy(blitz::fortranArray) { for_cons(M,check_order_storage);}
  ref (const blitz::Array<T,n> & M, bool check_order_storage=true): Mcopy(blitz::fortranArray){ for_cons(M,check_order_storage);}
  ~ref() {if (need_copy) (*Mref) = Mcopy; }
  operator T* (){ return (need_copy ?  Mcopy.dataFirst() :  Mref->dataFirst());}
  bool use_copy(){return(need_copy);}
private:
  inline void for_cons(const blitz::Array<T,n> & M, bool check_order){
    need_copy = (!(M.isStorageContiguous()));
    if (check_order) for (int i=0; i<n;i++) need_copy = (need_copy || (M.ordering(i)!=i));
#ifdef DEBUG_REF_WARNING
    if (need_copy) std::cout<<"WARNING : REF : COPY NEEDED. Performance will be degraded"<<std::endl;
#endif
    Mref = (blitz::Array<T,n> *)&M;
    // The copy has the same shape but is ordered like a fortran array
    if (need_copy) {Mcopy.resize(M.shape());Mcopy=M;}
  }
  blitz::Array<T,n> Mcopy;
  blitz::Array<T,n> *Mref;
  bool need_copy;
};


/*q* @defgroup blitzlapack Blitz and LAPACK
@{
*/

/**
   Matrix multiplication using the BLAS dgemm
*/
void  matmul_lapack(const blitz::Array<double,2> & A, const blitz::Array<double,2> & B,blitz::Array<double,2> & RES);

/**
   Matrix multiplication using the BLAS zgemm
*/
void  matmul_lapack(const blitz::Array<std::complex<double>,2> & A, const blitz::Array<std::complex<double>,2> & B, blitz::Array<std::complex<double>,2> & C);

/// Multiply M -> AMB on site if A, B are square and of correct size
template< typename VAL>
void matmul_A_M_B( const blitz::Array<VAL,2> & A, blitz::Array<VAL,2> & M, const blitz::Array<VAL,2> & B ) { 
  assert (A.extent(0) == A.extent(1));   
  assert (B.extent(0) == B.extent(1));   
  assert (A.extent(1) == M.extent(0));   
  assert (B.extent(0) == M.extent(1));
  blitz::Array<VAL,2> tmp(M.extent(0),B.extent(1), blitz::fortranArray);
  matmul_lapack(M, B ,tmp);
  matmul_lapack( A, tmp, M);
}

double Inverse_Matrix_with_lapack_and_Return_Determinant (double * refM, int n, int LDA, std::string message="");

inline double Inverse_Matrix_with_lapack_and_Return_Determinant (blitz::Array<double,2> & M, std::string message="")
{  
  assert (M.rows() ==M.cols());
  int n = M.rows();
  if (n==0) return 1;
  ref<double,2> refM(M,false);//do not check the ordering since inverse and transposition commutes
  double * p = refM;
  return Inverse_Matrix_with_lapack_and_Return_Determinant(p,n,n,message);
}


/**   
      Inversion using lapack
*/
bool Inverse_Matrix_with_lapack (blitz::Array<double,2> & M);

/**   
      Inversion using lapack
*/
bool Inverse_Matrix_with_lapack (blitz::Array<std::complex<double>,2> & M);

/**   
       \f$ A^{-1}_B \f$  using lapack
*/
bool A_inv_B_with_lapack  (blitz::Array<double,2> & A, blitz::Array<double,2> & B );

/**   
       \f$ A^{-1}_B \f$  using lapack
*/
bool A_inv_B_with_lapack  (blitz::Array<std::complex<double>,2> & A, blitz::Array<std::complex<double>,2> & B );


/**  Diagonalisation  using lapack
 */
int Diagonalise_with_lapack(blitz::Array<double,2> & M, blitz::Array<double,1> & ev , bool Compute_Eigenvectors);

/**  Diagonalisation  using lapack
 */
int Diagonalise_with_lapack(blitz::Array<std::complex<double>,2> & A, blitz::Array<double,1> & ev , bool Compute_Eigenvectors);

/**
 */
blitz::Array<std::complex<double>,1> Root_Polynome_w_Hessenberg_Matrix(const blitz::Array<std::complex<double>,1> & poly); 


/**
   
*/
bool svd1_lapack(blitz::Array<std::complex<double>,2> & M, blitz::Array<std::complex<double>,1> & ev,blitz::Array<std::complex<double>,2> & vl, blitz::Array<std::complex<double>,2> & vr);

bool svd_lapack(blitz::Array<std::complex<double>,2> & M, blitz::Array<double,1> & Sigma,blitz::Array<std::complex<double>,2> & U, blitz::Array<std::complex<double>,2> & VT);

  
/* @}
 */
  

#define C_F_default true

//-----------------------------------------------------------------------------------------------
// BINARY SAVE-LOAD

template<class T, int n>
void binary_save(const blitz::Array<T,n> & M,std::string filename)
{
  T z;assert(M.isStorageContiguous());
  std::stringstream f;M.dumpStructureInformation(f);
  std::ofstream out(filename.c_str(),std::ios::binary);
  std::string s = f.str();int i=int(s.length())+1;char c='1';
  out.write( (char*) &i,sizeof(i));
  out.write( (s.c_str()) ,i*sizeof(c));  
  out.write( (char*)M.dataFirst(),M.numElements()*sizeof(z));
}

template<class T, int n>
void binary_load(blitz::Array<T,n> & M,std::string filename)
{
  assert(M.isStorageContiguous());T z;
  std::ifstream out(filename.c_str(),std::ios::binary);
  int i;char c='1';
  out.read( (char*) &i,sizeof(i));
  char *st = new char[i+1];
  out.read( st ,i*sizeof(c)); std::string s(st); 
  std::stringstream f;  M.dumpStructureInformation(f);
  if (f.str() != s) FATAL("Can not load binary : array do not conform. Structure (file vs array)"<<s<<"----"<<f.str());
  out.read( (char*)M.dataFirst(),M.numElements()*sizeof(z));
}

//-----------------------------------------------------------------------------------------------

template<class T>
std::ostream & operator << (std::ostream &S, const std::vector<T> & V) {
  S<<V.size()<<std::endl;
  for (unsigned int i=0; i< V.size(); ++i) S<<V[i]<<std::endl;
  return S;
}

template<class T>
std::istream & operator >> (std::istream & S,  std::vector<T> & V) {
  unsigned int t; S>>t;V.resize(t);
  for (unsigned int i=0; i< t; ++i) S>>V[i];
  return S;
}
    

//-----------------------------------------------------------------------------------------------

/** Translation blitz <-> vector. The blitz array is resized
 */
/*template<class T, int N>
void  vector_to_blitz2(blitz::Array<T,N> & res, vector<T> & v){
  if (N==1) { int n= v.size(); res.resize(n);int sh =  res.lbound(blitz::firstDim);
    res =0; for (int i=0;i<n;i++) res(i+sh)=v[i];
  }
  else {blitz::Array<T,N-1> tmp(fortranArray); 
    for (int i=0;i<v.size;i++) 
      { vector_to_blitz2(tmp,v[i]);
}
*/

/** Translation blitz <-> vector. The blitz array is resized
 */
template<class T>
void  vector_to_blitz(blitz::Array<T,1> & res, const std::vector<T> & v){
  int n= int(v.size());//64bits
  res.resize(n);int sh =  res.lbound(blitz::firstDim);
  for (int i=0;i<n;i++) res(i+sh)=v[i];
}

/** Translation blitz <-> vector. The blitz array is resized
    res(j,i) = v[i][j]
 */
template<class T>
void vector_to_blitz(blitz::Array<T,2> &res, const std::vector<std::vector<T> > & v){
  int n1= v.size();int n2=0,tmp; for (int i=0;i<n1;i++) {tmp =v[i].size(); n2 = std::max(n2,tmp);}
  res.resize(n2,n1); int sh1 =  res.lbound(blitz::firstDim); int sh2 =  res.lbound(blitz::secondDim);
  for (int i=0;i<n1;i++)  for (int j=0;j<int(v[i].size());j++) res(j+sh1,i+sh2)=v[i][j];
}

/** Translation blitz <-> vector. The blitz array is resized
    res(k,j,i) = v[i][j][k]
 */
template <class T>
void vector_to_blitz(blitz::Array<T,3> & res,const std::vector<std::vector<std::vector<T> > > & v){
  int n1= v.size();int n2=0,tmp,tmp3; for (int i=0;i<n1;i++) {tmp =v[i].size(); n2 = std::max(n2,tmp);} 
  int n3=0; for (int i=0;i<n1;i++)   for (int j=0;j<int(v[i].size());j++) {tmp3 =v[i][j].size(); n3 = std::max(n3,tmp3);} 
  res.resize(n3,n2,n1); int sh1 =  res.lbound(blitz::firstDim); int sh2 =  res.lbound(blitz::secondDim);int sh3 =  res.lbound(blitz::thirdDim);
  for (int i=0;i<n1;i++)  for (int j=0;j<int(v[i].size());j++) for (int k=0;k<int(v[i][j].size());k++)  
	    res(k+sh1,j+sh2,i+sh3)=v[i][j][k];
}

/** Translation blitz <-> vector
 */
template <class T>
std::vector<T> blitz_to_vector(const blitz::Array<T,1> & a){
  std::vector<T> v;
  for (int i=a.lbound(blitz::firstDim);i<=a.ubound(blitz::firstDim);i++) v.push_back(a(i));
  return(v);
}

/** Translation blitz <-> vector
 */
template <class T>
std::vector<std::vector<T> > blitz_to_vector(const blitz::Array<T,2> & a){
  std::vector<std::vector<T> >v; blitz::Array<T,1> tmp(a.extent(blitz::secondDim));
  for (int i=a.lbound(blitz::firstDim);i<=a.ubound(blitz::firstDim);i++) 
    {tmp = a(i,ALL); v.push_back(blitz_to_vector(tmp));}
  return(v);
}

/** Translation blitz <-> vector
 */
template <class T>
std::vector<std::vector<std::vector<T> > > blitz_to_vector(const blitz::Array<T,3> & a){
  std::vector<std::vector<std::vector<T> > > v;
  for (int i=a.lbound(blitz::firstDim);i<=a.ubound(blitz::firstDim);i++) 
    v.push_back(blitz_to_vector(a(i,ALL,ALL)));
  return(v);
}

/** Translation blitz <-> vector
 */
template <class T>
std::vector<std::vector<std::vector<std::vector<T> > > > blitz_to_vector(const blitz::Array<T,4> & a){
  std::vector<std::vector<std::vector<std::vector<T> > > > v;
  for (int i=a.lbound(blitz::firstDim);i<=a.ubound(blitz::firstDim);i++) 
    v.push_back(blitz_to_vector(a(i,ALL,ALL,ALL)));
  return(v);
}

/** Translation blitz <-> vector
 */
template <class T>
std::vector<std::vector<std::vector<std::vector<std::vector<T> > > > >blitz_to_vector(const blitz::Array<T,5> & a){
  std::vector<std::vector<std::vector<std::vector<std::vector<T> > > > > v;
  for (int i=a.lbound(blitz::firstDim);i<=a.ubound(blitz::firstDim);i++) 
    v.push_back(blitz_to_vector(a(i,ALL,ALL,ALL,ALL)));
  return(v);
}


//-----------------------------------------------------------------------------------------------

/*
   Identity
*/
template<class T,class T2>
void set_Scalar( blitz::Array<T,2> & M, T2 x)
{
  assert(M.extent(blitz::secondDim) == M.extent(blitz::firstDim));
  M=0; for (int i=M.lbound(blitz::firstDim);i<=M.ubound(blitz::firstDim);i++)  M(i,i) = x;
}


//-----------------------------------------------------------------------------------------------
/**
 */
template<class T, int n>
inline blitz::Array<T,n> & operator <= (blitz::Array<T,n> & L, const blitz::Array<T,n> & R)
{L.resizeAndPreserve(R.shape()); L =R;return (L);} //{L.resize(R.shape()); L =R;return (L);}

/**
   Duplicate the array (two copies after one another)
*/
template<class T>
blitz::Array<T,1> duplicate(const blitz::Array<T,1> & M)
{
  int n= M.extent(0);
  blitz::Array<T,1> dup(2*n,blitz::fortranArray);
  dup(blitz::Range(1,n)) = M;
  dup(blitz::Range(n+1,2*n)) = M;
  return(dup);
}

//-----------------------------------------------------------------------------------------------

inline  int long_to_int(long x){ return int(x);}
//BZ_DECLARE_FUNCTION_RET( long_to_int, int);

//-----------------------------------------------------------------------------------------------

/**
   Return a constant array
*/
//template<class T,int n>
//Array<T,n> const_array(int i,T val) {blitz::Array<T,n> res(TinyVector<int,n>(i),fortranArray);res =val;return(res);}


// obsolete : to be removed ...
/**
   Return a constant array
*/
template<class T>
blitz::Array<T,1> const_array1(int i,T val) {blitz::Array<T,1> res(i,blitz::fortranArray);res =val;return(res);}

/**
   Return a constant array
*/
template<class T>
blitz::Array<T,2> const_array2(int i,T val) {blitz::Array<T,2> res(i,i,blitz::fortranArray);res =val;return(res);}

/**
   Return a constant array
*/
template<class T>
blitz::Array<T,3> const_array3(int i,T val) {blitz::Array<T,3> res(i,i,i,blitz::fortranArray);res =val;return(res);}

/**
   Return a scalar matrix
*/
template<class T>
blitz::Array<T,2> scalar_matrix(int n, T val) {blitz::Array<T,2> res(n,n,blitz::fortranArray);res =0;
  for (int i =1;i<=n;i++) res(i,i) = val; return(res);}

/**
   Return a diagonal matrix
*/
template<class T>
blitz::Array<T,2> diagonal_matrix(const blitz::Array<T,1> & diag) 
{int n=diag.extent(0); blitz::Array<T,2> res(n,n,blitz::fortranArray);res =0;
  for (int i =1;i<=n;i++) res(i,i) = diag(i); return(res);}

//-----------------------------------------------------------------------------------------------

/*
    Matrix multiplication 
*/
template<class T1,class T2, class T3>
inline void _for_matmul( const blitz::Array<T1,2> & M, const blitz::Array<T2,2> & B, blitz::Array<T3,2> & r)
{
  assert(M.extent(blitz::secondDim) == B.extent(blitz::firstDim));
  r=0;
  for (int i=M.lbound(blitz::firstDim),ir=r.lbound(blitz::firstDim);i<=M.ubound(blitz::firstDim);i++,ir++)
    for (int k=B.lbound(blitz::secondDim),kr=r.lbound(blitz::secondDim);k<=B.ubound(blitz::secondDim);k++,kr++)
      for (int j1=M.lbound(blitz::secondDim),j2=B.lbound(blitz::firstDim);j1<=M.ubound(blitz::secondDim);j1++,j2++) r(ir,kr)+= M(i,j1)*B(j2,k);
}

///
template<class T1,class T2, class T3>
inline void _for_matmul_v (const blitz::Array<T1,2> & M, const blitz::Array<T2,1> & V, blitz::Array<T3,1> & r)
{
  assert(M.extent(blitz::secondDim) == V.extent(blitz::firstDim));
  r=0;
  for (int i=M.lbound(blitz::firstDim),kr=r.lbound(blitz::firstDim);i<=M.ubound(blitz::firstDim);i++,kr++)
    for (int j1=M.lbound(blitz::secondDim),j2=V.lbound(blitz::firstDim);j1<=M.ubound(blitz::secondDim);j1++,j2++) r(kr)+= M(i,j1)*V(j2);
}


/**
    Matrix multiplication 
*/
template<class T>
inline blitz::Array<T,2> matmul(const blitz::Array<T,2> & M, const blitz::Array<T,2> & B, bool output_fortran_array=C_F_default) 
{
  blitz::GeneralArrayStorage<2> storage;if (output_fortran_array) storage = blitz::FortranArray<2>();
  blitz::Array<T,2> r(M.extent(blitz::firstDim),B.extent(blitz::secondDim),storage);
  if (M.extent(blitz::secondDim)>2) matmul_lapack(M,B,r); else _for_matmul(M,B,r);
  return(r);
}

/**
    Matrix multiplication : faster since no copy of result 
*/
template<class T>
inline void matmul(const blitz::Array<T,2> & M, const blitz::Array<T,2> & B, blitz::Array<T,2> & R) { matmul_lapack(M,B,R);}
  
/**
    Matrix multiplication 
*/
template<class T>
inline blitz::Array<std::complex<T>,2> matmul(const  blitz::Array<T,2> &  M, const blitz::Array<std::complex<T>,2> & B, bool output_fortran_array=C_F_default ) 
{ blitz::GeneralArrayStorage<2> storage;if (output_fortran_array) storage = blitz::FortranArray<2>();
  blitz::Array<std::complex<T>,2> cc(M.extent(blitz::firstDim),M.extent(blitz::secondDim),storage);cc = real_as_C(M);
  return(matmul(cc,B));}

/**
    Matrix multiplication 
*/
template<class T>
inline blitz::Array<std::complex<T>,2> matmul(const  blitz::Array<std::complex<T>,2> &  M, const blitz::Array<T,2> & B, bool output_fortran_array=C_F_default ) 
{ blitz::GeneralArrayStorage<2> storage;if (output_fortran_array) storage = blitz::FortranArray<2>();
  blitz::Array<std::complex<T>,2> cc(B.extent(blitz::firstDim),B.extent(blitz::secondDim),storage);cc = real_as_C(B);
  return(matmul(M,cc));}


///
template<class T>
inline blitz::Array<T,1> matmul(const blitz::Array<T,2> & M,const  blitz::Array<T,1> & X , bool output_fortran_array=C_F_default )
{ blitz::Array<T,1> r (M.extent(blitz::firstDim));  if (output_fortran_array) r.reindexSelf(1);
  _for_matmul_v(M,X,r);return(r);}


//-----------------------------------------------------------------------------------------------

/**
   Transposition : on site
*/
template<class T>
void transpose(blitz::Array<T,2> & M) { M.transposeSelf(blitz::secondDim,blitz::firstDim);}

/**
   Transposition
*/
template<class T>
void transpose(blitz::Array<T,3> & M) { M.transposeSelf(blitz::secondDim,blitz::firstDim,blitz::thirdDim);}


//-----------------------------------------------------------------------------------------------

/** 
    Scalar product
*/
template<class T>
inline T dot_product (const blitz::Array<T,1> & A, const blitz::Array<T,1> & B ) 
{  T r=0;for (int i=A.lbound(blitz::firstDim),j=B.lbound(0);i<=A.ubound(blitz::firstDim);i++,j++) 
	   r+= A(i)*B(j);return(r);}

/** 
    Scalar product
*/
template<class T>
inline T norm2 (const blitz::Array<T,1> & A) {return(sqrt(dot_product(A,A)));}


/*
inline double dot(const blitz::Array<double,1> & a, const blitz::Array<double,1> & b)
{
  int inc =1,N= a.size();
  return MYFORTRAN(ddot)(&N,ref<double,1>(a),&inc,ref<double,1>(b),&inc);
}

inline COMPLEX dot(const blitz::Array<COMPLEX,1> & a, const blitz::Array<COMPLEX,1> & b)
{
  int inc =1,N= a.size();
  return MYFORTRAN(zdot)(&N,ref<COMPLEX,1>(a),&inc,ref<COMPLEX,1>(b),&inc);
}
*/

//-----------------------------------------------------------------------------------------------

/** 
    Cross product. Dim 3 only
    If assert is active, it tests the dimension of the array
*/
template<class T>
inline blitz::Array<T,1> cross_product (const  blitz::Array<T,1> & A, const blitz::Array<T,1> & B, bool output_fortran_array=C_F_default)
{  
  blitz::Array<T,1> r(3);
  assert (A.extent(blitz::firstDim) == 3);
  assert (B.extent(blitz::firstDim) == 3);
  r=0;
  int ia = A.lbound(blitz::firstDim);int ib = B.lbound(blitz::firstDim);
  r(0) = A(ia+1)* B(ib+2) - B(ib+1) * A(ia+2);
  r(1) = - A(ia)* B(ib+2) + B(ib) * A(ia+2);
  r(2) = A(ia)*B(ib+1) - B(ib) * A(ia+1);
  if (output_fortran_array) r.reindexSelf(1);
  return(r);
}

//-----------------------------------------------------------------------------------------------

/** 
    Same as in Fortran
*/
int my_nint(double x);
//BZ_DECLARE_FUNCTION_RET(my_nint,int) 

/** 
    Real -> complex matrix
*/
std::complex<double> real_as_C (double x);
//BZ_DECLARE_FUNCTION_RET(real_as_C,std::complex<double>)


//-----------------------------------------------------------------------------------------------

/** 
    Inversion : 
    Uses explicit formula for dim =1,2 and lapack for higher dimension
*/
template<class T>  bool Inverse_Matrix (blitz::Array<T,2> & M)
{
  int i0 = M.lbound(blitz::firstDim),j0 = M.lbound(blitz::secondDim);
  T det,tmp;

  switch (M.rows()) {
  case(1) : 
    M(i0,j0) = T(1)/M(i0,j0);
    break;
  case(2) :
    det = M(i0,j0) * M(i0+1, j0+1) - M(i0+1,j0) * M(i0, j0+1);
    if (abs(det) < 1E-14) return(false);
    tmp = M(i0,j0); M(i0,j0) = M(i0+1, j0+1);   M(i0+1, j0+1) = tmp;
    M(i0+1,j0) *= -1;M(i0, j0+1) *=-1;
    M = M/det;
    break;
  default:
    return(Inverse_Matrix_with_lapack(M));
  }
  return (true);
}

/** Inverse but exit on failure */
template <class T>
inline void Inverse_Matrix_or_die (blitz::Array<T,2> & M, const char * s) { if (!Inverse_Matrix(M)) {std::cout << s;exit(1);} }


//-----------------------------------------------------------------------------------------------

double determinant (blitz::Array<double,2> & A);

//-----------------------------------------------------------------------------------------------

/** 
    Diagonalisation
    TO BE TESTED
    Input : M is the symmetric matrix
    Output : M is the eigenvectors and e the eigenvalues.
    It is a wrapper for dsyev if we use lapack.
*/
template<class T>  bool Diagonalise ( blitz::Array<T,2> & M, blitz::Array<T,1> & e, bool Compute_Eigenvectors=true)
{
  int i0 = M.lbound(blitz::firstDim),j0 = M.lbound(blitz::secondDim);
  T tmp,tmp2;

  switch (M.rows()) {
  case(1) : 
    e(i0) = M(i0,j0);
    M(i0,j0) = 1;
    break;
    /*
      not finished...
  case(2) : { 
      blitz::Array<T,2> MM(M);
    tmp2 = MM(1,1) - MM(0,0);
    tmp = sqrt(power(tmp2,2) + 4*MM(0,1)*MM(1,0) );
    e(0) = (MM(0,0) + MM(1,1) - tmp)/2;
    e(1) = (MM(0,0) + MM(1,1) + tmp)/2;
    MM(0,0) = tmp2 - tmp;
    MM(1,0) = 2*MM(0,1);
    MM(0,1) = tmp2 + tmp;
    MM(1,1) =  2*MM(0,1);}
    break;
    */
  default:
    return(  Diagonalise_with_lapack(M,e,Compute_Eigenvectors));
  }
 
  return (true);
}


//-----------------------------------------------------------------------------------------------


/** 
    A swap for matrices
*/
template<class T> void swap_matrix ( blitz::Array<T,2> & A,  blitz::Array<T,2> & B )
{
  assert(A.shape()==B.shape());
  T tmp;
  for (int i=A.lbound(blitz::firstDim),i2=B.lbound(blitz::firstDim);i<=A.ubound(blitz::firstDim);i++,i2++)
    for (int j=A.lbound(blitz::secondDim),j2=B.lbound(blitz::secondDim);j<=A.ubound(blitz::secondDim);j++,j2++)
      {tmp = A(i,j); A (i,j) = B(i2,j2); B(i2,j2) = tmp;}
}


/// 
template<int n>
bool operator==(const blitz::TinyVector<int, n> & a,const blitz::TinyVector<int, n> & b)
{ bool ok=true; for (int i=0;i<n;i++) ok= ok &&( a(i)==b(i));  return(ok);}

/// 
template<int n>
bool operator!=(const blitz::TinyVector<int, n> & a,const blitz::TinyVector<int, n> & b)
{ return (!(a==b));}


// Some custom storage

template <int n> class F_storage {};
template <> class F_storage<2> : public blitz::FortranArray<2>
{ public:   F_storage(int i, int j): blitz::FortranArray<2>() { ordering() = i-1,j-1;  }};
template <> class F_storage<3> : public blitz::FortranArray<3>
{ public:   F_storage(int i, int j, int k): blitz::FortranArray<3>() { ordering() = i-1,j-1,k-1;  }};
template <> class F_storage<4> : public blitz::FortranArray<4>
{ public:   F_storage(int i, int j, int k,int l): blitz::FortranArray<4>() { ordering() = i-1,j-1,k-1,l-1;  }};
template <> class F_storage<5> : public blitz::FortranArray<5>
{ public:   F_storage(int i, int j, int k,int l,int m): blitz::FortranArray<5>() { ordering() = i-1,j-1,k-1,l-1,m-1;  }};
template <> class F_storage<6> : public blitz::FortranArray<6>
{ public:   F_storage(int i, int j, int k,int l,int m,int n): blitz::FortranArray<6>() { ordering() = i-1,j-1,k-1,l-1,m-1,n-1;  }};
template <> class F_storage<7> : public blitz::FortranArray<7>
{ public:   F_storage(int i, int j, int k,int l,int m,int n,int p): blitz::FortranArray<7>() { ordering() = i-1,j-1,k-1,l-1,m-1,n-1,p-1;  }};


#ifdef HAVE_MPI
// --------------------------------------------
//  Additionnal MPI routine for some blitz arrays
// --------------------------------------------

/** */
template<class T,int n>
void  myMPI_bcast(blitz::Array<T,n> & a)
{ 
  int size = a.size();		 
  myMPI_barrier();
#ifdef MPI_BINDINGS_CXX
  MPI::COMM_WORLD.Bcast(ref<T,n>(a,false),size, mpi_datatype<T>::value, myMPI_master);
#else
  MPI_Bcast (ref<T,n>(a,false),size, mpi_datatype<T>::value, myMPI_master, MPI_COMM_WORLD); 
#endif
  myMPI_barrier();
}


/** */
template<class T,int n>
void  myMPI_reduce_sum(blitz::Array<T,n> & a,blitz::Array<T,n> & b)
{ 
  int size = a.size();assert(a.size()==b.size()); assert(a.ordering()==b.ordering()); 
  myMPI_barrier();
#ifdef MPI_BINDINGS_CXX
  MPI::COMM_WORLD.Reduce (ref<T,n>(a,false),ref<T,n>(b,false),size, mpi_datatype<T>::value, 
			  MPI_SUM, myMPI_master);
#else
  MPI_Reduce (ref<T,n>(a,false),ref<T,n>(b,false),size, mpi_datatype<T>::value, 
	      MPI_SUM, myMPI_master, MPI_COMM_WORLD); 
#endif
  myMPI_barrier();
}


/** */
template<class T,int n>
void  myMPI_reduce_sum_on_site(blitz::Array<T,n> & a)
{ 
  blitz::Array<T,n> b(a.copy());			  
  int size = a.size();
  myMPI_barrier();
#ifdef MPI_BINDINGS_CXX
  MPI::COMM_WORLD.Reduce (ref<T,n>(b,false),ref<T,n>(a,false),size, mpi_datatype<T>::value, 
			  MPI_SUM, myMPI_master);
#else
  MPI_Reduce (ref<T,n>(b,false),ref<T,n>(a,false),size, mpi_datatype<T>::value, 
	      MPI_SUM, myMPI_master, MPI_COMM_WORLD); 
#endif
  myMPI_barrier();
  //a=b; 
}

/** */
/*
template<class T,int n>
void myMPI_exchange(blitz::Array<T,n> & A,int from,blitz::Array<T,n> & B,int dest)
{ 
  assert(A.size()==B.size());					
  MPI_Status status; 
  if (myMPI_myid() == from)
    MPI_Send(ref<T,n>(A,false),A.size(), mpi_datatype<T>::value, dest, 99,MPI_COMM_WORLD);
  if (myMPI_myid() == dest)
    MPI_Recv(ref<T,n>(B,false),A.size(), mpi_datatype<T>::value, from, 99,MPI_COMM_WORLD,&status);
}
*/

/** */
template<class T,int n>
void myMPI_send(blitz::Array<T,n> & A, int dest)
{ 
  myMPI_barrier();
#ifdef MPI_BINDINGS_CXX
  MPI::COMM_WORLD.Send(ref<T,n>(A,false),A.size(), mpi_datatype<T>::value, dest, 99);
#else
  MPI_Send(ref<T,n>(A,false),A.size(), mpi_datatype<T>::value, dest, 99,MPI_COMM_WORLD);
#endif
}


/** */
template<class T,int n>
void myMPI_recv(blitz::Array<T,n> & A,int from)
{
  MPI_Status status; 
#ifdef MPI_BINDINGS_CXX
  MPI::COMM_WORLD.Recv(ref<T,n>(A,false),A.size(), mpi_datatype<T>::value, from, 99,&status);
#else
  MPI_Recv(ref<T,n>(A,false),A.size(), mpi_datatype<T>::value, from, 99,MPI_COMM_WORLD,&status);
#endif
}

/*
PROBLEM ...
   The ref is normally destroyed after the call.
   For non blocking communication, this can be a problem if
   copy has been used : when the recv is completed, the copy may no longer exists
   hence a segfault.
   So I assert that no copy has been made.
   IT IS ASSUMED THAT the array A will not be move in memory
*/

template<class T,int n>
void myMPI_irecv(blitz::Array<T,n> & A,int from,MPI_Request & request){
  ref<T,n> R(A,false);
  assert(!R.use_copy());
#ifdef MPI_BINDINGS_CXX
  request = MPI::COMM_WORLD.Irecv(R,A.size(), mpi_datatype<T>::value, from, 99);
#else
  MPI_Irecv(R,A.size(), mpi_datatype<T>::value, from, 99,MPI_COMM_WORLD,&request);
#endif
}



/** Scatter the array A to the array B on the nodes, by slicing the last dimension.
    FortranArray only.
    last dim of B must be declared blitz::Range(myMPI_slice_inf(1,NL),myMPI_slice_sup(1,NL)+1)
    because some node have one more point than other and scatter always dispatch
    j+1 
 */
template<class T,int n>
void myMPI_scatterv(const blitz::Array<T,n> & A,blitz::Array<T,n> & B)
{
  //assert(is_fortran_ordering(A));
  // Same formulas as for myMPI_slice.
  int imin=1,imax =A.extent(n-1);
  int j=(imax - imin + 1)/myMPI_number_proc();
  int i= imax - imin + 1 -myMPI_number_proc()*j;
  int *scounts = new int[myMPI_number_proc()];
  int *displs = new int[myMPI_number_proc()];
  for (int node=0; node<myMPI_number_proc(); ++node) 
    { 
      int inf=  node<=i-1 ?  imin + node*(j+1) :  imin + node*j + i;
      int sup = node<=i-1 ?  imin + (node+1)*(j+1) -1 : imin + (node+1)*j  + i - 1;
      displs[node] = inf-1; 
      scounts[node] = sup -inf +1; 
      //    cout<<"Node :"<<node<<" inf"<<inf<<"  sup"<<sup<<endl;
    } 
#ifdef MPI_BINDINGS_CXX
  MPI::COMM_WORLD.Scatterv(ref<T,n>(A), scounts, displs, mpi_datatype<T>::value,
			    ref<T,n>(B),j+1 , mpi_datatype<T>::value,  myMPI_master);
#else
  MPI_Scatterv(ref<T,n>(A), scounts, displs, mpi_datatype<T>::value,
	       ref<T,n>(B),j+1 , mpi_datatype<T>::value,  myMPI_master, MPI_COMM_WORLD);
#endif
  delete displs;delete scounts;
}

#endif
// end MPI 

/*@}*/

} // end namespace

#endif


