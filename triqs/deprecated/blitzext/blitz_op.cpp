
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

#include "blitz_op.hpp"
#include <triqs/utility/blas_headers.hpp>

#define COMPLEX std::complex<double>
using namespace std;
using namespace blitz;

namespace Blitz_OP {

int my_nint(double x){return(int(floor(x+0.5)));}
complex<double> real_as_C (double x){return(x);}

double sgn(complex<double> Z)
{
  if (real(Z)>0) return(1); else return(-1);
}

//---------------------------------------------------------------------------------

void matmul_lapack(const Array<double,2> & A, const Array<double,2> & B,Array<double,2> & RES)
{
  assert(A.extent(1) == B.extent(0));
  double alpha=1,beta =0;
  char c='N';
  int m1 = A.extent(0),m2 = A.extent(1),m3 = B.extent(1);
  MYFORTRAN(dgemm)(&c,&c,m1,m3,m2,alpha, ref<double,2>(A), m1, ref<double,2>(B), m2, beta, ref<double,2> (RES), m1);
}

//---------------------------------------------------------------------------------

void matmul_lapack(const Array<complex<double>,2> & A, const Array<complex<double>,2> & B, Array<complex<double>,2> & RES)
{
  assert(A.extent(1) == B.extent(0));
  complex<double> alpha=1,beta =0;
  char c='N';
  int m1 = A.extent(0),m2 = A.extent(1),m3 = B.extent(1);
  MYFORTRAN(zgemm)(&c,&c,m1,m3,m2,alpha, ref<complex<double>,2>(A), m1, ref<complex<double>,2> (B), m2, beta, ref<complex<double>,2> (RES), m1);
}

//---------------------------------------------------------------------------------


bool Inverse_Matrix_with_lapack (Array<double,2> & M)
{
  assert (M.rows() ==M.cols());
  int n = M.rows();
  ref<double,2> refM(M,false);//do not check the ordering since inverse and transposition commutes
  double *Mptr=refM;

  int *ipiv = new int[n*n];double *work = new double[n];int info;
  MYFORTRAN(dgetrf) (n,n,Mptr, n,ipiv,info);
  if (info) {
    cout<<"PB1 Inversion = "<<info<<endl<<flush;
    return (false);
  }
  MYFORTRAN(dgetri) (n,Mptr,n,ipiv,work,n,info);
  delete[] work; delete[] ipiv;
  if (info) cout<<"PB2 Inversion = "<<info<<endl<<flush;

  return ((info==0));
}


double Inverse_Matrix_with_lapack_and_Return_Determinant (double * refM, int n, int LDA, string message)
{
  int *ipiv = new int[n*n];double *work = new double[n];int info;
  MYFORTRAN(dgetrf) (n,n,refM, LDA,ipiv,info);
  if (info) {delete[] work; delete[] ipiv; FATAL(message<<"\nInverse_Matrix_and_Compute_Determinant : dgetrf pb  = "<<info);}
 
  // compute the det 
  double det =1;
  for (int i = 0; i<n; ++i) det *= refM[i*LDA +i];// the diagonal M[i,i]
  // compute the sign of the permutation
  bool flip=false;
  for (int i=0; i<n; i++)  if (ipiv[i]!=i+1) flip = !(flip);
 
  // finish the computation of the inverse
  MYFORTRAN(dgetri) (n,refM,LDA,ipiv,work,n,info);

  delete[] work; delete[] ipiv;
  if (info) FATAL(message<<"\nInverse_Matrix_and_Compute_Determinant : dgetrf pb  = "<<info);
  
  return (flip ? - det : det) ;
}


bool Inverse_Matrix_with_lapack ( Array<complex<double>,2> & M)
{
  assert (M.rows() ==M.cols());
  int n = M.rows();
  ref<complex<double>,2> refM(M,false);  //do not check the ordering since inverse and transposition commutes
  complex<double> *Mptr = refM;
 
  int *ipiv = new int[n*n];complex<double> *work = new complex<double>[n];int info;
  MYFORTRAN(zgetrf) (n,n,Mptr, n,ipiv,info);
  if (info) return (false);
  MYFORTRAN(zgetri) (n,Mptr ,n,ipiv,work,n,info);
  delete[] work; delete[] ipiv;
  return ((info==0));
}


/*
  bool Inverse_Matrix_with_lapack (double *Mptr, int n)
  {
  int *ipiv = new int[n*n];double *work = new double[n];int info;
  MYFORTRAN(dgetrf) (n,n,Mptr, n,ipiv,info);
  if (info) return (false);
  MYFORTRAN(dgetri) (n,Mptr,n,ipiv,work,n,info);
  delete[] work; delete[] ipiv;
  return ((info==0));
  }

  bool Inverse_Matrix_with_lapack (complex<double> *Mptr,int n)
  {
  int *ipiv = new int[n*n];complex<double> *work = new complex<double>[n];int info;
  MYFORTRAN(zgetrf) (n,n,Mptr, n,ipiv,info);
  if (info) return (false);
  MYFORTRAN(zgetri) (n,Mptr ,n,ipiv,work,n,info);
  delete[] work; delete[] ipiv; if (del) delete[] refM;
  return ((info==0));
  }
*/
   
   
//---------------------------------------------------------------------------------

bool A_inv_B_with_lapack  (Array<double,2> & A, Array<double,2> & B )
{
  assert(A.rows()==A.cols());
  int n = A.rows(); 
  int *ipiv = new int[n*n]; 
  char c='N';int info;
  ref<double,2> refA(A),refB(B);
  MYFORTRAN(dgetrf) (n,n,refA, n,ipiv,info);
  if (info) return (false);
  MYFORTRAN(dgetrs) (&c,n,n,refA,n,ipiv,refB,n,info);
  delete[] ipiv;
  return ((info==0)); 
}


bool A_inv_B_with_lapack  (Array<complex<double>,2> & A, Array<complex<double>,2> & B )
{
  assert(A.rows()==A.cols());
  int n = A.rows(); 
  int *ipiv = new int[n*n]; 
  char c='N';int info;
  ref<complex<double>,2> refA(A),refB(B);
  MYFORTRAN(zgetrf) (n,n,refA, n,ipiv,info);
  if (info) return (false);
  MYFORTRAN(zgetrs) (&c,n,n,refA,n,ipiv,refB,n,info);
  delete[] ipiv;
  return ((info==0));
}

//---------------------------------------------------------------------------------

double determinant (Array<double,2> & A)
{
  assert(A.rows()==A.cols());
  int n = A.rows(); 
  if (n==0) return 1;
  int *ipiv = new int[n]; 
  //char c='N';
  int info;
  ref<double,2> refA(A);

  MYFORTRAN(dgetrf) (n,n,refA, n,ipiv,info);
  assert(!info);

  double det =1;
  for (int i = A.lbound(0); i<=A.ubound(0); i++)
    det *= A(i,i);
  
  // compute the sign of the permutation
  bool flip=false;

  for (int i=0; i<n; i++)
    if (ipiv[i]!=i+1) flip = !(flip);
    
  delete[] ipiv;
  return (flip ? - det : det) ;
}


//---------------------------------------------------------------------------------


int Diagonalise_with_lapack(Array<double,2> & A, Array<double,1> & ev , bool Compute_Eigenvectors)
{
  assert(A.rows()==A.cols());
  int dim = A.rows(); 
  double *work;
  int info,lwork = 64*dim;
  char uplo='U',compz=(Compute_Eigenvectors ? 'V' : 'N') ;
  work=new double [lwork];
  ref<double,2> refA(A);
  ref<double,1> refD(ev);
  
  MYFORTRAN(dsyev) (&compz,&uplo,dim,refA,dim,refD,work,lwork,info);
  
  delete[] work;
  
  return (info);
}

//---------------------------------------------------------------------------------


int Diagonalise_with_lapack(Array<COMPLEX,2> & A, Array<double,1> & ev , bool Compute_Eigenvectors)
{
  assert(A.rows()==A.cols());
  int dim = A.rows(); 
  int info,lwork = 64*dim;
  char uplo='U',compz=(Compute_Eigenvectors ? 'V' : 'N') ;

  COMPLEX *work=new COMPLEX [lwork];
  double *work2=new double [lwork];
  
  //cout<<A<<ev<<endl;

  //  Array<double,2> t(1,1); t=8.2; Array<double,1> ee(1);
  //Diagonalise_with_lapack(t,ee,true);
  //cout<<t<<ee<<endl;
  
  MYFORTRAN(zheev) (&compz,&uplo,dim,ref<COMPLEX,2>(A),dim,ref<double,1>(ev),work,lwork,work2,info);
  
  //cout<<A<<ev<<endl;
  // assert(0);
  
  delete[] work;delete[] work2;

  if (info) cout<<"Error code zheev : "<<info<<endl;
  return (info);
}
  
  
//--------------------------------------------------------------

Array<COMPLEX,1> Root_Polynome_w_Hessenberg_Matrix(const Array<COMPLEX,1> & poly)  
{
  Array<COMPLEX,1> rac(poly.size()-1);
  char job='E',compz='N';
  int ilo =1 , n=poly.size()-1,info;
  int ihi=n,ldh=n,ldz=1,lwork=n*n;
  COMPLEX *Z=NULL;
  COMPLEX *work=new COMPLEX [lwork];
  
  Array<COMPLEX,2> M(n,n,fortranArray);
  M=0;
  for (int i = 1; i<n; i++) M(i+1,i) =1;
  for (int i = 1; i<=n; i++) M(i,n) = - poly(i)/poly(n+1);

  cout<<poly<<endl;
  cout<<M<<endl;

  MYFORTRAN(zhseqr)(&job,&compz,n,ilo,ihi,ref<COMPLEX,2>(M),ldh,ref<COMPLEX,1>(rac),Z,ldz,work,lwork,info);

  cout<<rac<<endl;

  
  delete[] work;
  return(rac);
}


//--------------------------------------------------------------

bool svd1_lapack(Array<COMPLEX,2> & M, Array<COMPLEX,1> & ev,Array<COMPLEX,2> & vl, Array<COMPLEX,2> & vr)
{
  assert (M.rows() ==M.cols());
  int n = M.rows(),info,lwork = 4 * n;
  COMPLEX *work = new COMPLEX[lwork];
  double *rwork = new double[2*n];
  char c='V';
  
  MYFORTRAN(zgeev)(&c,&c, n, ref<COMPLEX,2>(M), n, ref<COMPLEX,1>(ev), ref<COMPLEX,2>(vl), n, ref<COMPLEX,2>(vr), n,
		 work, lwork, rwork, info );

  vr = vl;
  Inverse_Matrix_with_lapack(vr);
  vr.transposeSelf(secondDim,firstDim);

  delete[] work; delete[] rwork;
  return ((info==0));
}


bool svd_lapack(Array<COMPLEX,2> & M, Array<double,1> & Sigma,Array<COMPLEX,2> & U, Array<COMPLEX,2> & VT)
{
  assert (M.rows() ==M.cols());
  int n = M.rows(),info,lwork = 6*n*n,lrwork = 10*n*n ;//4 * n;
  COMPLEX *work = new COMPLEX[lwork];
  double *rwork = new double[lrwork];
  int *iwork = new int[8*n];
  char c='A';
  
  // call and return the lwork. Refinement

  MYFORTRAN(zgesdd)(&c,n, n, ref<COMPLEX,2>(M), n, ref<double,1>(Sigma),ref<COMPLEX,2>(U), n, ref<COMPLEX,2>(VT), n, 
		    work, lwork, rwork, iwork, info );

  delete[] work; delete[] rwork; delete[] iwork;
  if (info) std::cout << "ZGESDD error :" << info << std::endl;
  return ((info==0));
}

}
