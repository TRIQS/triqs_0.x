
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

#ifndef MYBLAS_LAPACK_H
#define MYBLAS_LAPACK_H

// Will allow to adapt easily to compilers....
// emacs regepx : \(\w*\)_\(\W\) into MYFORTRAN(\1)\2
// http://jamesthornton.com/emacs/node/emacs_97.html

#include <triqs/utility/fortran_mangling.hpp>
#include <complex>

// Declaration of the lapack routines.
extern "C" 
{
 
  void MYFORTRAN(dgetrf) (int &,int &,double[], int & ,int [],int &);
  void MYFORTRAN(dgetri) (int & , double [],  int & ,int [],double [],int & ,int & );
  void MYFORTRAN(dgetrs) (char *,int & , int & , double [],  int & ,int [],double [],int & ,int & );
  
  void MYFORTRAN(zgetrf) (int &,int &,std::complex<double>[], int & ,int [],int &);
  void MYFORTRAN(zgetri) (int & , std::complex<double> [],  int & ,int [],std::complex<double> [],int & ,int & );
  void MYFORTRAN(zgetrs) (char *,int & , int & , std::complex<double> [],  int & ,int [],std::complex<double> [],int & ,int & );

  // C <- alpha AB + beta C
  void MYFORTRAN(dgemm) (char *, char *,  //TRANSA, TRANSB  
		       int & , //
		       int & , //
		       int & , //
		       double &, //ALPHA, 
		       double[], //  A, 
		       int &, // LDA 
		       double[], // B
		       int &, // LDB,
		       double &, // BETA, 
		       double[], // C
		       int & //  LDC
		       );

  void MYFORTRAN(zgemm) (char *, char *,  //TRANSA, TRANSB  
		       int & , //
		       int & , //
		       int & , //
		       std::complex<double> &, //ALPHA, 
		       std::complex<double>[], //  A, 
		       int &, // LDA 
		       std::complex<double>[], // B
		       int &, // LDB,
		       std::complex<double> &, // BETA, 
		       std::complex<double>[], // C
		       int & //  LDC
		       );

  void MYFORTRAN(zhseqr)( char *,char*, //JOB, COMPZ
			int &, // N
			int &, //ILO
			int &,// IHI 
			std::complex<double>[], // H
			int &, // LDH
			std::complex<double>[], //W
			std::complex<double>[], //Z
			int &, // LDZ
			std::complex<double>[], //WORK
			int &, // LWORK
			int & // INFO
			);

  void MYFORTRAN(dsyev)(char*,char*,        // JOBZ and UPLO
		      int &,              // Matrix Size
		      double[],            // matrix
		      int&,               // LDA of the matrix
		      double[],           // Eigenvalues array
		      double[],            // WORK
		      int&,               // LWORK
		      int &               // INFO
		      );


  void MYFORTRAN(zheev)(char*,char*,        // JOBZ and UPLO
		      int &,              // Matrix Size
		      std::complex<double> [],            // matrix
		      int&,               // LDA of the matrix
		      double[],           // Eigenvalues array
		      std::complex<double>[],            // WORK
		      int &,               // LWORK
		      double[],  // WORK2
		      int &               // INFO
		      );

  void MYFORTRAN(zgeev)(char*,char*,        //  JOBVL, JOBVR
		      int &,              // Matrix Size
		      std::complex<double> [],            // matrix
		      int&,               // LDA of the matrix
		      std::complex<double>[],   // Eigenvalues array
		      std::complex<double>[],   // VL
		      int &,               //LDVL
		      std::complex<double>[],   // VR
		      int &,               //LDVR
		      std::complex<double>[],   // WORK
		      int &,               // LWORK
		      double[],  // RWORK
		      int &               // INFO
		      );


  void MYFORTRAN(zgesdd)(char*,               //  JOBZ
			 int &,               // M
			 int &,               // N
			 std::complex<double> [],  // matrix
			 int&,                // LDA of the matrix
			 double[],            // Singular values array
			 std::complex<double>[],   // U
			 int &,               //LDU
			 std::complex<double>[],   // VT
			 int &,               //LDVT
			 std::complex<double>[],   // WORK
			 int &,               // LWORK
			 double[],            // RWORK
			 int [],              //IWORK
			 int &                // INFO
			 );

  double MYFORTRAN(dswap)(const int & N , double x[], const int& incx, double y[], const int& incy);
  double MYFORTRAN(dcopy)(const int & N , double x[], const int& incx, double y[], const int& incy);
  double MYFORTRAN(dscal)(const int & N , double & alpha,double x[], const int& incx);


  double MYFORTRAN(ddot)(const int & N, double  x[], const int& incx, double y[], const int & incy);
  double MYFORTRAN(zdot)(const int & N, std::complex<double>  x[], const int& incx, std::complex<double> y[], const int& incy);
  
  void MYFORTRAN(dtrmv)(const char & UPLO, const char & TRANS, const char & DIAG, 
			int & N, double A[], int & LDA, double X[], int & INCX);
  void MYFORTRAN(dgemv)(const char* trans, const int & m, const int & n, const double & alpha, const double A[], int & lda,
			const double x[], const int & incx, const double & beta, double y[], const int & incy);
  void MYFORTRAN(dger)(const int &, const int &, const double &, double [], const int &, double [], const int &, double [], const int &);

}

#endif

