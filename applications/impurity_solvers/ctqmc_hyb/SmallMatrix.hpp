
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

#ifndef SMALL_MATRIX_H_jkj
#define SMALL_MATRIX_H_jkj
 
#define TRIQS_FORTRAN( id ) id##_


/** 
    Class of small matrices. 
    Optimised for calculation with small matrices.
*/


// There are two possible orderings: by lines or by columns
enum OrderTypeSM {ByLines,ByColumns};


// A generic SmallMatrix Class
template<typename VAL, OrderTypeSM ORDER>
class SmallMatrix {};


// The template specialization for Line Ordering
template<typename VAL>
class SmallMatrix<VAL, ByLines> {

  friend class SmallMatrix<VAL,ByColumns>;

public:
  const int n1, n2;

protected:
  VAL * restrict data;
  const bool is11;

public:

  // Constructors
  SmallMatrix(int n1_, int n2_, VAL * restrict data_) : n1(n1_), n2(n2_), data(data_), is11((n1==1) && (n2==1)) {};

  SmallMatrix(int n1_, int n2_, std::vector<VAL> & data_) : n1(n1_), n2(n2_), data(&data_[0]), is11((n1==1) && (n2==1)) {};

  SmallMatrix(const SmallMatrix & M) : n1(M.n1), n2(M.n2), data(M.data), is11(M.is11) {}

  // Access to the i,j element
  inline VAL operator() (int i, int j) const { return data[i*n2 + j]; }
  inline VAL & operator() (int i, int j) { return data[i*n2 + j]; } 
  inline const VAL * restrict data_() const {return data;}

  // The = operator (makes a copy of M in *this)
  SmallMatrix<VAL,ByLines> & operator = (const SmallMatrix<VAL,ByLines> & M) {
    assert (n1==M.n1); assert(n2==M.n2);
    data = M.data;
    return *this;
  }

  // Set to the product of two small matrices
  void setTo_AB(const SmallMatrix<VAL,ByLines> &, const SmallMatrix<VAL,ByColumns> &);

  void setTo_ADB(const SmallMatrix<VAL,ByLines> &, const std::vector<VAL> &, const SmallMatrix<VAL,ByColumns> &);

  void setFrom_A(const SmallMatrix<VAL,ByLines> & A) {
    assert(n1==A.n1); assert(n2==A.n2);
    for (int i = 0; i < n1*n2; ++i) data[i] = A.data[i];
  }

  void setFrom_DA(const SmallMatrix<VAL,ByLines> & A, const std::vector<VAL> & D) {
    assert(n1==A.n1); assert(n2==A.n2);
    for (int i = 0; i < n1; ++i)
      for (int j = 0; j < n2; ++j)
          data[i*n2+j] = A.data[i*n2+j]*D[i];
  }

  void setTo_AB_blas(const SmallMatrix<VAL,ByLines> &, const SmallMatrix<VAL,ByColumns> &);

};

// multiplication with a loop A \times B
template<typename VAL>
void SmallMatrix<VAL, ByLines>::setTo_AB(const SmallMatrix<VAL,ByLines> & A, const SmallMatrix<VAL,ByColumns> & B) {

  if (A.is11 && B.is11)
    data[0] = A.data[0] * B.data[0];
  else {
    const int nn1(n1), nn2(n2), nn3(A.n2);
    int idx(0), k1(0), k2(0);
    for (int i=0; i<nn1; ++i) {
      for (int j=0; j<nn2; ++j) {
        VAL tmp(0);
        k2=j*nn3;
        for (int k1k=k1, k2k=k2; k1k<k1+nn3; ++k1k, ++k2k) { 
          tmp += A.data[k1k]*B.data[k2k];
        }
        data[idx] = tmp;
        ++idx;
      }
      k1+=nn3;
    }
  }
}

// multiplication with a loop A \times D \times B
template<typename VAL>
void SmallMatrix<VAL, ByLines>::setTo_ADB(const SmallMatrix<VAL,ByLines> & A, const std::vector<VAL> & D, const SmallMatrix<VAL,ByColumns> & B) {

  const VAL * restrict pD(&D[0]);

  if (A.is11 && B.is11)
    data[0] = A.data[0] * D[0] * B.data[0];
  else {
    const int nn1(n1), nn2(n2), nn3(A.n2);
    int idx(0), k1(0), k2(0);
    for (int i=0; i<nn1; ++i) {
      for (int j=0; j<nn2; ++j) {
        VAL tmp(0);
        k2=j*nn3;
        for (int k=0, k1k=k1, k2k=k2; k<nn3; ++k, ++k1k, ++k2k) { 
          tmp += A.data[k1k]*pD[k]*B.data[k2k];
        }
        data[idx] = tmp;
        ++idx;
      }
      k1+=nn3;
    }
  }
}


// multiplication using blas
template<typename VAL>
void SmallMatrix<VAL, ByLines>::setTo_AB_blas(const SmallMatrix<VAL,ByLines> & A, const SmallMatrix<VAL,ByColumns> & B) {
  VAL alpha=1, beta =0;
  char c='N', ct='T';
  int m1 = n1, m2 = n2;
  int m3 = A.n2;
  TRIQS_FORTRAN(dgemm)(&ct,&c,m2,m1,m3,alpha,(VAL *)B.data,m3,(VAL *)A.data,m3,beta,data,m2);
}




// The template specialization for Column Ordering
template<typename VAL>
class SmallMatrix<VAL, ByColumns> {

  friend class SmallMatrix<VAL,ByLines>;

public:
  const int n1, n2;

protected:
  VAL * restrict data;
  const bool is11;

public:

  // Constructors
  SmallMatrix(int n1_, int n2_, VAL * restrict data_) : n1(n1_), n2(n2_), data(data_), is11((n1==1) && (n2==1)) {};

  SmallMatrix(int n1_, int n2_, std::vector<VAL> & data_) : n1(n1_), n2(n2_), data(&data_[0]), is11((n1==1) && (n2==1)) {};

  SmallMatrix(const SmallMatrix & M) : n1(M.n1), n2(M.n2), data(M.data), is11(M.is11) {}

  // Access to the i,j element
  inline VAL operator() (int i, int j) const { return data[j*n1 + i]; }
  inline VAL & operator() (int i, int j) { return data[j*n1 + i]; }

  // The = operator (makes a copy of M in *this)
  SmallMatrix<VAL,ByColumns> & operator = (const SmallMatrix<VAL,ByColumns> & M) {
    assert (n1==M.n1); assert(n2==M.n2);
    data = M.data;
    return *this;
  }

  // Set to the product of two small matrices
  void setTo_AB(const SmallMatrix<VAL,ByLines> &, const SmallMatrix<VAL,ByColumns> &);

  void setTo_ADB(const SmallMatrix<VAL,ByLines> &, const std::vector<VAL> &, const SmallMatrix<VAL,ByColumns> &);

  void setFrom_A(const SmallMatrix<VAL,ByColumns> & A) {
    assert(n1==A.n1); assert(n2==A.n2);
    for (int i = 0; i < n1*n2; ++i) data[i] = A.data[i];
  }

  void setFrom_DA(const SmallMatrix<VAL,ByColumns> & A, const std::vector<VAL> & D) {
    assert(n1==A.n1); assert(n2==A.n2);
    for (int j = 0; j < n2; ++j)
      for (int i = 0; i < n1; ++i)
          data[j*n1+i] = A.data[j*n1+i]*D[i];
  }

  void setTo_AB_blas(const SmallMatrix<VAL,ByLines> &, const SmallMatrix<VAL,ByColumns> &);

};


// multiplication with a loop A \times B
template<typename VAL>
void SmallMatrix<VAL, ByColumns>::setTo_AB(const SmallMatrix<VAL,ByLines> & A, const SmallMatrix<VAL,ByColumns> & B) {

  if (A.is11 && B.is11)
    data[0] = A.data[0] * B.data[0];
  else {
    const int nn1(n1), nn2(n2), nn3(A.n2);
    int idx(0), k1(0), k2(0);
    for (int j=0; j<nn2; ++j) {
      for (int i=0; i<nn1; ++i) {
        VAL tmp(0);
        k1=i*nn3;
        for (int k1k=k1, k2k=k2; k1k<k1+nn3; ++k1k, ++k2k) { 
          tmp += A.data[k1k]*B.data[k2k];
        }
        data[idx] = tmp;
        ++idx;
      }
      k2+=nn3;
    }
  }
}


// multiplication with a loop A \times D \times B
template<typename VAL>
void SmallMatrix<VAL, ByColumns>::setTo_ADB(const SmallMatrix<VAL,ByLines> & A, const std::vector<VAL> & D, const SmallMatrix<VAL,ByColumns> & B) {

  const VAL * restrict pD(&D[0]);

  if (A.is11 && B.is11)
    data[0] = A.data[0]  * D[0] * B.data[0];
  else {
    const int nn1(n1), nn2(n2), nn3(A.n2);
    int idx(0), k1(0), k2(0);
    for (int j=0; j<nn2; ++j) {
      for (int i=0; i<nn1; ++i) {
        VAL tmp(0);
        k1=i*nn3;
        for (int k=0, k1k=k1, k2k=k2; k<nn3; ++k, ++k1k, ++k2k) { 
          tmp += A.data[k1k]*pD[k]*B.data[k2k];
        }
        data[idx] = tmp;
        ++idx;
      }
      k2+=nn3;
    }
  }
}


// multiplication using blas
template<typename VAL>
void SmallMatrix<VAL, ByColumns>::setTo_AB_blas(const SmallMatrix<VAL,ByLines> & A, const SmallMatrix<VAL,ByColumns> & B) {
  VAL alpha=1, beta =0;
  char c='N', ct='T';
  int m1 = n1, m2 = n2;
  int m3 = A.n2;
  TRIQS_FORTRAN(dgemm)(&ct,&c,m1,m2,m3,alpha,(VAL *)A.data,m3,(VAL *)B.data,m3,beta,data,m1);
}

#endif
