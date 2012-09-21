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
#ifndef TRIQS_DETMANIP_H
#define TRIQS_DETMANIP_H
#include <vector>
#include <iterator>
#include <triqs/arrays/matrix.hpp>
#include <triqs/arrays/linalg/inverse.hpp>
#include <triqs/arrays/linalg/a_x_ty.hpp>
#include <triqs/arrays/linalg/matmul.hpp>
#include <triqs/arrays/linalg/mat_vec_mul.hpp>
#include <boost/numeric/bindings/blas/level1/dot.hpp>
#include <triqs/arrays/proto/matrix_algebra.hpp>

namespace triqs { namespace det_manip { 

 /**
 */ 
 template<typename FunctionTypeArg>
  class det_manip {
   public:

    typedef typename boost::unwrap_reference<FunctionTypeArg>::type  FunctionType;
    typedef typename FunctionType::result_type                       value_type;
    typedef typename FunctionType::argument_type                     xy_type;
    typedef triqs::arrays::vector<value_type>                        vector_type;
    typedef triqs::arrays::matrix<value_type>                        matrix_type;
    typedef triqs::arrays::matrix_view<value_type>                   matrix_view_type;

   protected: // the data
    typedef std::ptrdiff_t int_type;
    typedef triqs::arrays::range range;

    FunctionType f;
    value_type det,newdet;
    size_t Nmax,N, last_try;
    std::vector<size_t> row_num,col_num;
    std::vector<xy_type> x_values,y_values;
    int sign,newsign;
    matrix_type mat_inv;

   private:
    struct work_data_type { // temporary work data, not saved, serialized, etc....  
     xy_type x, y, x1, y1;
     vector_type MB,MC, col, row;
     matrix_type MB2,MC2, B2, C2, ksi2;
     size_t i0,j0,ireal0,jreal0, i1,j1,ireal1,jreal1; 
     value_type ksi;

     void reserve(size_t s) { 
      col.resize(s); row.resize(s);
      MB.resize(s); MC.resize(s); MB2.resize(s,2); MC2.resize(2,s); B2.resize(s,2), C2.resize(2,s);
      ksi2.resize(2,2); 
      MB()=0; MC()=0; MB2() = 0; MC2() = 0; 
     }
    };

    work_data_type w;

   public:

    /** 
     * Like for std::vector, reserve memory for a bigger size.
     * Preserves only the matrix, not the temporary working vectors/matrices, so do NOT use it 
     * between a try_XXX and acomplete_operation
     */
    void reserve (size_t new_size) { 
     if (new_size <= Nmax) return;
     matrix_type Mcopy(mat_inv);
     size_t N0 = Nmax; Nmax = new_size;
     mat_inv.resize(Nmax,Nmax); mat_inv(range(0,N0), range(0,N0)) = Mcopy; // keep the content of mat_inv ---> into the lib ?
     row_num.reserve(Nmax);col_num.reserve(Nmax); x_values.reserve(Nmax);y_values.reserve(Nmax);
     w.reserve(Nmax);
    }

    /** 
     * \brief Constructor.
     * 
     * \param F         The function (NB : a copy is made of the F object in this class). // CHANGE this with boost::ref ? useful ?
     * \param init_size The maximum size of the matrix before a resize (like reserve in std::vector).
     *                  Like std::vector, resize is automatic (by a factor 2) but can yield a performance penalty
     *                  if it happens too often.
     */
    det_manip(FunctionTypeArg F,size_t init_size):
     f(boost::unwrap_ref(F)), Nmax(0) , N(0){ 
      this->reserve(init_size);
      mat_inv()=0;
      last_try=0; det = 1; sign =1;
     }

    /** \brief Constructor.
     * 
     * \param F         The function (NB : a copy is made of the F object in this class). // CHANGE this with boost::ref ? useful ?
     * \tparam ArgumentContainer 
     * \param X, Y : container for X,Y. 
     */
    template<typename ArgumentContainer1, typename ArgumentContainer2>
     det_manip(FunctionTypeArg F, ArgumentContainer1 const & X, ArgumentContainer2 const & Y) : f(boost::unwrap_ref(F)) { 
      if (X.size() != Y.size()) TRIQS_RUNTIME_ERROR<< " X.size != Y.size";
      last_try=0; sign =1;
      N =X.size(); 
      if (N==0) { det = 1; reserve(30); return;} 
      if (N>Nmax) reserve(2*N); // put some margin..
      std::copy(X.begin(),X.end(), std::back_inserter(x_values));
      std::copy(Y.begin(),Y.end(), std::back_inserter(y_values));
      mat_inv()=0;
      for (size_t i=0; i<N; ++i) { 
       row_num.push_back(i);col_num.push_back(i); 
       for (size_t j=0; j<N; ++j)
	mat_inv(i,j) = this->f(x_values[i],y_values[j]);
      }
      mat_inv = inverse(mat_inv);
      det = determinant(mat_inv);
     }

    /// Put to size 0 : like a vector 
    void clear () { 
     N = 0; sign = 1;det =1; last_try = 0;
     row_num.clear(); col_num.clear(); x_values.clear(); y_values.clear(); 
    }

    //----------------------- READ ACCESS TO DATA ----------------------------------

    /// Current size of the matrix
    size_t size() const { return N;}

    /// Returns the i-th values of x
    xy_type const & get_x(size_t i) const { return x_values[row_num[i]];}

    /// Returns the j-th values of y
    xy_type const & get_y(size_t j) const { return y_values[col_num[j]];}

    /** det M of the current state of the matrix.  */
    value_type determinant() const {return sign*det;}

    /** Returns M^{-1}(i,j) */
    value_type inverse_matrix(size_t i,size_t j) const {return mat_inv(col_num[i],row_num[j]);} // warning : need to invert the 2 permutations.

    /// Returns the inverse matrix. Warning : this is slow, since it create a new copy, and reorder the lines/cols
    matrix_view_type inverse_matrix() const {
     matrix_type res(N,N);
     for (size_t i=0; i<N;i++)
      for (size_t j=0; j<N;j++)
       res(i,j) = inverse_matrix(i,j);
     return res;
    }

    /// Rebuild the matrix. Warning : this is slow, since it create a new matrix and re-evaluate the function. 
    matrix_view_type matrix() const {
     matrix_type res(N,N);
     for (size_t i=0; i<N;i++)
      for (size_t j=0; j<N;j++)
       res(i,j) = f(get_x(i), get_y(j));
     return res;
    }

   private:
    //  ------------     BOOST Serialization ------------
    friend class boost::serialization::access;
    template<class Archive>
     void serialize(Archive & ar) {
      using boost::serialization::make_nvp;
      ar & make_nvp("Minv",mat_inv) & make_nvp("Nmax",Nmax) & make_nvp("row_num",row_num) & make_nvp("col_num",col_num) 
       & make_nvp("x_values",x_values) & make_nvp("y_values",y_values); 
     }

   public:

    // ------------------------- OPERATIONS -----------------------------------------------

    /**
     * Insert operation at colum j0 and row i0.
     *
     * The operation consists in adding :
     *
     *    * a column  f(x_i,    y_{j0}) 
     *    * and a row f(x_{i0}, x_j)
     *
     * The new colum/row will be at col j0, row i0.
     *
     * 0 <= i0,j0 <= N, where N is the current size of the matrix.
     * The current column j0 (resp. row i0) will become column j0+1 (resp. row i0+1).
     * Inserting at N simply add the new col at the end.

     * Returns the ratio of det Minv_new / det Minv.
     *
     * This routine does NOT make any modification. It has to be completed with complete_operation().
     */
    value_type try_insert(size_t i0, size_t j0, xy_type const & x, xy_type const & y) {

     // check input and store it for complete_operation
     w.i0=i0;w.j0=j0; last_try = 1;
     w.x = x;w.y = y;
     if (N==Nmax) reserve(2*Nmax);
     assert(i0<=N);  assert(j0<=N); assert(i0>=0); assert(j0>=0);

     // treat empty matrix separately (blas on vectors of size 0 ??)
     if (N==0) { newdet = f(x,y); newsign = 1; return newdet; }

     // I add the row and col and the end. If the move is rejected,
     // no effect since N will not be changed : Minv(i,j) for i,j>=N 
     // has no meaning.
     for (size_t k= 0; k< N; k++) {
      w.col(k) = f(x_values[k],y);
      w.row(k) = f(x, y_values[k]);
     }

     triqs::arrays::range R(0,N);
     w.MB(R) = mat_inv(R,R) * w.col(R);
     w.MB(N) = -1;

     // w.ksi = Delta(tau,tauP) - Cw.MB using BLAS
     w.ksi = f(x,y);
     w.row(N) = w.ksi;  
     w.ksi -= boost::numeric::bindings::blas::dot( w.row(R) , w.MB(R) );

     // compute the newdet
     newdet = det*w.ksi;
     newsign = ((i0 + j0)%2==0 ? sign : -sign);   // since N-i0 + 1 + N-j0 +1 = i0+j0 [2]
     return (newdet/det)*(newsign*sign);          // sign is unity, hence 1/sign == sign
    } 

    //------------------------------------------------------------------------------------------
   private : 

    void complete_insert () {
     assert(x_values.size()==N); assert(y_values.size()==N);
     assert( row_num.size()==N); assert( col_num.size()==N); 
     // store the new value of x,y. They are seen through the same permutations as rows and cols resp.
     x_values.push_back(w.x); y_values.push_back(w.y);
     row_num.push_back(0); col_num.push_back(0);

     // special empty case again
     if (N==0) { N=1; mat_inv(0,0) = 1/newdet; return; }

     triqs::arrays::range R1(0,N);  
     w.MC(R1) = mat_inv(R1,R1).transpose() * w.row(R1); 
     w.MC(N) = -1;

     N++;

     // keep the real position of the row/col
     // since we insert a col/row, we have first to push the col at the right
     // and then say that col w.i0 is stored in N, the last col.
     // same for rows
     for (int_type i =N-2; i>=int_type(w.i0); i--) row_num[i+1]= row_num[i];
     row_num[w.i0] = N-1;
     for (std::ptrdiff_t i =N-2; i>=int_type(w.j0); i--) col_num[i+1]= col_num[i];
     col_num[w.j0] = N-1;

     // Minv is ok, we need to complete 
     w.ksi = 1/w.ksi;

     // compute the change to the inverse
     // M += w.ksi w.MB w.MC with BLAS. first put the 0 
     triqs::arrays::range R(0,N);
     mat_inv(R,N-1) = 0;
     mat_inv(N-1,R) = 0;
     mat_inv(R,R) += triqs::arrays::a_x_ty(w.ksi, w.MB(R) ,w.MC(R)) ;//mat_inv(R,R) += w.ksi* w.MB(R) * w.MC(R)
    }

   public : 
    //------------------------------------------------------------------------------------------

    /**
     * Double Insert operation at colum j0,j1 and row i0,i1.
     *
     * The operation consists in adding :
     *    * a column  f(x_i,    y_{j0}) 
     *    * and a row f(x_{i0}, x_j)
     * The new colum/row will be at col j0, row i0.
     *
     * 0 <= i0,i1,j0,j1 <= N+1, where N is the current size of the matrix.

     * Returns the ratio of det Minv_new / det Minv.
     * This routine does NOT make any modification. It has to be completed with complete_operation().
     */

    value_type try_insert2(size_t i0, size_t i1, size_t j0, size_t j1, xy_type const & x0, xy_type const & x1, xy_type const & y0, xy_type const & y1) {
     // check input and store it for complete_operation
     w.i0=i0;w.i1 = i1;
     w.j0=j0;w.j1 = j1;
     last_try = 10;
     w.x = x0;w.y = y0; w.x1 = x1;w.y1 = y1;
     if (N < Nmax) reserve(2*Nmax); // check this resize ... we add 2 lines
     assert(i0<=N);  assert(j0<=N); assert(i0>=0); assert(j0>=0);
     assert(i1<=N+1);  assert(j1<=N+1); assert(i1>=0); assert(j1>=0);

     // w.ksi = Delta(tau,tauP) - Cw.MB using BLAS
     w.ksi2(0,0) = f(x0,y0);
     w.ksi2(0,1) = f(x0,y1);
     w.ksi2(1,0) = f(x1,y0);
     w.ksi2(1,1) = f(x1,y1);

     // treat empty matrix separately (blas on vectors of size 0 ??)
     if (N==0) {
      newdet = w.ksi2(0,0) * w.ksi2(1,1) - w.ksi2(1,0)* w.ksi2(0,1);
      newsign = 1;
      return newdet;
     }

     // I add the row and col and the end. If the move is rejected,
     // no effect since N will not be changed : Minv(i,j) for i,j>=N 
     // has no meaning.
     for (size_t k= 0; k< N; k++) {
      w.B2(k,0) = f(x_values[k],y0);
      w.B2(k,1) = f(x_values[k],y1);
      w.C2(0,k) = f(x0, y_values[k]);
      w.C2(1,k) = f(x1, y_values[k]);
     }

     triqs::arrays::range R(0,N);
     w.MB2(R,range(0,2)) = mat_inv(R,R) * w.B2(R,range(0,2)); 
     w.MB2(range(N,N+2), range(N, N+2) ) = -1; // identity matrix 

     w.ksi2 -= w.C2 (range(0,2), R) * w.MB2(R, range(0,2));

     // compute the newdet
     newdet = det * ( w.ksi2(0,0) * w.ksi2(1,1) - w.ksi2 (1,0)* w.ksi2(0,1));

     newsign = ((i0 + j0 + i1 + j1)%2==0 ? sign : -sign); // since N-i0 + 1 + N-j0 +1 = i0+j0 [2]
     return (newdet/det)*(newsign*sign); // sign is unity, hence 1/sign == sign
    } 

    //------------------------------------------------------------------------------------------
   private:
    void complete_insert2 () {
     assert(x_values.size()==N); assert(y_values.size()==N);
     assert( row_num.size()==N); assert( col_num.size()==N); 
     // store the new value of x,y. They are seen through the same permutations as rows and cols resp.
     x_values.push_back(w.x);  y_values.push_back(w.y);
     x_values.push_back(w.x1); y_values.push_back(w.y1);
     row_num.push_back(0); col_num.push_back(0);
     row_num.push_back(0); col_num.push_back(0);

     // special empty case again
     if (N==0) { N=2; mat_inv(range(0,2), range(0,2)) = inverse( w.ksi2); return; }

     triqs::arrays::range R1(0,N);  
     w.MC2(range(0,2),R1) = w.C2(range(0,2),R1) * mat_inv(R1,R1);

     w.MC2(range(N,N+2), range(N, N+2) ) = -1; // identity matrix 

     // keep the real position of the row/col
     // since we insert a col/row, we have first to push the col at the right
     // and then say that col w.i0 is stored in N, the last col.
     // same for rows
     N++;
     for (int_type i =N-2; i>=int_type(w.i0); i--) row_num[i+1]= row_num[i];
     row_num[w.i0] = N-1;
     for (int_type i =N-2; i>=int_type(w.j0); i--) col_num[i+1]= col_num[i];
     col_num[w.j0] = N-1;

     N++;
     for (int_type i =N-2; i>=int_type(w.i1); i--) row_num[i+1]= row_num[i];
     row_num[w.i1] = N-1;
     for (int_type i =N-2; i>=int_type(w.j1); i--) col_num[i+1]= col_num[i];
     col_num[w.j1] = N-1;

     // Minv is ok, we need to complete 
     w.ksi2 = inverse (w.ksi2);

     // compute the change to the inverse
     // M += w.ksi w.MB w.MC with BLAS. first put the 0 
     triqs::arrays::range R(0,N);
     mat_inv(R,range(N-2,N)) = 0;
     mat_inv(range(N-2,N),R) = 0;
     mat_inv(R,R) += w.MB2 * (w.ksi2 * w.MC2); 
    }

   public:
    //------------------------------------------------------------------------------------------

    /** 
     * Consider the removal the colj0 and row i0 from the matrix. 
     *
     * Returns the ratio of det Minv_new / det Minv.
     * This routine does NOT make any modification. It has to be completed with complete_operation().
     */
    value_type try_remove(size_t i0, size_t j0){
     assert(i0<N);  assert(j0<N); assert(i0>=0); assert(j0>=0);
     w.i0=i0;w.j0=j0;last_try = 2;
     w.jreal0 = col_num[w.j0];
     w.ireal0 = row_num[w.i0];
     // compute the newdet
     // first we resolve the w.ireal0,w.jreal0, with the permutation of the Minv, then we pick up what
     // will become the 'corner' coefficient, if the move is accepted, after the exchange of row and col.
     // See complete_operation
     w.ksi = mat_inv(w.jreal0,w.ireal0);
     newdet = det*w.ksi;
     newsign = ((i0 + j0)%2==0 ? sign : -sign);
     return (newdet/det)*(newsign*sign); // sign is unity, hence 1/sign == sign
    }
    //------------------------------------------------------------------------------------------
   private:
    void complete_remove() {
     if (N==1) {
      N=0;
      row_num.pop_back(); col_num.pop_back();
      x_values.pop_back(); y_values.pop_back();
      det =1;
      return;
     }

     // repack the matrix _Minv
     // swap the rows w.ireal0 and N, w.jreal0 and N in _Minv
     // Remember that for M row/col is interchanged by inversion, transposition.
     {
      triqs::arrays::range R(0,N);
      if (w.jreal0 !=N-1){
       triqs::arrays::swap( mat_inv(w.jreal0,R), mat_inv(N-1,R));
       y_values[w.jreal0] = y_values[N-1]; 
      }

      if (w.ireal0 !=N-1){
       triqs::arrays::swap (mat_inv(R,w.ireal0),  mat_inv(R,N-1));
       x_values[w.ireal0] = x_values[N-1];
      }
     }

     N--;

     // M <- a - d^-1 b c with BLAS
     w.ksi = - 1/mat_inv(N,N);
     triqs::arrays::range R(0,N);

     mat_inv(R,R) += triqs::arrays::a_x_ty(w.ksi,mat_inv(R,N),mat_inv(N,R));

     // modify the permutations
     // see case 1
     for (size_t k =w.i0; k<N; k++) {row_num[k]= row_num[k+1];}
     for (size_t k =w.j0; k<N; k++) {col_num[k]= col_num[k+1];}
     for (size_t k =0; k<N; k++) 
     {
      if (col_num[k]==N) col_num[k]=w.jreal0;
      if (row_num[k]==N) row_num[k]=w.ireal0;
     }
     row_num.pop_back(); col_num.pop_back();
     x_values.pop_back(); y_values.pop_back();
    }

   public:
    //------------------------------------------------------------------------------------------

    /** 
     * Double Removal operation of cols j0,j1 and rows i0,i1
     *
     * Returns the ratio of det Minv_new / det Minv.
     * This routine does NOT make any modification. It has to be completed with complete_operation().
     */
    value_type try_remove2(size_t i0, size_t i1, size_t j0, size_t j1) { 
     assert(i0<N);  assert(j0<N); assert(i0>=0); assert(j0>=0);
     w.i0=i0; w.i1=i1;
     w.j0=j0; w.j1=j1;
     last_try =11;
     w.w.jreal0 = col_num[w.j0];
     w.w.jreal1 = col_num[w.j1];
     w.ireal0 = row_num[w.i0];
     w.ireal1 = row_num[w.i1];
     // compute the newdet
     // first we resolve the w.ireal0,w.jreal0, with the permutation of the Minv, then we pick up what
     // will become the 'corner' coefficient, if the move is accepted, after the exchange of row and col.
     // See complete_operation
     for (size_t u =0; u<2; ++u)
      for (size_t v =0; v<2; ++v)
       w.ksi2(u,v) = mat_inv(w.jreal0_2[u],w.ireal0_2[v]);

     newdet = det * ( w.ksi2(0,0) * w.ksi2(1,1) - w.ksi2 (1,0)* w.ksi2(0,1));
     newsign = ((i0 + j0+ i1 + j1)%2==0 ? sign : -sign);
     return (newdet/det)*(newsign*sign); // sign is unity, hence 1/sign == sign
    }
    //------------------------------------------------------------------------------------------
   private:
    void complete_remove2() {
     if (N==2) { N=0; }

     else { 

      size_t i_real_max =std::max(w.ireal0,w.ireal1);  
      size_t i_real_min =std::min(w.ireal0,w.ireal1);  
      size_t j_real_max =std::max(w.jreal0,w.jreal1);  
      size_t j_real_min =std::min(w.jreal0,w.jreal1); 

      y_values[ j_real_max ] = y_values[N-1];
      y_values[ j_real_min ] = y_values[N-2];

      x_values[ i_real_max ] = x_values[N-1];
      x_values[ i_real_min ] = x_values[N-2];

      N =-2;

      // M <- a - d^-1 b c with BLAS
      triqs::arrays::range Rn(0,N), R2(N-2,N);
      w.ksi2 = - inverse( mat_inv(R2,R2)); 

      mat_inv(Rn,Rn) += mat_inv(Rn,R2) * w.ksi2 * mat_inv(R2,Rn); 

      // modify the permutations
      // see case 1
      size_t i_sh=0, j_sh =0;
      for (size_t k =0; k<N; k++) {
       if ((k==w.i0) || (k==w.i1)) i_sh++;
       if ((k==w.j0) || (k==w.j1)) j_sh++;  
       row_num[k]= row_num[k+i_sh];
       col_num[k]= col_num[k+j_sh];
       if (col_num[k]==N)   col_num[k]=j_real_max;
       if (col_num[k]==N-1) col_num[k]=j_real_min;
       if (row_num[k]==N)   row_num[k]=i_real_max;
       if (row_num[k]==N-1) row_num[k]=i_real_min;
      }

      for (int u=0; u<2; ++u) { row_num.pop_back(); col_num.pop_back(); x_values.pop_back(); y_values.pop_back(); } 
     }
    }
    //------------------------------------------------------------------------------------------

    /**
     * Consider the change the column j0 and the corresponding y.
     *
     * Returns the ratio of det Minv_new / det Minv.
     * This routine does NOT make any modification. It has to be completed with complete_operation().
     */
    value_type try_change_col(size_t j0, xy_type const & y) {
     assert(j0<N); assert(j0>=0);
     w.j0=j0;last_try = 3;
     w.jreal0 = col_num[j0];
     w.y = y;

     // Compute the col B.
     for (size_t i= 0; i<N;i++) w.MC(i) = f(x_values[i] , w.y) - f(x_values[i], y_values[w.jreal0]);
     triqs::arrays::range R(0,N);
     w.MB(R) = mat_inv(R,R) * w.MC(R);

     // compute the newdet
     w.ksi = (1+w.MB(w.jreal0));
     newdet = det*w.ksi;
     newsign = sign;
     return (newdet/det)*(newsign*sign); // sign is unity, hence 1/sign == sign
    }
    //------------------------------------------------------------------------------------------
   private:
    void complete_change_col() {
     triqs::arrays::range R(0,N);
     y_values[w.jreal0] = w.y;

     // modifying M : Mij += w.ksi Bi Mnj 
     // using Shermann Morrison formula.
     // implemented in 2 times : first Bn=0 so that Mnj is not modified ! and then change Mnj
     // Cf notes : simply multiply by -w.ksi
     w.ksi = - 1/(1+ w.MB(w.jreal0));
     w.MB(w.jreal0) = 0;
     mat_inv(R,R) += triqs::arrays::a_x_ty(w.ksi,w.MB(R), mat_inv(w.jreal0,R));
     mat_inv(w.jreal0,R)*= -w.ksi; 
    }

    //------------------------------------------------------------------------------------------
   public:

    /**
     * Consider the change the row i0 and the corresponding x.
     *
     * Returns the ratio of det Minv_new / det Minv.
     * This routine does NOT make any modification. It has to be completed with complete_operation().
     */
    value_type try_change_row(size_t i0, xy_type const & x) {
     assert(i0<N); assert(i0>=0);
     w.i0=i0;last_try = 4;
     w.ireal0 = row_num[i0];
     w.x = x;


     // Compute the col B.
     for (size_t i= 0; i<N;i++) w.MB(i) = f(w.x, y_values[i] ) -  f(x_values[w.ireal0], y_values[i] ); 
     triqs::arrays::range R(0,N);
     w.MC(R) = mat_inv(R,R).transpose() * w.MB(R); 

     // compute the newdet
     w.ksi = (1+w.MC(w.ireal0));
     newdet = det*w.ksi;
     newsign = sign;
     return (newdet/det)*(newsign*sign); // sign is unity, hence 1/sign == sign
    }
    //------------------------------------------------------------------------------------------
   private:
    void complete_change_row() {
     triqs::arrays::range R(0,N);
     x_values[w.ireal0] = w.x;

     // modifying M : M ij += w.ksi Min Cj
     // using Shermann Morrison formula.
     // impl. Cf case 3
     w.ksi = - 1/(1+ w.MC(w.ireal0));
     w.MC(w.ireal0) = 0;
     mat_inv(R,R) += triqs::arrays::a_x_ty(w.ksi,mat_inv(R,w.ireal0),w.MC);
     mat_inv(R,w.ireal0) *= -w.ksi;
    }
    //------------------------------------------------------------------------------------------
   public:
    /**
     *  Finish the move of the last try_xxx called.
     *  Throws if no try_xxx has been done or if the last operation was complete_operation.
     */
    void complete_operation() {
     switch(last_try){
      case(1):
       complete_insert();
       break;
      case(2): 
       complete_remove();
       break;
      case(3): 
       complete_change_col();
       break;
      case(4): 
       complete_change_row();
       break;
      case(10): 
       complete_insert2();
       break;
      case(11): 
       complete_remove2();
       break;
      case(0): 
       break; // double call of complete_operation... 
      default: 
       TRIQS_RUNTIME_ERROR<< "Misuing det_manip";
     }
     det = newdet;
     sign = newsign;
     last_try =0;
    }

    ///
    enum RollDirection {None,Up, Down,Left,Right};

    /**
     * "Cyclic Rolling" of the determinant.
     *
     * Right : Move the Nth col to the first col cyclically.
     * Left  : Move the first col to the Nth, cyclically.
     * Up    : Move the first row to the Nth, cyclically.
     * Down  : Move the Nth row to the first row cyclically.
     *
     * Returns -1 is the roll changes the sign of the det, 1 otherwise
     * NB : this routine is not a try_xxx : it DOES make the modification and does not need to be completed...
     * WHY is it like this ???? : try_roll : return det +1/-1.
     */
    int roll_matrix(RollDirection roll) {
     size_t tmp;
     const int_type NN=N;
     switch (roll) {
      case(None) : 
       return 1;
      case(Down) : 
       tmp = row_num[N-1];
       for (int_type i =NN-2; i>=0; i--) row_num[i+1]= row_num[i];
       row_num[0] = tmp;
       break;
      case(Up) : 
       tmp = row_num[0];
       for (int_type i =0; i<N-1; i++) row_num[i]= row_num[i+1];
       row_num[N-1] = tmp;
       break;
      case(Right) : 
       tmp = col_num[N-1];
       for (int_type i =NN-2; i>=0; i--) col_num[i+1]= col_num[i];
       col_num[0] = tmp;
       break;
      case(Left): 
       tmp = col_num[0];
       for (int_type i =0; i<N-1; i++) col_num[i]= col_num[i+1];
       col_num[N-1] = tmp;
       break;
      default:
       assert(0);
     }
     // signature of the cycle of order N : (-1)^(N-1)
     if ((N-1)%2==1) { sign *=-1; return -1;}
     return 1;
    }
  };
}}
#endif
