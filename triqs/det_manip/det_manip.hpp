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
#ifndef TRIQS_DETMANIP_H
#define TRIQS_DETMANIP_H
#include <vector>
#include <iterator>
#include <triqs/arrays/matrix.hpp>
#include <triqs/arrays/mapped_functions.hpp>
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
    
    // serialized data
    value_type det;
    size_t Nmax,N, last_try;
    std::vector<size_t> row_num,col_num;
    std::vector<xy_type> x_values,y_values;
    int sign;
    matrix_type mat_inv;
    long long n_opts, n_opts_max_before_check;

   private:
    // temporary work data, not saved, serialized, etc....  
    struct work_data_type1 { 
     xy_type x, y;
     vector_type MB,MC, B, C;
     value_type ksi;
     size_t i,j,ireal,jreal;
     void reserve(size_t s) { B.resize(s); C.resize(s); MB.resize(s); MC.resize(s); MB()=0; MC()=0; }
    };
    
    struct work_data_type2 { 
     xy_type x[2], y[2];
     matrix_type MB,MC, B, C,ksi;
     size_t i[2],j[2],ireal[2],jreal[2];
     void reserve(size_t s) { MB.resize(s,2); MC.resize(2,s); B.resize(s,2), C.resize(2,s); ksi.resize(2,2); MB() = 0; MC() = 0; }
     value_type det_ksi() const { return ksi(0,0) * ksi(1,1) - ksi(1,0)* ksi(0,1);}
    };

    work_data_type1 w1;
    work_data_type2 w2;
    value_type newdet;
    int newsign;

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
     w1.reserve(Nmax); w2.reserve(Nmax);
    }

   private:
    void _construct_common() { 
     last_try=0; sign =1;
     n_opts=0; n_opts_max_before_check = 100;
    }

   public:
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
      reserve(init_size);
      mat_inv()=0;
      det = 1; 
     _construct_common();
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
      _construct_common();
      N =X.size(); 
      if (N==0) { det = 1; reserve(30); return;} 
      if (N>Nmax) reserve(2*N); // put some margin..
      std::copy(X.begin(),X.end(), std::back_inserter(x_values));
      std::copy(Y.begin(),Y.end(), std::back_inserter(y_values));
      mat_inv()=0;
      for (size_t i=0; i<N; ++i) { 
       row_num.push_back(i);col_num.push_back(i); 
       for (size_t j=0; j<N; ++j)
	mat_inv(i,j) = f(x_values[i],y_values[j]);
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
      ar & make_nvp("Nmax",Nmax) & make_nvp("N",N) 
       & make_nvp("n_opts",n_opts) & make_nvp("n_opts_max_before_check",n_opts_max_before_check) 
       & make_nvp("det",det) & make_nvp("sign",sign) 
       & make_nvp("Minv",mat_inv) 
       & make_nvp("row_num",row_num) & make_nvp("col_num",col_num) 
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
    value_type try_insert(size_t i, size_t j, xy_type const & x, xy_type const & y) {

     // check input and store it for complete_operation
     assert(i<=N);  assert(j<=N); assert(i>=0); assert(j>=0);
     if (N==Nmax) reserve(2*Nmax);
     last_try = 1;
     w1.i=i; w1.j=j; w1.x=x; w1.y = y;

     // treat empty matrix separately 
     if (N==0) { newdet = f(x,y); newsign = 1; return newdet; }

     // I add the row and col and the end. If the move is rejected,
     // no effect since N will not be changed : Minv(i,j) for i,j>=N 
     // has no meaning.
     for (size_t k= 0; k< N; k++) {
      w1.B(k) = f(x_values[k],y);
      w1.C(k) = f(x, y_values[k]);
     }
     range R(0,N);
     w1.MB(R) = mat_inv(R,R) * w1.B(R);
     w1.ksi = f(x,y) - boost::numeric::bindings::blas::dot( w1.C(R) , w1.MB(R) );
     newdet = det*w1.ksi;
     newsign = ((i + j)%2==0 ? sign : -sign);   // since N-i0 + N-j0  = i0+j0 [2]
     return (newdet/det)*(newsign*sign);          // sign is unity, hence 1/sign == sign
    } 

    //------------------------------------------------------------------------------------------
   private : 

    void complete_insert () {
     // store the new value of x,y. They are seen through the same permutations as rows and cols resp.
     x_values.push_back(w1.x); y_values.push_back(w1.y);
     row_num.push_back(0); col_num.push_back(0);

     // special empty case again
     if (N==0) { N=1; mat_inv(0,0) = 1/newdet; return; }

     range R1(0,N);  
     w1.MC(R1) = mat_inv(R1,R1).transpose() * w1.C(R1); 
     w1.MC(N) = -1;
     w1.MB(N) = -1;

     N++;

     // keep the real position of the row/col
     // since we insert a col/row, we have first to push the col at the right
     // and then say that col w1.i is stored in N, the last col.
     // same for rows
     for (int_type i =N-2; i>=int_type(w1.i); i--) row_num[i+1]= row_num[i];
     row_num[w1.i] = N-1;
     for (int_type i =N-2; i>=int_type(w1.j); i--) col_num[i+1]= col_num[i];
     col_num[w1.j] = N-1;

     // Minv is ok, we need to complete 
     w1.ksi = 1/w1.ksi;

     // compute the change to the inverse
     // M += w1.ksi w1.MB w1.MC with BLAS. first put the 0 
     range R(0,N);
     mat_inv(R,N-1) = 0;
     mat_inv(N-1,R) = 0;
     mat_inv(R,R) += triqs::arrays::a_x_ty(w1.ksi, w1.MB(R) ,w1.MC(R)) ;//mat_inv(R,R) += w1.ksi* w1.MB(R) * w1.MC(R)
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
     assert(i0!=i1); assert(j0!=j1);assert(i0<=N);  assert(j0<=N); assert(i0>=0); assert(j0>=0);
     assert(i1<=N+1);  assert(j1<=N+1); assert(i1>=0); assert(j1>=0);

     if (N >= Nmax) reserve(2*Nmax); // check this resize ... we add 2 lines
     last_try = 10;
     w2.i[0]=i0;w2.i[1]=i1; w2.j[0]=j0;w2.j[1]=j1;
     w2.x[0] = x0;w2.y[0] = y0; w2.x[1] = x1;w2.y[1] = y1;

     // w1.ksi = Delta(tau,tauP) - Cw.MB using BLAS
     w2.ksi(0,0) = f(x0,y0);
     w2.ksi(0,1) = f(x0,y1);
     w2.ksi(1,0) = f(x1,y0);
     w2.ksi(1,1) = f(x1,y1);

     // treat empty matrix separately 
     if (N==0) {
      newdet = w2.det_ksi(); 
      newsign = 1;
      return newdet;
     }

     // I add the rows and cols and the end. If the move is rejected,
     // no effect since N will not be changed : inv_mat(i,j) for i,j>=N has no meaning.
     for (size_t k= 0; k< N; k++) {
      w2.B(k,0) = f(x_values[k],y0);
      w2.B(k,1) = f(x_values[k],y1);
      w2.C(0,k) = f(x0, y_values[k]);
      w2.C(1,k) = f(x1, y_values[k]);
     }
     range R(0,N), R2(0,2);
     w2.MB(R,R2) = mat_inv(R,R) * w2.B(R,R2); 
     w2.ksi -= w2.C (R2, R) * w2.MB(R, R2);
     newdet = det * w2.det_ksi();
     newsign = ((i0 + j0 + i1 + j1)%2==0 ? sign : -sign); // since N-i0 + N-j0 + N + 1 -i1 + N+1 -j1 = i0+j0 [2]
     return (newdet/det)*(newsign*sign); // sign is unity, hence 1/sign == sign
    } 

    //------------------------------------------------------------------------------------------
   private:
    void complete_insert2 () {
     // store the new value of x,y. They are seen through the same permutations as rows and cols resp.
     for (int k=0; k<2; ++k) { 
      x_values.push_back(w2.x[k]); y_values.push_back(w2.y[k]);
      row_num.push_back(0); col_num.push_back(0);
     } 

     range R2(0,2); 
     if (N==0) { N=2; mat_inv(R2, R2) = inverse( w2.ksi); return; }

     range Ri(0,N);   
     w2.MC(R2,Ri) = w2.C(R2,Ri) * mat_inv(Ri,Ri);
     w2.MC(R2, range(N, N+2) ) = -1; // identity matrix 
     w2.MB(range(N,N+2), R2 ) = -1; // identity matrix ! 

     // keep the real position of the row/col
     // since we insert a col/row, we have first to push the col at the right
     // and then say that col w2.i[0] is stored in N, the last col.
     // same for rows
     for (int k =0; k<2; ++k) { 
      N++;
      for (int_type i =N-2; i>=int_type(w2.i[k]); i--) row_num[i+1]= row_num[i];
      row_num[w2.i[k]] = N-1;
      for (int_type i =N-2; i>=int_type(w2.j[k]); i--) col_num[i+1]= col_num[i];
      col_num[w2.j[k]] = N-1;
     }
     w2.ksi = inverse (w2.ksi);
     range R(0,N);
     mat_inv(R,range(N-2,N)) = 0;
     mat_inv(range(N-2,N),R) = 0;
     mat_inv(R,R) += w2.MB(R,R2) * (w2.ksi * w2.MC(R2,R)); 
    }

   public:
    //------------------------------------------------------------------------------------------

    /** 
     * Consider the removal the colj0 and row i0 from the matrix. 
     *
     * Returns the ratio of det Minv_new / det Minv.
     * This routine does NOT make any modification. It has to be completed with complete_operation().
     */
    value_type try_remove(size_t i, size_t j){
     assert(i<N);  assert(j<N); assert(i>=0); assert(j>=0);
     w1.i=i;w1.j=j;last_try = 2;
     w1.jreal = col_num[w1.j];
     w1.ireal = row_num[w1.i];
     // compute the newdet
     // first we resolve the w1.ireal,w1.jreal, with the permutation of the Minv, then we pick up what
     // will become the 'corner' coefficient, if the move is accepted, after the exchange of row and col.
     w1.ksi = mat_inv(w1.jreal,w1.ireal);
     newdet = det*w1.ksi;
     newsign = ((i + j)%2==0 ? sign : -sign);
     return (newdet/det)*(newsign*sign); // sign is unity, hence 1/sign == sign
    }
    //------------------------------------------------------------------------------------------
   private:
    void complete_remove() {
     if (N==1) { clear(); return; }

     // repack the matrix inv_mat
     // swap the rows w1.ireal and N, w1.jreal and N in inv_mat
     // Remember that for M row/col is interchanged by inversion, transposition.
     {
      range R(0,N);
      if (w1.jreal !=N-1){
       triqs::arrays::swap( mat_inv(w1.jreal,R), mat_inv(N-1,R));
       y_values[w1.jreal] = y_values[N-1]; 
      }

      if (w1.ireal !=N-1){
       triqs::arrays::swap (mat_inv(R,w1.ireal),  mat_inv(R,N-1));
       x_values[w1.ireal] = x_values[N-1];
      }
     }

     N--;

     // M <- a - d^-1 b c with BLAS
     w1.ksi = - 1/mat_inv(N,N);
     range R(0,N);

     mat_inv(R,R) += triqs::arrays::a_x_ty(w1.ksi,mat_inv(R,N),mat_inv(N,R));

     // modify the permutations
     for (size_t k =w1.i; k<N; k++) {row_num[k]= row_num[k+1];}
     for (size_t k =w1.j; k<N; k++) {col_num[k]= col_num[k+1];}
     for (size_t k =0; k<N; k++) 
     {
      if (col_num[k]==N) col_num[k]=w1.jreal;
      if (row_num[k]==N) row_num[k]=w1.ireal;
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

     assert(N>=2); assert(i0!=i1); assert(j0!=j1);
     assert(i0<N);  assert(j0<N); assert(i0>=0); assert(j0>=0); 
     assert(i1<N+1);  assert(j1<N+1);assert(i1>=0); assert(j1>=0);  

     last_try =11;

     w2.i[0]=std::min(i0,i1); 
     w2.i[1]=std::max(i0,i1); 
     w2.j[0]=std::min(j0,j1); 
     w2.j[1]=std::max(j0,j1); 
     w2.ireal[0] = row_num[w2.i[0]]; 
     w2.ireal[1] = row_num[w2.i[1]];
     w2.jreal[0] = col_num[w2.j[0]]; 
     w2.jreal[1] = col_num[w2.j[1]];

     // compute the newdet
     w2.ksi(0,0) = mat_inv(w2.jreal[0],w2.ireal[0]);
     w2.ksi(1,0) = mat_inv(w2.jreal[1],w2.ireal[0]);
     w2.ksi(0,1) = mat_inv(w2.jreal[0],w2.ireal[1]);
     w2.ksi(1,1) = mat_inv(w2.jreal[1],w2.ireal[1]);

     newdet = det * w2.det_ksi();
     newsign = ((i0 + j0+ i1 + j1)%2==0 ? sign : -sign);

     return (newdet/det)*(newsign*sign); // sign is unity, hence 1/sign == sign
    }
    //------------------------------------------------------------------------------------------
   private:
    void complete_remove2() {
     if (N==2) { clear(); return;} // put the sign to 1 also .... Change complete_remove... 

     size_t i_real_max =std::max(w2.ireal[0],w2.ireal[1]);  
     size_t i_real_min =std::min(w2.ireal[0],w2.ireal[1]);  
     size_t j_real_max =std::max(w2.jreal[0],w2.jreal[1]);  
     size_t j_real_min =std::min(w2.jreal[0],w2.jreal[1]); 

     range R(0,N);

     if (j_real_max != N-1) { 
      triqs::arrays::swap( mat_inv(j_real_max,R), mat_inv(N-1,R));
      y_values[ j_real_max ] = y_values[N-1];
     }
     if (j_real_min != N-2) { 
      triqs::arrays::swap( mat_inv(j_real_min,R), mat_inv(N-2,R));
      y_values[ j_real_min ] = y_values[N-2];
     }
     if (i_real_max != N-1) { 
      triqs::arrays::swap (mat_inv(R,i_real_max),  mat_inv(R,N-1));
      x_values[ i_real_max ] = x_values[N-1];
     }
     if (i_real_min != N-2) { 
      triqs::arrays::swap (mat_inv(R,i_real_min),  mat_inv(R,N-2));
      x_values[ i_real_min ] = x_values[N-2];
     }

     N -= 2;

     // M <- a - d^-1 b c with BLAS
     range Rn(0,N), Rl(N,N+2);
     //w2.ksi = mat_inv(Rl,Rl);
     //w2.ksi = inverse( w2.ksi);
     w2.ksi =  inverse( mat_inv(Rl,Rl));

     // write explicitely the second product on ksi for speed ?
     mat_inv(Rn,Rn) -= mat_inv(Rn,Rl) * (w2.ksi * mat_inv(Rl,Rn));

     // modify the permutations
     for (size_t k =w2.i[0]; k<w2.i[1]-1; k++)   row_num[k] = row_num[k+1];
     for (size_t k =w2.i[1]-1; k<N; k++)         row_num[k] = row_num[k+2];
     for (size_t k =w2.j[0]; k<w2.j[1]-1; k++)   col_num[k] = col_num[k+1];
     for (size_t k =w2.j[1]-1; k<N; k++)         col_num[k] = col_num[k+2];
     for (size_t k =0; k<N; k++) {
      if (col_num[k]==N+1)   col_num[k]=j_real_max;
      if (col_num[k]==N)     col_num[k]=j_real_min;
      if (row_num[k]==N+1)   row_num[k]=i_real_max;
      if (row_num[k]==N)     row_num[k]=i_real_min;
     }

     for (int u=0; u<2; ++u) { row_num.pop_back(); col_num.pop_back(); x_values.pop_back(); y_values.pop_back(); } 
    }
    //------------------------------------------------------------------------------------------

    /**
     * Consider the change the column j and the corresponding y.
     *
     * Returns the ratio of det Minv_new / det Minv.
     * This routine does NOT make any modification. It has to be completed with complete_operation().
     */
    value_type try_change_col(size_t j, xy_type const & y) {
     assert(j<N); assert(j>=0);
     w1.j=j;last_try = 3;
     w1.jreal = col_num[j];
     w1.y = y;

     // Compute the col B.
     for (size_t i= 0; i<N;i++) w1.MC(i) = f(x_values[i] , w1.y) - f(x_values[i], y_values[w1.jreal]);
     range R(0,N);
     w1.MB(R) = mat_inv(R,R) * w1.MC(R);

     // compute the newdet
     w1.ksi = (1+w1.MB(w1.jreal));
     newdet = det*w1.ksi;
     newsign = sign;
     return (newdet/det)*(newsign*sign); // sign is unity, hence 1/sign == sign
    }
    //------------------------------------------------------------------------------------------
   private:
    void complete_change_col() {
     range R(0,N);
     y_values[w1.jreal] = w1.y;

     // modifying M : Mij += w1.ksi Bi Mnj 
     // using Shermann Morrison formula.
     // implemented in 2 times : first Bn=0 so that Mnj is not modified ! and then change Mnj
     // Cf notes : simply multiply by -w1.ksi
     w1.ksi = - 1/(1+ w1.MB(w1.jreal));
     w1.MB(w1.jreal) = 0;
     mat_inv(R,R) += triqs::arrays::a_x_ty(w1.ksi,w1.MB(R), mat_inv(w1.jreal,R));
     mat_inv(w1.jreal,R)*= -w1.ksi; 
    }

    //------------------------------------------------------------------------------------------
   public:

    /**
     * Consider the change the row i and the corresponding x.
     *
     * Returns the ratio of det Minv_new / det Minv.
     * This routine does NOT make any modification. It has to be completed with complete_operation().
     */
    value_type try_change_row(size_t i, xy_type const & x) {
     assert(i<N); assert(i>=0);
     w1.i=i;last_try = 4;
     w1.ireal = row_num[i];
     w1.x = x;

     // Compute the col B.
     for (size_t i= 0; i<N;i++) w1.MB(i) = f(w1.x, y_values[i] ) -  f(x_values[w1.ireal], y_values[i] ); 
     range R(0,N);
     w1.MC(R) = mat_inv(R,R).transpose() * w1.MB(R); 

     // compute the newdet
     w1.ksi = (1+w1.MC(w1.ireal));
     newdet = det*w1.ksi;
     newsign = sign;
     return (newdet/det)*(newsign*sign); // sign is unity, hence 1/sign == sign
    }
    //------------------------------------------------------------------------------------------
   private:
    void complete_change_row() {
     range R(0,N);
     x_values[w1.ireal] = w1.x;

     // modifying M : M ij += w1.ksi Min Cj
     // using Shermann Morrison formula.
     // impl. Cf case 3
     w1.ksi = - 1/(1+ w1.MC(w1.ireal));
     w1.MC(w1.ireal) = 0;
     mat_inv(R,R) += triqs::arrays::a_x_ty(w1.ksi,mat_inv(R,w1.ireal),w1.MC);
     mat_inv(R,w1.ireal) *= -w1.ksi;
    }
    //------------------------------------------------------------------------------------------
   private: 

    bool check_mat_inv (double precision) const { 
     if (N==0) return true;
     const bool relative = true;
     matrix_type res(N,N);
     for (size_t i=0; i<N;i++)
      for (size_t j=0; j<N;j++)
       res(i,j) = f(x_values[i], y_values[j]);
     res = inverse(res); 
     value_type r = max_element (abs(res - mat_inv(N,N)));
     value_type r2= max_element (abs(res + mat_inv(N,N)));
     return ( r >  (relative ? precision * r2 : precision));
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
       last_try=0; return; 
       break; // double call of complete_operation... 
      default: 
       TRIQS_RUNTIME_ERROR<< "Misuing det_manip";
     }
     det = newdet;
     sign = newsign;
     last_try =0;
     ++n_opts;
     if (n_opts > n_opts_max_before_check) { 
      if (!check_mat_inv(1.e-8))
       TRIQS_RUNTIME_ERROR << "Deviation too large ";
     }
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
