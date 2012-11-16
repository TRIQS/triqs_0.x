/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2012 by M. Ferrero, O. Parcollet
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
#ifndef TRIQS_GF_LOCAL_TAIL_H
#define TRIQS_GF_LOCAL_TAIL_H
#include "../tools.hpp"

namespace triqs { namespace gf { namespace local {

 namespace tqa= triqs::arrays; namespace tql= triqs::clef;  namespace mpl= boost::mpl; 
 typedef std::complex<double> dcomplex;

 class tail;                            // the value class
 class tail_view;                       // the view class

 template<typename G> struct LocalTail  : mpl::false_{};  // a boolean trait to identify the objects modelling the concept LocalTail
 template<> struct LocalTail<tail >     : mpl::true_{};
 template<> struct LocalTail<tail_view >: mpl::true_{};

 // a trait to find the scalar of the algebra i.e. the true scalar and the matrix ...
 template <typename T> struct is_scalar_or_element : mpl::or_< tqa::ImmutableMatrix<T>, tup::is_in_ZRC<T> > {};

 // ---------------------- implementation --------------------------------

 /// A common implementation class. Idiom : ValueView
 template<bool IsView> class tail_impl  {  
  public : 
   typedef void has_view_type_tag; // Idiom : ValueView 
   typedef tail_view view_type;
   typedef tail      non_view_type;

   typedef arrays::Option::C storage_order;
   typedef arrays::array      <dcomplex,3,storage_order>                         data_non_view_type;
   typedef arrays::array_view <dcomplex,3,storage_order>                         data_view_type;
   typedef typename mpl::if_c<IsView, data_view_type, data_non_view_type>::type  data_type;

   typedef arrays::matrix_view<dcomplex,       storage_order>  mv_type;
   typedef arrays::matrix_view<const dcomplex, storage_order>  const_mv_type;

   data_view_type data_view()             { return data;} 
   const data_view_type data_view() const { return data;}

   int order_min() const {return omin;}
   int order_max() const {return omin + size() -1;}
   size_t size()   const {return data.shape()[2];} 

   typedef tqa::mini_vector<size_t,2> shape_type;
   shape_type shape() const { return shape_type(data.shape()[0], data.shape()[1]);} 
   size_t shape(int i) const { return data.shape()[i];} 

   bool is_decreasing_at_infinity() const { return (order_min() >=1);}

  protected:
   int omin;
   data_type data;

   tail_impl(): omin(0),data() {} // all arrays of zero size (empty)
   tail_impl (size_t N1, size_t N2, int order_min, size_t size_) : omin(order_min), data(tqa::make_shape(N1,N2,size_)){data()=0;}
   tail_impl( data_type const &d, int order_min ): omin(order_min), data(d){}
   tail_impl(tail_impl          const & x): omin(x.omin), data(x.data){}
   tail_impl(tail_impl<!IsView> const & x): omin(x.omin), data(x.data){}

   friend class tail_impl<!IsView>;
  public:

   void operator = (tail_impl          const & x) { omin = x.omin; data = x.data;}
   void operator = (tail_impl<!IsView> const & x) { omin = x.omin; data = x.data;}

   mv_type operator() (int n)       { 
    if (n>this->order_max()) TRIQS_RUNTIME_ERROR<<" n > Max Order. n= "<<n <<", Max Order = "<<order_max() ; 
    if (n<this->order_min()) TRIQS_RUNTIME_ERROR<<" n < Min Order. n= "<<n <<", Min Order = "<<order_min() ;
    return this->data(tqa::range(), tqa::range(), n- this->order_min());
   }

   const_mv_type operator() (int n) const {  
    if (n>this->order_max()) TRIQS_RUNTIME_ERROR<<" n > Max Order. n= "<<n <<", Max Order = "<<order_max() ;
    if (n<this->order_min())  { mv_type::non_view_type r(this->shape()); r()=0; return r;} 
    return this->data(tqa::range(), tqa::range(), n- this->order_min());
   }

   operator freq_infty() const { return freq_infty();}

   /// Save in txt file : doc the format  ? ---> prefer serialization or hdf5 !
   void save(std::string file,  bool accumulate=false) const {}

   /// Load from txt file : doc the format ?
   void load(std::string file){}

   /// 
   friend void h5_write (tqa::h5::group_or_file fg, std::string subgroup_name, tail_impl const & t) {
    BOOST_AUTO( gr , fg.create_group(subgroup_name) );
    // Add the attribute
    h5_write(gr,"omin",t.omin);
    h5_write(gr,"data",t.data);
   }

   friend void h5_read  (tqa::h5::group_or_file fg, std::string subgroup_name, tail_impl & t){
    BOOST_AUTO( gr,  fg.open_group(subgroup_name) );
    // Check the attribute or throw 
    h5_read(gr,"omin",t.omin);
    h5_read(gr,"data",t.data);
   }

   //  BOOST Serialization
   friend class boost::serialization::access;
   template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
     ar & boost::serialization::make_nvp("omin",omin);
     ar & boost::serialization::make_nvp("data",data);
    }

   friend std::ostream & operator << (std::ostream & out, tail_impl const & x) { 
    out <<" tail/tail_view : Ordermin/max = "<< x.order_min() << "  "<< x.order_max();
    for (int u =x.order_min();u<= x.order_max(); ++u) out <<"\n ...  Order "<<u << " = " << x(u);
    return out; 
   } 

   /*  friend tail operator * (dcomplex a, tail_view const& r);
       friend tail operator * (tail_view const & r, dcomplex a);
       friend tail operator * (tail_view const & l, tail_view const& r);
       friend tail operator / (tail_view const & l, tail_view const& r);
       friend tail operator / (tail_view const & r, dcomplex a);
       friend tail operator + (tail_view const & l, tail_view const& r);
       friend tail operator - (tail_view const & l, tail_view const& r);
       */
 };

 // -----------------------------
 ///The View class of GF
 class tail_view : public tail_impl <true> { 
  typedef tail_impl <true>  B;
#ifdef TRIQS_ALLOW_EMPTY_VIEW
  public:
  tail_view ():B(){}
#endif
  public :
  template<bool V> tail_view(tail_impl<V> const & t): B(t){}
  tail_view(B::data_type const &d, int order_min): B(d,order_min){}
  void rebind( tail_view const &X) { omin = X.omin; data.rebind(X.data);}

  using B::operator=; // import operator = from impl. class or the default = is synthetized 
  tail_view & operator=(const tail_view & rhs) { 
   if (this->data.is_empty()) rebind(rhs); else B::operator=(rhs); return *this; 
  }
  tail_view & operator=(std::complex<double> const & x) { this->data = x; return *this;} 
  using B::operator(); // import all previously defined operator() for overloading
  TRIQS_CLEF_ADD_LAZY_CALL_WITH_VIEW(1,tail_view); // add lazy call, using the view in the expression tree.
  friend std::ostream & triqs_nvl_formal_print(std::ostream & out, tail_view const & x) { return out<<"tail_view";}

  void print_me() const { std::cout  << *this << std::endl ; }
 };

 // -----------------------------
 ///The regular class 
 class tail : public tail_impl <false> { 
  typedef tail_impl <false>  B;
  public : 
  tail():B() {} 
  typedef tqa::mini_vector<size_t,2> shape_type;
  tail(size_t N1, size_t N2,  int order_min=-1, size_t size_ =5):B(N1,N2,order_min,size_) {} 
  tail(shape_type const & sh,  int order_min=-1, size_t size_=5 ):B(sh[0],sh[1],order_min,size_) {} 
  tail(tail const & g): B(g){}
  tail(tail_view const & g): B(g){} 
  //template<typename GfType> tail(GfType const & x): B() { *this = x;} // to maintain value semantics

  using B::operator=;
  using B::operator();
  TRIQS_CLEF_ADD_LAZY_CALL_WITH_VIEW(1,tail_view); // add lazy call, using the view in the expression tree.

  /// The simplest tail corresponding to  : omega
  static tail_view omega(size_t N1, size_t N2, size_t size_) { tail t(N1,N2,-1, size_); t(t.order_min())=1; return t; }

  /// The simplest tail corresponding to  : omega, constructed from a shape for convenience
  static tail_view omega(tail::shape_type const & sh, size_t size_) { return omega(sh[0],sh[1],size_);}

 };

 /// Slice in orbital space
 template<bool V> tail_view slice_target(tail_impl<V> const & t, tqa::range R1, tqa::range R2) { 
  return tail_view(t.data_view()(R1,R2,tqa::range()),t.order_min());
 } 

 // various operations in a simple, non proto way
 inline tail inverse(tail_view const & t ) {
  // find in t
  int omin1 = - t.order_min();
  int si = t.size();
  tail t_inv(t.shape(), omin1, si);
  const int omin = t_inv.order_min(); //const int omax = t_inv.order_max();
  t_inv(omin) = inverse(t(t.order_min()));
  for (int n=1; n< int(t_inv.size());n++) {
  //for (size_t n=1; n< t_inv.size();n++) {
   // 0<= n-p < t.size() ---> p <= n, ok  and p > n- t.size() 
   const int pmin = std::max(0, n - int(t.size())+1 );
   // workaround on the bug on calling lapack with views ....
   // ?? fix this in the lib ??
   for (int p=pmin; p< n ; p++) { t_inv(omin + n) -= triqs::arrays::matrix< dcomplex >(t_inv(omin + p)*t(t.order_min() + n-p)); }
   t_inv(omin + n) *= t_inv(omin);
  }
  return t_inv;
 }

 inline tail mult_impl(tail_view const & l, tail_view const& r) {
  if (l.shape()[1] != r.shape()[0]) TRIQS_RUNTIME_ERROR<< "tail multiplication : shape mismatch";
  tail res(l.shape()[0], r.shape()[1],l.order_min() + r.order_min(), std::min(l.size(),r.size()));
  res.data_view() =0;
  for (int n =res.order_min(); n<=res.order_max(); ++n) {
   // sum_{p}^n a_p b_{n-p}. p <= a.n_max, p >= a.n_min and n-p <=b.n_max and n-p >= b.n_min
   // hence p <= min ( a.n_max, n-b.n_min ) and p >= max ( a.n_min, n- b.n_max)  
   const int pmin = std::max(l.order_min(), n - r.order_max() );
   const int pmax = std::min(l.order_max(), n - r.order_min() );
   for (int p = pmin; p <= pmax; ++p)  { res(n) += l(p) * r(n-p);} 
  }
  return res;
 }

 template<typename T1, typename T2> 
  TYPE_ENABLE_IF(tail,mpl::and_<LocalTail<T1>, LocalTail<T2>>)
  operator* (T1 const & a, T2 const & b) { return mult_impl(a,b); }

 template<typename T1, typename T2> TYPE_ENABLE_IF(tail,mpl::and_<tqa::ImmutableMatrix<T1>, LocalTail<T2>>)
  operator* (T1 const & a, T2 const & b) { 
   tail res(b); for (size_t i=0; i<res.size(); ++i) res(i)=a*res(i); return res;
  }

 template<typename T1, typename T2> TYPE_ENABLE_IF(tail,mpl::and_<LocalTail<T1>, tqa::ImmutableMatrix<T2>>)
  operator* (T1 const & a, T2 const & b) { 
   tail res(a); for (size_t i=0; i<res.size(); ++i) res(i)=res(i)*b; return res;
  }

 inline tail operator * (dcomplex a, tail_view const & r) { tail res(r); res.data_view()*=a; return res;}
 inline tail operator * (tail_view const & r, dcomplex a) { return a*r; }

 template<typename T1, typename T2> TYPE_ENABLE_IF(tail,mpl::and_<LocalTail<T1>, LocalTail<T2>>)
  operator/ (T1 const & a, T2 const & b) { return a *inverse(b); }

 inline tail operator / (tail_view const & r, dcomplex a) { tail res(r); res.data_view() /=a; return res;}
 inline tail operator / (dcomplex a, tail_view const & r) { return a * inverse(r); }

 template<typename T1, typename T2> TYPE_ENABLE_IF(tail,mpl::and_<LocalTail<T1>, LocalTail<T2>>)
  operator + (T1 const & l, T2 const& r) { 
   using arrays::range;
   if (l.shape() != r.shape()) TRIQS_RUNTIME_ERROR<< "tail addition : shape mismatch";
   tail res(l.shape(), std::min(l.order_min() ,r.order_min()), std::min(l.size(),r.size()));
   for (int i = res.order_min(); i<res.order_max(); ++i) res(i) = l(i) + r(i); 
   return res;
  }

 template<typename T1, typename T2> TYPE_ENABLE_IF(tail,mpl::and_<LocalTail<T1>, LocalTail<T2>>)
  operator - (T1 const & l, T2 const& r) { 
   using arrays::range;
   if (l.shape() != r.shape()) TRIQS_RUNTIME_ERROR<< "tail addition : shape mismatch";
   tail res(l.shape(), std::min(l.order_min() ,r.order_min()), std::min(l.size(),r.size()));
   for (int i = res.order_min(); i<res.order_max(); ++i) res(i) = l(i) - r(i); 
   return res;
  }

 template<typename T1, typename T2> TYPE_ENABLE_IF(tail,mpl::and_<is_scalar_or_element<T1>, LocalTail<T2>>)
  operator + (T1 const & a, T2 const & t) { 
   tail res(t.shape(), std::min(0 ,t.order_min()), t.size());
   for (int i = res.order_min(); i<res.order_max(); ++i) res(i) = t(i); 
   res(0) += a;
   return res;
  }

 template<typename T1, typename T2> TYPE_ENABLE_IF(tail,mpl::and_<LocalTail<T1>, is_scalar_or_element<T2>>)
  operator + (T1 const & t, T2 const & a) { return a+t;}

 template<typename T1, typename T2> TYPE_ENABLE_IF(tail,mpl::and_<is_scalar_or_element<T1>, LocalTail<T2>>)
  operator - (T1 const & a, T2 const & t) { return (-a) + t;}

 template<typename T1, typename T2> TYPE_ENABLE_IF(tail,mpl::and_<LocalTail<T1>, is_scalar_or_element<T2>>)
  operator - (T1 const & t, T2 const & a) { return (-a) + t;}

}}}
#endif
