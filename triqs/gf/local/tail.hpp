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
#include <triqs/arrays.hpp>
#include <triqs/arrays/algorithms.hpp>
#include <triqs/gf/tools.hpp>
#include <triqs/python_tools/cython_proxy.hpp>

namespace triqs { namespace gf { namespace local {

 namespace details {
   static double small = 1.e-10;
 }

 namespace tqa= triqs::arrays; namespace tql= triqs::clef;  namespace mpl= boost::mpl;
 typedef std::complex<double> dcomplex;

 class tail;                            // the value class
 class tail_view;                       // the view class

 template<typename G> struct LocalTail  : mpl::false_{};  // a boolean trait to identify the objects modelling the concept LocalTail
 template<> struct LocalTail<tail >     : mpl::true_{};
 template<> struct LocalTail<tail_view >: mpl::true_{};
 template<> struct LocalTail<python_tools::cython_proxy<tail_view>>: mpl::true_{};

 // a trait to find the scalar of the algebra i.e. the true scalar and the matrix ...
 template <typename T> struct is_scalar_or_element : mpl::or_< tqa::ImmutableMatrix<T>, utility::is_in_ZRC<T> > {};

 // ---------------------- implementation --------------------------------

 /// A common implementation class. Idiom : ValueView
 template<bool IsView> class tail_impl  {
  public :
   typedef void has_view_type_tag; // Idiom : ValueView
   typedef tail_view view_type;
   typedef tail      non_view_type;

   typedef arrays::array      <dcomplex,3>                         data_non_view_type;
   typedef arrays::array_view <dcomplex,3>                         data_view_type;
   typedef typename mpl::if_c<IsView, data_view_type, data_non_view_type>::type  data_type;

   typedef arrays::array<long,2>                                   mask_non_view_type;
   typedef arrays::array_view<long,2>                              mask_view_type;
   typedef typename mpl::if_c<IsView, mask_view_type, mask_non_view_type>::type  mask_type;

   typedef arrays::matrix_view<dcomplex>  mv_type;
   typedef arrays::matrix_view<dcomplex>  const_mv_type;
   //typedef arrays::matrix_view<const dcomplex>  const_mv_type;

   data_view_type data_view()             { return data;}
   const data_view_type data_view() const { return data;}
   mask_view_type mask_view()             { return mask;}
   const mask_view_type mask_view() const { return mask;}

   long order_min() const {return omin;}
   long order_max() const {return min_element(mask);}
   size_t size() const {return data.shape()[0];}
   long smallest_nonzero() const {
     long om = omin;
     while ((om < this->order_max()) && (max_element(abs(data(om-omin,tqa::range(),tqa::range()))) < details::small)) om++;
     return om;
   }

   typedef tqa::mini_vector<size_t,2> shape_type;
   shape_type shape() const { return shape_type(data.shape()[1], data.shape()[2]);}
   size_t shape(int i) const { return data.shape()[i];}

   bool is_decreasing_at_infinity() const { return (smallest_nonzero() >=1);}

  protected:

   long omin;
   mask_type mask;
   data_type data;

   // All constructors
   tail_impl(): omin(0), mask(), data() {} // all arrays of zero size (empty)
   tail_impl(size_t N1, size_t N2, size_t size_, long order_min):
     omin(order_min), mask(tqa::make_shape(N1,N2)), data(tqa::make_shape(size_,N1,N2)) {
     mask() = order_min+size_-1;
     data() = 0;
   }
   tail_impl(data_type const &d, long order_min, mask_type const &om): omin(order_min), mask(om), data(d) {}
 //  tail_impl(tail_impl          const & x): omin(x.omin), mask(x.mask), data(x.data){}
   tail_impl(tail_impl<!IsView> const & x): omin(x.omin), mask(x.mask), data(x.data){}
 tail_impl(tail_impl const &) = default;
  tail_impl(tail_impl &&) = default;
 
   friend class tail_impl<!IsView>;
  public:

   mv_type operator() (int n) {
     if (n>this->order_max()) TRIQS_RUNTIME_ERROR<<" n > Max Order. n= "<<n <<", Max Order = "<<order_max() ;
     if (n<this->order_min()) TRIQS_RUNTIME_ERROR<<" n < Min Order. n= "<<n <<", Min Order = "<<order_min() ;
     return this->data(n-omin, tqa::range(), tqa::range());
   }

   const_mv_type operator() (int n) const {
     if (n>this->order_max()) TRIQS_RUNTIME_ERROR<<" n > Max Order. n= "<<n <<", Max Order = "<<order_max() ;
     if (n<this->order_min())  { mv_type::non_view_type r(this->shape()); r()=0; return r;}
     return this->data(n-omin,tqa::range(), tqa::range());
   }

   operator freq_infty() const { return freq_infty();}

   /// Save in txt file : doc the format  ? ---> prefer serialization or hdf5 !
   void save(std::string file,  bool accumulate=false) const {}

   /// Load from txt file : doc the format ?
   void load(std::string file){}

   ///
   friend void h5_write (h5::group fg, std::string subgroup_name, tail_impl const & t) {
    auto  gr = fg.create_group(subgroup_name);
    // Add the attribute
    h5_write(gr,"omin",t.omin);
    h5_write(gr,"mask",t.mask);
    h5_write(gr,"data",t.data);
   }

   friend void h5_read  (h5::group fg, std::string subgroup_name, tail_impl & t){
    auto gr = fg.open_group(subgroup_name);
    // Check the attribute or throw
    h5_read(gr,"omin",t.omin);
    h5_read(gr,"mask",t.mask);
    h5_read(gr,"data",t.data);
   }

   //  BOOST Serialization
   friend class boost::serialization::access;
   template<class Archive>
    void serialize(Archive & ar, const unsigned int version) {
     ar & boost::serialization::make_nvp("omin",omin);
     ar & boost::serialization::make_nvp("mask",mask);
     ar & boost::serialization::make_nvp("data",data);
    }

   friend std::ostream & operator << (std::ostream & out, tail_impl const & x) {
    out <<"tail/tail_view: min/smallest/max = "<< x.order_min() << " " << x.smallest_nonzero() << " "<< x.order_max();
    for (long u = x.order_min(); u <= x.order_max(); ++u) out <<"\n ...  Order "<<u << " = " << x(u);
    return out;
   }

 };

 // -----------------------------
 ///The View class of GF
 class tail_view : public tail_impl <true> {
  typedef tail_impl <true>  B;
  friend class tail;

  public :
  template<bool V> tail_view(tail_impl<V> const & t): B(t){}
  tail_view(B::data_type const &d, long order_min, B::mask_type const &om): B(d, order_min, om){}
  tail_view(tail_view const &) = default;
  tail_view(tail_view &&) = default;
  void rebind( tail_view const &X) {
    omin = X.omin;
    mask.rebind(X.mask);
    data.rebind(X.data);
  }
  inline void rebind( tail const &X);

  //using B::operator=; // import operator = from impl. class or the default = is synthetized
  //tail_view & operator=(const tail_view & rhs) { if (this->data.is_empty()) rebind(rhs); else B::operator=(rhs); return *this; }

  // operator = for views
  tail_view & operator = (const tail_view & rhs) {
    if (this->data.is_empty()) rebind(rhs);
    else {
      if (rhs.omin < omin) TRIQS_RUNTIME_ERROR<<"rhs has too small omin";
      if ((data.shape()[1] != rhs.data.shape()[1]) || (data.shape()[2] != rhs.data.shape()[2]))
        TRIQS_RUNTIME_ERROR<<"rhs has different shape";
      for (size_t i=0; i<mask.shape()[0]; ++i)
        for (size_t j=0; j<mask.shape()[1]; ++j)
          mask(i,j) = std::min(rhs.mask(i,j), long(omin+size()-1));
      for (size_t n=0; n<std::min(size(), size_t(rhs.size()-omin+rhs.omin)); ++n)
        if (n < rhs.omin-omin) data(n,tqa::range(),tqa::range()) = 0.0;
        else data(n,tqa::range(),tqa::range()) = rhs.data(n-rhs.omin+omin,tqa::range(),tqa::range());
    }
    return *this;
  }
  inline tail_view & operator=(const tail & rhs);

  tail_view & operator=(std::complex<double> const & x) {
    if (omin > 0) TRIQS_RUNTIME_ERROR<<"lhs has too large omin";
    for (size_t n=0; n<size(); ++n) data(n, tqa::range(), tqa::range()) = 0.0;
    data(-omin, tqa::range(), tqa::range()) = x;
    mask() = omin+size()-1;
    return *this;
  }

  using B::operator(); // import all previously defined operator() for overloading
  friend std::ostream & triqs_nvl_formal_print(std::ostream & out, tail_view const & x) { return out<<"tail_view";}
  
  void print_me() const { std::cout  << *this << std::endl ; }
 };

 // -----------------------------
 ///The regular class
 class tail : public tail_impl <false> {
  typedef tail_impl <false>  B;
  friend class tail_view;
  public :
  tail():B() {}
  typedef tqa::mini_vector<size_t,2> shape_type;
  tail(size_t N1, size_t N2, size_t size_ = 10, long order_min=-1): B(N1,N2,size_,order_min) {}
  tail(shape_type const & sh, size_t size_ = 10, long order_min=-1): B(sh[0],sh[1],size_,order_min) {}
  tail(tail const & g): B(g){}
  tail(tail_view const & g): B(g){}
  tail(tail &&) = default;

  // operator = for values
  tail & operator = (tail_view const & rhs) {
    omin = rhs.omin;
    mask = rhs.mask;
    data = rhs.data;
    return *this;
  }
  tail & operator = (tail const & rhs) {
    omin = rhs.omin;
    mask = rhs.mask;
    data = rhs.data;
    return *this;
  }

  using B::operator();

  /// The simplest tail corresponding to  : omega
  static tail_view omega(size_t N1, size_t N2, size_t size_) {
    tail t(N1, N2, size_, -1);
    t(-1) = 1;
    return t;
  }

  /// The simplest tail corresponding to  : omega, constructed from a shape for convenience
  static tail_view omega(tail::shape_type const & sh, size_t size_) { return omega(sh[0],sh[1],size_);}

 };

  template<typename RHS> void assign_from_expression(tail_view & t,RHS const & rhs) { t = rhs( tail::omega(t.shape(),t.size())); }
 
  inline void tail_view::rebind( tail const &X) {
   omin = X.omin;
   mask.rebind(X.mask);
   data.rebind(X.data);
 }
 inline tail_view & tail_view::operator = (const tail & rhs) {
    if (this->data.is_empty()) rebind(rhs);
    else {
      if (rhs.omin < omin) TRIQS_RUNTIME_ERROR<<"rhs has too small omin";
      if ((data.shape()[1] != rhs.data.shape()[1]) || (data.shape()[2] != rhs.data.shape()[2]))
        TRIQS_RUNTIME_ERROR<<"rhs has different shape";
      for (size_t i=0; i<mask.shape()[0]; ++i)
        for (size_t j=0; j<mask.shape()[1]; ++j)
          mask(i,j) = std::min(rhs.mask(i,j), long(omin+size()-1));
      for (size_t n=0; n<std::min(size(),size_t(rhs.size()-omin+rhs.omin)); ++n)
        if (n < rhs.omin-omin) data(n,tqa::range(),tqa::range()) = 0.0;
        else data(n,tqa::range(),tqa::range()) = rhs.data(n-rhs.omin+omin, tqa::range(), tqa::range());
    }
    return *this;
  }

 /// Slice in orbital space
 template<bool V> tail_view slice_target(tail_impl<V> const & t, tqa::range R1, tqa::range R2) {
  return tail_view(t.data_view()(tqa::range(),R1,R2),t.order_min(),t.mask_view()(R1,R2));
 }

 inline tail inverse(tail_view const & t) {

  // find in t
  long omin1 = - t.smallest_nonzero();
  long omax1 = t.order_max() + 2*omin1;
  size_t si = omax1-omin1+1;
  tail t_inv(t.shape(), si, omin1);

  t_inv(omin1) = inverse(t(-omin1));

  for (size_t n=1; n<si; n++) {
   for (size_t p=0; p<n; p++) {
     t_inv(omin1 + n) -= t(n-omin1-p) * t_inv(omin1+p);
   }
   t_inv(omin1 + n) = t_inv(omin1) * make_clone(t_inv(omin1 + n));
  }
  return t_inv;
 }

 inline tail mult_impl(tail_view const & l, tail_view const& r) {
  if (l.shape()[1] != r.shape()[0]) TRIQS_RUNTIME_ERROR<< "tail multiplication : shape mismatch";
  long omin1 = l.smallest_nonzero() + r.smallest_nonzero();
  long omax1 = std::min(r.order_max()+l.smallest_nonzero(), l.order_max()+r.smallest_nonzero());
  size_t si = omax1-omin1+1;
  tail res(l.shape()[0], r.shape()[1], si, omin1);
  res.data_view() =0;
  for (long n=res.order_min(); n<=res.order_max(); ++n) {
   // sum_{p}^n a_p b_{n-p}. p <= a.n_max, p >= a.n_min and n-p <=b.n_max and n-p >= b.n_min
   // hence p <= min ( a.n_max, n-b.n_min ) and p >= max ( a.n_min, n- b.n_max)
   const long pmin = std::max(l.smallest_nonzero(), n - r.order_max() );
   const long pmax = std::min(l.order_max(), n - r.smallest_nonzero() );
   for (long p = pmin; p <= pmax; ++p)  { res(n) += l(p) * r(n-p);}
  }
  return res;
 }

 template<typename T1, typename T2>
  TYPE_ENABLE_IF(tail,mpl::and_<LocalTail<T1>, LocalTail<T2>>)
  operator* (T1 const & a, T2 const & b) { return mult_impl(a,b); }

 template<typename T1, typename T2> TYPE_ENABLE_IF(tail,mpl::and_<tqa::ImmutableMatrix<T1>, LocalTail<T2>>)
  operator* (T1 const & a, T2 const & b) {
   tail res(b); for (long n=res.order_min(); n<=res.order_max(); ++n) res(n)=a*res(n); return res;
  }

 template<typename T1, typename T2> TYPE_ENABLE_IF(tail,mpl::and_<LocalTail<T1>, tqa::ImmutableMatrix<T2>>)
  operator* (T1 const & a, T2 const & b) {
   tail res(a); for (long n=res.order_min(); n<=res.order_max(); ++n) res(n)=res(n)*b; return res;
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
   long omin1 = std::min(l.smallest_nonzero(),r.smallest_nonzero());
   long omax1 = std::min(l.order_max(),r.order_max());
   size_t si = omax1-omin1+1;
   tail res(l.shape(), si, omin1);
   for (long i = res.order_min(); i<=res.order_max(); ++i) res(i) = l(i) + r(i);
   return res;
  }

 template<typename T1, typename T2> TYPE_ENABLE_IF(tail,mpl::and_<LocalTail<T1>, LocalTail<T2>>)
  operator - (T1 const & l, T2 const& r) {
   using arrays::range;
   if (l.shape() != r.shape()) TRIQS_RUNTIME_ERROR<< "tail substraction : shape mismatch";
   long omin1 = std::min(l.smallest_nonzero(),r.smallest_nonzero());
   long omax1 = std::min(l.order_max(),r.order_max());
   size_t si = omax1-omin1+1;
   tail res(l.shape(), si, omin1);
   for (long i = res.order_min(); i<=res.order_max(); ++i) res(i) = l(i) - r(i);
   return res;
  }

 template<typename T1, typename T2> TYPE_ENABLE_IF(tail,mpl::and_<is_scalar_or_element<T1>, LocalTail<T2>>)
  operator + (T1 const & a, T2 const & t) {
   tail res(t.shape(), t.size(), std::min(long(0), t.smallest_nonzero()));
   for (long i = res.order_min(); i<=res.order_max(); ++i) res(i) = t(i);
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
