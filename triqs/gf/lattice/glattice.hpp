#ifndef glattice_hpp
#define glattice_hpp

#include "boost/fusion/include/vector.hpp"
#include "boost/mpl/bool.hpp"
#include "boost/mpl/if.hpp"
#include "../local/gf.hpp"

namespace triqs { namespace gf { namespace lattice {

 template<typename MeshType, typename GfLocal> class glattice;
 template<typename MeshType, typename GfLocal> class glattice_view;
 template<typename E> class glattice_expr;

 template<typename G> struct LatticeGf : boost::mpl::false_{};  // a boolean trait to identify the objects modelling the concept LatticeGf
 template<typename M,typename Gloc> struct LatticeGf<glattice<M,Gloc> >     : boost::mpl::true_{};
 template<typename M,typename Gloc> struct LatticeGf<glattice_view<M,Gloc> >: boost::mpl::true_{};
 template<typename E> struct LatticeGf<glattice_expr<E> >: boost::mpl::true_{};



 template<typename MeshType,typename GfLocal, bool IsView> class glattice_impl{

  public:
   typedef MeshType space_mesh_type;
   typedef typename MeshType::domain_type space_domain_type;
   typedef typename GfLocal::mesh_type time_mesh_type;
   typedef typename GfLocal::domain_type time_domain_type;

   typedef typename boost::fusion::vector< space_domain_type, time_domain_type > domain_type;
   typedef typename boost::fusion::vector< space_mesh_type, time_mesh_type > mesh_type;
   

   //typedef typename bf::result_of::at<domain_type,mpl::_int_<0> >::type space_domain_type2;
   //space_domain_type space_domain() const { return bf::at<mpl::int_<0> >(domain());}

   typedef arrays::vector<GfLocal> data_view_type;
   typedef arrays::vector_view<GfLocal> data_non_view_type;
   typedef typename boost::mpl::if_c<IsView, data_view_type,data_non_view_type>::type data_type;


   typedef glattice_view<MeshType,GfLocal> view_type;
   typedef glattice<MeshType,GfLocal> non_view_type;

   space_mesh_type const & space_mesh() const {return my_space_mesh;}



  protected:
   space_mesh_type my_space_mesh;
   data_type data;

   //constructors
   glattice_impl();
   glattice_impl(space_mesh_type const & m, data_view_type const & dat): my_space_mesh(m),data(dat){}
   glattice_impl(glattice_impl const & x) : my_space_mesh(x.my_space_mesh),data(x.data){}


  public:
   template<typename RHS> void operator = (RHS const & rhs){
    my_space_mesh = rhs.space_mesh();
    //BOOST_AUTO(r0,rhs(space_mesh()[0]));
    //tqa::resize_or_check_if_view(data,tqa::make_shape(r0.shape()[0],r0.shape[1],my_space_mesh.size()));

   }
   TRIQS_NVL_HAS_AUTO_ASSIGN();

   //call operators
   typedef typename GfLocal::mv_type mv_type;
   typedef typename GfLocal::const_mv_type const_mv_type;


   //operator () (bf::vector<space_domain_type::element_type, time_domain_type::element_type> const & ) {}
   // appel G( bf::make_vector( k, omega) 

   mv_type operator() (meshes::mesh_pt<space_mesh_type> const & x,meshes::mesh_pt<time_mesh_type> const & y) {
    return (data(x.i))(arrays::range(),arrays::range(),y.i);
   }
   const_mv_type operator() (meshes::mesh_pt<space_mesh_type> const & x,meshes::mesh_pt<time_mesh_type> const & y) const {
    return (data(x.i))(arrays::range(),arrays::range(),y.i);
   }


   typedef typename space_domain_type::point_type arg0_type;
   typedef typename time_domain_type::point_type arg1_type;

   mv_type operator()(arg0_type const & i,arg1_type const & j)  {

    
    return NULL;//(this->space_mesh().interpolate(*this,i))->mesh().interpolate(*this.j) ;
   }
   const_mv_type operator()(arg0_type const & i,arg1_type const & j) const {
    return NULL;//(this->space_mesh().interpolate(*this,i))->mesh().interpolate(*this.j) ;
   }

   TRIQS_LAZY_ADD_LAZY_CALL_WITH_VIEW(2,view_type);

   void save(std::string file,bool accumulate = false) const{}
   void load(std::string file){}


   friend std::ostream & operator << (std::ostream & out,glattice_impl const & x) { return out << "glattice_view";}

 };


 template<typename MeshType,typename GfLocal> class glattice : public glattice_impl<MeshType,GfLocal,false>{
  typedef glattice_impl<MeshType,GfLocal,false> B;
  public:
  glattice():B(){}
  glattice(glattice const & g) : B(g){}
  glattice(glattice_view<MeshType,GfLocal> const & g) : B(g){}
  template<typename GfType> glattice(GfType const & x): B() { *this = x;}


  using B::operator=; 

 };


 template<typename MeshType,typename GfLocal> class glattice_view : public glattice_impl<MeshType,GfLocal,true>{
  typedef glattice_impl<MeshType,GfLocal,true> B;
  protected:
  friend class glattice_impl<MeshType, GfLocal, true>; friend class glattice_impl<MeshType, GfLocal, false>;
  public:
  glattice_view(glattice_view const & g):B(g){}


 };

 // -------------------------------   Expression template for gf  --------------------------------------------------

 template< typename T> struct scalar_identification_trait : local::is_scalar_or_element<T> {}; // is_scalar_or_element defined in tail
 template< typename T> struct gf_identification_trait : LatticeGf<T> {};
 static const int gf_arity=2;
 // this include in intentionnally IN the namespace
#include "../gf_proto.hpp"

}}}
#endif

