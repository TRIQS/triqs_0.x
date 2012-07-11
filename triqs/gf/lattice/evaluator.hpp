#ifndef EVALUATOR_HPP
#define EVALUATOR_HPP

template<typename F, typename TAG=void> class evaluator;

template<typename F, typename TAG, typename SpaceMeshType,typename TimeMeshType, bool IsView>  
class evaluator<glattice_impl<SpaceMeshType,TimeMeshType,IsView>, void> {
 F const & f;
 public:
 evaluator(F const & f_):f(f_){}
 //typename F::result_type operator() (typename F::arg0_type const & i, typename F::arg1_type const & j) const {  return f(i,j); }

 template<typename Ktype, typename OmegaType>
 auto operator() (Ktype const & k, OmegaType const & om) const ->decltype(f(k,om))  {  return f(k,om); }

};

template<typename F, typename TAG> class evaluator<glattice_impl, interpolator_tag> {
 F const & f;
 public:
 evaluator(F const & f_):f(f_){}
 typename F::result_type operator() (typename F::arg0_type const & i, typename F::arg1_type const & j) const {  return 2*f(i,j); }
};


#endif
