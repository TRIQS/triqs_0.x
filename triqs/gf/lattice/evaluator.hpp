#ifndef EVALUATOR_HPP
#define EVALUATOR_HPP

template<typename Tag, typename F> class evaluator{
  F const & f;
  public:
    evaluator(F const & f_):f(f_){}

    F::result_type operator() (F::arg0_type const i,F::arg1_type const j) const {  return f(i,j); }
    F::result_type operator<Default,glattice_impl>() (F::arg0_type const i,F::arg1_type const j) const {  return f(i,j); }
}


#endif
