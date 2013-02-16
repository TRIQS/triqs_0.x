#include "./common.hpp"

using std::cout; using  std::endl;

struct myexpr {

 TRIQS_CLEF_IS_EXPRESSION();
 typedef bf::set<> ph_set; 
 template<typename BfMapType> struct call_rtype { typedef double type; };

 template<typename BfMapType> typename call_rtype<BfMapType>::type  
  eval_on_map (BfMapType const & ph_value_dict) const { return 78;} 
 
 friend std::ostream & triqs_nvl_formal_print(std::ostream & out, myexpr const & x) { return out<<"myexpr";}
};

int main() { 
 myexpr e1;
 TEST(y_ + e1);
 TEST(tql::eval(e1 ,y_=3));
}

