/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2013 by O. Parcollet
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
#include <triqs/utility/tuple_tools.hpp>
//#include <functional>
#include <iostream>
#include <cmath>
#include <stdexcept>

//double f(int i, double x, double y, int k) { std::cout   << i << x << y <<k <<std::endl ; return k + i+ x + 2*y;}
struct fun { 
 double operator()(int i, double x, double y, int k) {  return 6*k + i - 1.3*x + 2*y;}
};

struct fun2 { 
 double operator()(double x, double y) {  return x+y;}
};

int main(int argc, char **argv) {

 auto t = std::make_tuple(1,2.3,4.3,8);
 auto t2 = std::make_tuple(1,2,3,4);
 auto t1 = std::make_tuple(1,2.3,4.3,8);
  
 {
  auto res = triqs::tuple::apply(fun(),t);
  std::cerr  << " f(t) =" << res << std::endl ;
  if ( std::abs((res -  fun()(1,2.3,4.3,8))) > 1.e-13) throw std::runtime_error(" ");
 }

 {
 auto r = triqs::tuple::apply_on_zip(fun2(),t1,t2);
  std::cerr  << " [f(a,b) for (a,b) in zip(t1,t2)] =" 
   << std::get<0>(r) << " "
   << std::get<1>(r) << " "
   << std::get<2>(r) << " "
   << std::get<3>(r) << std::endl;
 }

 {
  auto res = triqs::tuple::fold([](double x,double r) { return x+r;}, t, 0);
  std::cerr  << " " << res << std::endl ;
  if ( std::abs((res -  15.6)) > 1.e-13) throw std::runtime_error(" ");
 }

 {
  auto res = triqs::tuple::fold_on_zip([](double x, double y, double r) { return x+ 2*y +r;}, t1,t2, 0);
  std::cerr  << " " << res << std::endl ;
  if ( std::abs((res -  35.6)) > 1.e-13) throw std::runtime_error(" ");
 }


}




