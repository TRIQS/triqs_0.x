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
#include <triqs/utility/apply_on_tuple.hpp>
//#include <functional>
#include <iostream>
#include <cmath>

//double f(int i, double x, double y, int k) { std::cout   << i << x << y <<k <<std::endl ; return k + i+ x + 2*y;}
struct fun { 
 double operator()(int i, double x, double y, int k) {  return 6*k + i - 1.3*x + 2*y;}
};

int main(int argc, char **argv) {
 //std::function<double(int,double,double,int)> F(f);
 auto t = std::make_tuple(1,2.3,4.3,8);
 fun F;
 std::cerr  << " f(t) =" << triqs::apply_on_tuple(F,t)<< std::endl ;
 if ( std::abs((triqs::apply_on_tuple(F,t) -  F(1,2.3,4.3,8))) > 1.e-13) throw std::runtime_error(" ");
}




