/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2012 by O. Parcollet
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
#ifndef TRIQS_CLEF_MATH_H 
#define TRIQS_CLEF_MATH_H
#include <math.h>
#include "./core.hpp"

namespace triqs { namespace clef { 

#define TRIQS_CLEF_MATH_FNT1 (cos)(sin)(tan)(cosh)(sinh)(tanh)(acos)(asin)(atan)(exp)(log)(sqrt)(abs)(floor)
#define TRIQS_CLEF_MATH_FNT2 (pow)

#define AUX1(r, data, elem) using std::elem;
#define AUX(r, data, elem) TRIQS_CLEF_MAKE_FNT_LAZY(1,elem)
 BOOST_PP_SEQ_FOR_EACH(AUX1, nil , TRIQS_CLEF_MATH_FNT1);
 BOOST_PP_SEQ_FOR_EACH(AUX , nil , TRIQS_CLEF_MATH_FNT1);
#undef AUX  

#define AUX(r, data, elem) TRIQS_CLEF_MAKE_FNT_LAZY(2,elem)
 BOOST_PP_SEQ_FOR_EACH(AUX1, nil , TRIQS_CLEF_MATH_FNT2);
 BOOST_PP_SEQ_FOR_EACH(AUX, nil , TRIQS_CLEF_MATH_FNT2);
#undef AUX  
#undef AUX1  

#undef TRIQS_CLEF_MATH_FNT1
#undef TRIQS_CLEF_MATH_FNT2

}}

#endif

