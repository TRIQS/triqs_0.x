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
#ifndef TRIQS_UTILITY_CONCEPT_TOOLS_H
#define TRIQS_UTILITY_CONCEPT_TOOLS_H
#include <boost/preprocessor/punctuation/comma_if.hpp>
#include <boost/preprocessor/seq/for_each_i.hpp>
#include <boost/type_traits/is_base_of.hpp>
#include <boost/type_traits/is_complex.hpp>

#define TRIQS_CONCEPT_TAGNAME(MyBeautifulConcept) MyBeautifulConcept##__concept_tag

#define TRIQS_MODEL_CONCEPT(MyBeautifulConcept) TRIQS_CONCEPT_TAGNAME(MyBeautifulConcept)

#define TRIQS_DEFINE_CONCEPT_AND_ASSOCIATED_TRAIT(MyBeautifulConcept) \
 struct TRIQS_CONCEPT_TAGNAME(MyBeautifulConcept) {};\
 template<typename T> struct MyBeautifulConcept : boost::is_base_of<TRIQS_CONCEPT_TAGNAME(MyBeautifulConcept) , T> {};

#define TRIQS_DEFINE_CONCEPT_AND_ASSOCIATED_TRAIT_R_AUX(r, data, i, elem) BOOST_PP_COMMA_IF(i) TRIQS_CONCEPT_TAGNAME(elem) 

#define TRIQS_DEFINE_CONCEPT_AND_ASSOCIATED_TRAIT_R(MyBeautifulConcept,Rs) \
 struct TRIQS_CONCEPT_TAGNAME(MyBeautifulConcept) : BOOST_PP_SEQ_FOR_EACH_I (TRIQS_DEFINE_CONCEPT_AND_ASSOCIATED_TRAIT_R_AUX,nil,Rs) {};\
 template<typename T> struct MyBeautifulConcept : boost::is_base_of<TRIQS_CONCEPT_TAGNAME(MyBeautifulConcept), T> {};

#ifdef TRIQS_COMPILE_TIME_DEBUG
#define TRIQS_ASSERT_MODEL_CONCEPT(MyBeautifulConcept,T)  BOOST_CONCEPT_ASSERT((BCC_##MyBeautifulConcept<T>));
#else
#define TRIQS_ASSERT_MODEL_CONCEPT(MyBeautifulConcept,T)
#endif

namespace triqs { namespace utility { 
}}
#endif

