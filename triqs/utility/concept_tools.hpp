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
#include <boost/preprocessor/seq/enum.hpp>

//#define TRIQS_CONCEPT_TAG(MyBeautifulConcept) __TRIQS_models_concept_##MyBeautifulConcept##_tag
#define TRIQS_MODEL_CONCEPT(MyBeautifulConcept) concept_tags::MyBeautifulConcept

#define TRIQS_DEFINE_CONCEPT_AND_ASSOCIATED_TRAIT(MyBeautifulConcept) \
 namespace concept_tags { struct MyBeautifulConcept{};}\
 template<typename T> struct MyBeautifulConcept : boost::is_base_of<concept_tags::MyBeautifulConcept, T> {};

#define TRIQS_DEFINE_CONCEPT_AND_ASSOCIATED_TRAIT_R(MyBeautifulConcept,Rs) \
 namespace concept_tags { struct MyBeautifulConcept: BOOST_PP_SEQ_ENUM(Rs) {};}\
 template<typename T> struct MyBeautifulConcept : boost::is_base_of<concept_tags::MyBeautifulConcept, T> {};

//#define TRIQS_DEFINE_CONCEPT_ASSOCIATED_TRAIT1(MyBeautifulConcept) \
// template<typename T, typename Enable=void> struct MyBeautifulConcept : mpl::false_{};\
// template<typename T> struct MyBeautifulConcept<T, typename T::TRIQS_CONCEPT_TAG(MyBeautifulConcept)> : mpl::true_{};

//#define TRIQS_REFINE_CONCEPT(C1, C2) template<typename T> struct C2<T, typename T::TRIQS_CONCEPT_TAG(C1)> : mpl::true_{};

#ifdef TRIQS_COMPILE_TIME_DEBUG
#define TRIQS_ASSERT_MODEL_CONCEPT(MyBeautifulConcept,T)  BOOST_CONCEPT_ASSERT((BCC_##MyBeautifulConcept<T>));
#else
#define TRIQS_ASSERT_MODEL_CONCEPT(MyBeautifulConcept,T)
#endif

namespace triqs { namespace utility { 
}}
#endif

