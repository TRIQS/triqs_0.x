/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011-2013 by O. Parcollet
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
#ifndef TRIQS_ARRAYS_FLAGS_H
#define TRIQS_ARRAYS_FLAGS_H
#include "../indexmaps/permutation.hpp"
namespace triqs { namespace arrays {

 typedef unsigned long long ull_t;

 // Flags is a 64 bit unsigned int.
 // 0 is the default option.
 // Meaning of the bits : 
 // 0   -> Boundcheck
 // 1,2 -> Predefined order : 
 //    0 : None, 1 : C, 2 : Fortran 
 // 3,4 -> Init 
 //    0 : Noinit, 1 : NanInit, 2 : DefaultInit

 constexpr ull_t BoundCheck   = 1ull << 0;
 //constexpr ull_t COrder       = 1ull << 1;
 //constexpr ull_t FortranOrder = 1ull << 2;
 constexpr ull_t NanInit      = 1ull << 3;
 constexpr ull_t DefaultInit  = 1ull << 4;

 // NB : flags MUST be insensitive to slicing ...
 // i.e. when I slice, the flags does not change.
 
 namespace flags { 
  constexpr ull_t get(ull_t f, int a)   { return  (f & (1ull<<a)) >> a;}
  constexpr ull_t get2(ull_t f, int a)  { return  (f & (((1ull<<a) + (1ull<< (a+1))) >> a) );}

#ifdef TRIQS_ARRAYS_ENFORCE_BOUNDCHECK
  constexpr bool bound_check      (ull_t f) { return true;}
#else
  constexpr bool bound_check      (ull_t f) { return get (f, 0);}
#endif

  //constexpr bool c_order          (ull_t f) { return get (f, 1);}
  //constexpr bool fortran_order    (ull_t f) { return get (f, 2);}

  constexpr ull_t memory_order (int r, ull_t f, ull_t mo) { return mo + get(f,1) * permutations::identity(r) + get(f,2) * permutations::ridentity(r);}

  //or return an int 0,1,2
  //constexpr bool nan_init         (ull_t f) { return get (f, 3);}
  //constexpr bool default_init     (ull_t f) { return get (f, 4);}
  //constexpr bool no_init          (ull_t f) { return (!( default_init(f) || no_init(f))); } 

  constexpr ull_t init_mode        (ull_t f) { return get2 (f,3);}

  template<ull_t F> struct init_tag1;
  template<> struct init_tag1<0> { typedef Tag::no_init type;};
  template<> struct init_tag1<1> { typedef Tag::nan_inf_init type;};
  template<> struct init_tag1<2> { typedef Tag::default_init type;};

  // for the init_tag, we pass the *whole* option flag.
  template<ull_t F> struct init_tag : init_tag1 < init_mode(F)> {};

  /*template<ull_t F, ull_t mo> struct assert_make_sense {
   static_assert ( (!( c_order(F) && fortran_order(F))), "You asked C and Fortran order at the same time...");
   static_assert ( (!( (c_order(F) || fortran_order(F)) && mo )), "You asked C or Fortran order and gave a memory order ...");
   static_assert ( (init_mode (F) != 3), "You asked nan and default init at the same time...");
   //static_assert ( (!( nan_init(F) && default_init(F))), "You asked nan and default init at the same time...");
  };
  */
 }

}}//namespace triqs::arrays 
#endif


