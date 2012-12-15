/*******************************************************************************
 *
 * TRIQS: a Toolbox for Research in Interacting Quantum Systems
 *
 * Copyright (C) 2011 by O. Parcollet
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
#ifndef TRIQS_ARRAY_FLAGS_H
#define TRIQS_ARRAY_FLAGS_H
namespace triqs { namespace arrays {

 typedef unsigned long long ull_t;

 // Flags is a 64 bit unsigned int.
 // Meaning of the bits : 
 // 0   -> Boundcheck
 // 1,2 -> Predefined order : 
 //    0 : None, 1 : C, 2 : Fortran 
 // 3,4 -> Init 
 //    0 : Noinit, 1 : NanInit, 2 : DefaultInit

 constexpr ull_t BoundCheck   = 1ull << 0;
 constexpr ull_t COrder       = 1ull << 1;
 constexpr ull_t FortranOrder = 1ull << 2;
 constexpr ull_t NanInit      = 1ull << 3;
 constexpr ull_t DefaultInit  = 1ull << 4;

 namespace flags { 
  constexpr ull_t get(ull_t f, int a)   { return  (f & (1ull<<a)) >> a;}
  constexpr ull_t get2(ull_t f, int a)  { return  (f & (1ull<<a + 1ull<< (a+1)) >> a );}

  constexpr bool bound_check      (ull_t f) { return get (f, 0);}

  constexpr bool c_order          (ull_t f) { return get (f, 1);}
  constexpr bool fortran_order    (ull_t f) { return get (f, 2);}

  //or return an int 0,1,2
  //constexpr bool nan_init         (ull_t f) { return get (f, 3);}
  //constexpr bool default_init     (ull_t f) { return get (f, 4);}
  //constexpr bool no_init          (ull_t f) { return (!( default_init(f) || no_init(f))); } 

  constexpr ull_t init_mode        (ull_t f) { return get2 (f,3);}

  template<ull_t F> struct assert_make_sense {
   static_assert ( (!( c_order(F) && fortran_order(F))), "You asked C and Fortran order at the same time...");
   static_assert ( (init_mode (F) != 3), "You asked nan and default init at the same time...");
   //static_assert ( (!( nan_init(F) && default_init(F))), "You asked nan and default init at the same time...");
  };
 }

 // move this out ....
 template < class V, int R, class Opt, class ViewTag > struct ViewFactory;

}}//namespace triqs::arrays 
#endif


