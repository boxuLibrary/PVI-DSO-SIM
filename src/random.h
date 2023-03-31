/////////////////////////////////////////////////////////////////////////////////
//
//  Multilayer Feature Graph (MFG), version 1.0
//  Copyright (C) 2011-2015 Yan Lu, Madison Treat, Dezhen Song
//  Netbot Laboratory, Texas A&M University, USA
//
//  This program is free software; you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation; either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//
/////////////////////////////////////////////////////////////////////////////////

/********************************************************************************
 * Random number generator
 ********************************************************************************/

/*
 * Our implementation of 2 Psuedo-Random Number Generators, both variants of the
 * XOR Shift PRNG algorithm.  Both variants here were copied directly from
 * Wikipedia (http://en.wikipedia.org/wiki/Xorshift), which copied them directly
 * from "An experimental exploration of Marsaglia's xorshift generators,
 * scrambled" by Sebastiano Vigna (2014), located online at
 * http://arxiv.org/abs/1402.6246
 */

#ifndef RANDOM_H
#define RANDOM_H

#ifdef	__cplusplus
extern "C" {
#endif

#ifdef _MSC_VER
typedef unsigned __int64 uint64_t;
#define UINT64_C(val) (val##ui64)
#else // _MSC_VER
#include <stdint.h>
// if UINT64_C is already defined in stdint.h, there is no need to redefine it
#ifndef UINT64_C
#define UINT64_C(val) (val##ULL)
#endif // UINT64_C

#endif // _MSC_VER

/*
 * This 64-bit XOR Shift* algorithm is used to provide a single 64-bit seed for
 * the 1024-bit XOR Shift* algorithm, which passes more tests than this variant.
 */
uint64_t xorshift64star(void);

/*
 * This 1024-bit XOR Shift* algorithm is our implementation of a PRNG to provide
 * repeatable results from MFG calculations.
 */
uint64_t xorshift1024star(void);

/*
 * A seed function and simple function wrapper for xorshift1024star().
 */
void seed_xrand(uint64_t val);
uint64_t xrand(void);

#ifdef	__cplusplus
}
#endif

#endif	/* RANDOM_H */

