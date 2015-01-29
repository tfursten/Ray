#pragma once
#ifndef REXP_H
#define REXP_H

#include <cstdint>
#include <cmath>
#include <cfloat>

#include "xorshift64.h"

extern const double ew[256];
extern const double ef[256];
extern const int64_t ek[256];

inline double rand_exp_inv(xorshift64 &rng) { return -log(rng.get_double52()); }

inline double rand_exp_zig(xorshift64 &rng) {
	const double r = 7.69711747013104972;
	uint64_t u = rng.get_uint64();
	// use the top 8 high bits for b
	uint64_t b = u >> 56;
	// use the rest for a
	int64_t a = static_cast<int64_t>(u & UINT64_C(0x00ffffffffffffff)); 
	while( a > ek[b] ) {
		if(b == 0)
			return r+rand_exp_inv(rng);
		double x = a*ew[b];
		// we can cache ef[b-1]-ef[b], but it should be minor
		if(ef[b]+(ef[b-1]-ef[b])*rng.get_double52() < exp(-x) )
			return x;
		u = rng.get_uint64();
		b = u >> 56;
		a = static_cast<int64_t>(u & UINT64_C(0x00ffffffffffffff)); 
	}
	
	return a*ew[b];
}

inline double rand_exp(xorshift64 &rng, double rate = 1.0) { return rand_exp_zig(rng)/rate; }

inline double rand_exp_mean(xorshift64 &rng, double mu = 1.0) { return rand_exp_zig(rng)*mu; }

inline double rand_exp_trunc(xorshift64 &rng, double lim, double rate = 1.0) {
	double u = rng.get_double52();
	return -log1p(u*(expm1(-rate*lim)))/rate;
}

#endif
