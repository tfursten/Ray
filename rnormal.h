#pragma once
#ifndef RNORMAL_H
#define RNORMAL_H

#include <cstdint>
#include <cmath>

#include "xorshift64.h"
#include "rexp.h"

/*
 *     Ziggurat method as implemented in GSL
 *     George Marsaglia, Wai Wan Tsang
 *     The Ziggurat Method for Generating Random Variables
 *     Journal of Statistical Software, vol. 5 (2000), no. 8
 *     http://www.jstatsoft.org/v05/i08/
 */

extern const double ytab[256];
extern const double ktab[256];
extern const double wtab[256];

inline double rand_normal(xorshift64 &rng, double mu, double sigma) {
	const double R = 3.6554204190269415;
	double x, y;

	for(;;) {
		uint64_t u = rng.get_uint64();
		// use the top 8 high bits for b
		uint64_t b = u >> 56;
		// use the rest for a
		int64_t a = static_cast<int64_t>(u << 8);
		double aa = static_cast<double>(a);

		x = aa * wtab[b];

		if (std::abs(aa) < ktab[b])
			break;

		if(b < 255) {
			y = ytab[b + 1] + (ytab[b] - ytab[b + 1]) * rng.get_double52();
		} else {
			x = R + rand_exp(rng,R);
			y = exp(-R * (x - 0.5 * R)) * rng.get_double52();
			// Use aa to pick the sign of x
			x = copysign(x,aa);
		}


			break;
    }
	return (mu + sigma * x);
}

inline double rand_abs_normal(xorshift64 &rng, double mu, double sigma) {
	const double R = 3.6554204190269415;
	double x, y;

	for(;;) {
		uint64_t u = rng.get_uint64();
		// use the top 8 high bits for b
		uint64_t b = u >> 56;
		// use the rest for a
		int64_t a = static_cast<int64_t>((u << 8) & UINT64_C(0x7FFFFFFFFFFFFFFF));
		double aa = static_cast<double>(a);

		x = aa * wtab[b]; //abs of aa ??

		if (aa < ktab[b])
			break;

		if(b < 255) {
			y = ytab[b + 1] + (ytab[b] - ytab[b + 1]) * rng.get_double52();
		} else {
			x = R + rand_exp(rng,R);
			y = exp(-R * (x - 0.5 * R)) * rng.get_double52();
		}

		if (y < exp(-0.5 * x * x))
			break;
    }
	return (mu + sigma * x);
}

#endif
