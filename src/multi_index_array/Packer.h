
// Functions for packing floats and double into IEEE-754 format to try to keep portability across platforms when saving lattice and MIA data
// The following code was copied from Brian "Beej Jorgensen" Hall's website, which has been supplied to the public domain.
// http://beej.us/guide/bgnet/examples/ieee754.c
// Original copyright notice: http://beej.us/guide/bgnet/output/html/singlepage/bgnet.html#copyright


#ifndef PACKER_H_INCLUDED
#define PACKER_H_INCLUDED


#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>

#define pack754_32(f) (pack754_LibMIA((f), 32, 8))
#define pack754_64(f) (pack754_LibMIA((f), 64, 11))
#define unpack754_32(i) (unpack754_LibMIA((i), 32, 8))
#define unpack754_64(i) (unpack754_LibMIA((i), 64, 11))

typedef long double pack754_type;


inline uint64_t pack754_LibMIA(pack754_type f, unsigned bits, unsigned expbits)
{
    pack754_type fnorm;
    int shift;
    long long sign, exp, significand;
    unsigned significandbits = bits - expbits - 1; // -1 for sign bit

    if (f == 0.0) return 0; // get this special case out of the way

    // check sign and begin normalization
    if (f < 0) { sign = 1; fnorm = -f; }
    else { sign = 0; fnorm = f; }

    // get the normalized form of f and track the exponent
    shift = 0;
    while(fnorm >= 2.0) { fnorm /= 2.0; shift++; }
    while(fnorm < 1.0) { fnorm *= 2.0; shift--; }
    fnorm = fnorm - 1.0;

    // calculate the binary form (non-float) of the significand data
    significand = fnorm * ((1LL<<significandbits) + 0.5f);

    // get the biased exponent
    exp = shift + ((1<<(expbits-1)) - 1); // shift + bias

    // return the final answer
    return (sign<<(bits-1)) | (exp<<(bits-expbits-1)) | significand;
}

inline pack754_type unpack754_LibMIA(uint64_t i, unsigned bits, unsigned expbits)
{
    pack754_type result;
    long long shift;
    unsigned bias;
    unsigned significandbits = bits - expbits - 1; // -1 for sign bit

    if (i == 0) return 0.0;

    // pull the significand
    result = (i&((1LL<<significandbits)-1)); // mask
    result /= (1LL<<significandbits); // convert back to float
    result += 1.0f; // add the one back on

    // deal with the exponent
    bias = (1<<(expbits-1)) - 1;
    shift = ((i>>significandbits)&((1LL<<expbits)-1)) - bias;
    while(shift > 0) { result *= 2.0; shift--; }
    while(shift < 0) { result /= 2.0; shift++; }

    // sign it
    result *= (i>>(bits-1))&1? -1.0: 1.0;

    return result;
}

#endif
