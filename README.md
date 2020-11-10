# MinSSE: A minimal 2/3D Heisenberg lattice simulator using the stochastic series expansion.

Minimal compilation:

`clang++ sse_driver.cpp --std=c++17 -Ofast -march=native -mtune=native -g -fno-unroll-loops -malign-double -DNDEBUG DRANDOM_GENERATOR` where
`RANDOM_GENERATOR` can be `XORSHIFT128`, `XORSHIFT32` (don't really recommend this...), `SQUARES`, `AES4`, or `LCG`.

Minimal run with:

`./a.out 4 4 4 2.0 0.0 100 1000 1000 dummy.dat`

The above line runs a 4x4x4 lattice at inverse temperature (beta) of 2, 0.0 dilution fraction, gathering 100 bins of data, each of 1000 Monte Carlo sweeps. 1000 equilibriation sweeps are carried out prior to taking mesurements.

The run signature is:

`./a.out nx ny nz beta dilution nbins msteps isteps outfilename`
# Dependencies:

1. This project uses [`SmallPRNG`](https://github.com/DKenefake/SmallPRNG) as fast and reliable PRNG source.
