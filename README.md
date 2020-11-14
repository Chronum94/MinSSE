# MinSSE: A minimal 2/3D Heisenberg lattice simulator using the stochastic series expansion.

Minimal compilation:

`clang++ sse_driver.cpp --std=c++17 -Ofast -DNDEBUG -DRANDOM_GENERATOR` where
`RANDOM_GENERATOR` can be `XORSHIFT128`, `XORSHIFT32` (don't really recommend this...), `SQUARES`, `AES4`, or `LCG`. You can also replace `clang++` with `g++` and it should compile.

Not-so-minimal compilation:

`clang++ sse_driver.cpp --std=c++17 -Ofast -march=native -mtune=native -g -fno-unroll-loops -malign-double -DNDEBUG -DRANDOM_GENERATOR`

The general run signature is:

`./a.out nx ny nz beta dilution nbins msteps isteps int32_seed`

`nx`, `ny`, `nz`: Lattice sizes along each direction. Periodic boundaries are implemented in all directions. If `nz=1`, you have a 2D lattice.\
`beta`: Inverse temperature.\
`dilution`: Percentage of vacancies in lattice.\
`nbins`: Number of bins in which to bin the run data.\
`msteps`: Number of Monte-Carlo sweeps per bin. A measurement is done after each sweep.\
`isteps`: Number of equilibriation sweeps. This is to bring the number of operators to an equilibriated number.\
`int32_seed`: An int32 number, which SmallPRNG reads in, mixes a bit, and use _that_ as a seed. It's still all deterministic.

Minimal run with (assuming you're currently in the directory with the `a.out` executable):

```
mkdir TestDir
cd TestDir
../a.out 4 4 1 2.0 0.0 100 1000 1000 234545234
```

The above line runs a 4x4x1 lattice at inverse temperature (beta) of 2, 0.0 dilution fraction, gathering 100 bins of data, each of 1000 Monte Carlo sweeps. 1000 equilibriation sweeps are carried out prior to taking measurements.

The above run will output a few files:

`spinmask.out`: A 'mask' of the spins. 1 corresponds to occupied, 0 to a vacancy.\
`spinbonds.out`: The number of bonds on a spin.\
`bonds.out`: Bonds and the index of the spins it connects.\
`out.log`: Some basic calculation details.\
`locsusc.dat`: The local susceptibility for each spin in the lattice. Each row is a bin.\
`results.dat`: The uniform susceptibility and total energy of the lattice for each bin.

# Dependencies:

1. This project uses [`SmallPRNG`](https://github.com/DKenefake/SmallPRNG) as fast and reliable PRNG source.

# References:

1. [Sandvik, Anders W. "Computational studies of quantum spin systems." AIP Conference Proceedings. Vol. 1297. No. 1. American Institute of Physics, 2010.](https://aip.scitation.org/doi/abs/10.1063/1.3518900?casa_token=eM1OTCjSqtYAAAAA:8cc3nsljVlsP3A3Fus4rd4VmW6zXLqtJ_G-CKE-9Gp7XC2oECrho5_Z2bpRkGu8KGqFHFmdsgdBG) Section 5.
