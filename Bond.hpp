#pragma once

template <class IntegerType> struct Bond {
  IntegerType s1;
  IntegerType s2;
  void invert_spins() {
    s1 *= -1;
    s2 *= -1;
  }

  inline bool are_spins_parallel() { return s1 & s2; }
};