#pragma once
#include <fstream>
#include <iostream>
#include <string>

// TODO: Don't be bad.
template <class IntegerType, class FloatType> struct SimulationInput {
  uint16_t nx;
  uint16_t ny;
  FloatType beta;
  FloatType dilution;

  IntegerType nbins;
  IntegerType msteps;
  IntegerType isteps;
  long unsigned seed;
  IntegerType max_sites = nx * ny;
  FloatType beta_n = max_sites * beta;
  IntegerType beta_n_integer = static_cast<IntegerType>(beta_n);
  // friend std::ostream& operator<<(std::ostream& out, const
  // SimulationInput<IntegerType, FloatType>& in); template <class IntegerType,
  // class FloatType>
  friend std::ostream &
  operator<<(std::ostream &out,
             const SimulationInput<IntegerType, FloatType> &in) {
    out << "========================================\n";
    out << " _ _          _  _  _\n";
    out << "| | | | |\\ | |_ |_ |_\n";
    out << "|   | | | \\|  _| _||_\n";
    out << "Lattice: " << in.nx << " x " << in.ny << " = " << in.max_sites
        << " sites\n";
    out << "Beta: " << in.beta << "\n";
    out << "Dilution: " << in.dilution << "\n";
    out << "Bin size: " << in.nbins << "\n";
    out << "MC sweeps in each bin: " << in.msteps << "\n";
    out << "Equilibriation sweeps: " << in.isteps << "\n";
    out << "seed: " << in.seed << "\n";
    out << "========================================\n";
    return out;
  }
};
