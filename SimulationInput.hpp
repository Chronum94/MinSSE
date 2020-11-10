#pragma once
#include <fstream>
#include <iostream>
#include <string>

// TODO: Don't be bad.
template <class IntegerType, class FloatType> struct SimulationInput {
  SimulationInput() {
    nx = 4;
    ny = 4;
    nz = 4;
    beta = 1;
    dilution = 0.0;

    nbins = 10;
    msteps = 10;
    isteps = 10;
    max_sites = nx * ny * nz;
  }

  SimulationInput(uint16_t _nx, uint16_t _ny, uint16_t _nz, FloatType _beta,
                  FloatType _dilution, IntegerType _nbins = 10,
                  IntegerType _msteps = 1000, IntegerType _isteps = 1000) {
    nx = _nx;
    ny = _ny;
    nz = _nz;
    beta = _beta;
    dilution = _dilution;

    nbins = _nbins;
    msteps = _msteps;
    isteps = _isteps;
    max_sites = nx * ny * nz;
    n_active_sites = static_cast<IntegerType>((1.0 - dilution) * max_sites);

    seed = 0;

    if (1.0 - dilution < 2.0 / max_sites) {
      std::cout << "You have a VERY diluted lattice!\nThere maybe no "
                   "bonds/sites.\nThis not NOT meaningful!";
    }
  }

  uint16_t nx;
  uint16_t ny;
  uint16_t nz;
  FloatType beta;
  FloatType dilution;

  IntegerType nbins;
  IntegerType msteps;
  IntegerType isteps;
  long unsigned seed;
  IntegerType max_sites;
  IntegerType n_active_sites;

  std::string resultsfile{"minsse_results.dat"};
  // friend std::ostream& operator<<(std::ostream& out, const
  // SimulationInput<IntegerType, FloatType>& in); template <class IntegerType,
  // class FloatType>
  friend std::ostream &
  operator<<(std::ostream &out,
             const SimulationInput<IntegerType, FloatType> &in) {
    out << "========================================\n";
    out << "Lattice: " << in.nx << " x " << in.ny << " x " << in.nz << " = "
        << in.max_sites << " max sites\n";
    out << "Beta: " << in.beta << "\n";
    out << "Dilution: " << in.dilution << ", exactly " << in.n_active_sites
        << " active sites.\n";
    out << "Bin count: " << in.nbins << "\n";
    out << "MC sweeps in each bin: " << in.msteps << "\n";
    out << "Equilibriation sweeps: " << in.isteps << "\n";
    out << "seed: " << in.seed << "\n";
    out << "========================================\n";
    out << "TotSites: " << in.max_sites << "\n";
    out << "ActSites: " << in.n_active_sites << "\n";
    return out;
  }
};
