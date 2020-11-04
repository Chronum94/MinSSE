#include <numeric>
#include <vector>

template <class IType, class FType, class PRNGState>
std::vector<IType> make_lattice_points(IType nx, IType ny, IType nz,
                                       FType dilution, PRNGState prng) {
  auto total_spins = nx * ny * nz;

  std::vector<IType> dilution_mask(total_spins, 1);
  std::vector<IType> site_indices(total_spins);

  std::iota(site_indices.begin(), site_indices.end(), 0);
}

std::vector<IType> dilute_lattice(std::vector<IType> spin_indices, std_vector)