#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>

template <class StringType> Simulation readinput(StringType filename) {
  std::ifstream inputfile(filename);
  std::string line;
  unsigned int num_lines = 0;
  while (std::getline(inputfile, line)) {
    if (line.length() == 0) {
      continue;
    }
    std::cout << line << '\n';
    auto linestream = std::stringstream{line};
    std::string current_token;
    switch
      num_lines

          ++ num_lines;
  }
  std::cout << num_lines;
  if (num_lines != 7) {
    std::cout << "Invalid input file! Input files need 7 valid lines.\n";

    std::cout << "nx ny nz # 3 integers; lattice size along each dimension.\n";
    std::cout << "beta # Float/double; dimensionless inverse temperature in "
                 "units of J.\n";
    std::cout << "nbins # Positive integer for the binning when calculating "
                 "expectation values.\n";
    std::cout << "isteps # Positive integer for the ???.\n";
    std::cout << "seed # Positive integer, seed for the pseudorandom number "
                 "generator.\n";
    throw;
  }

  return Simulation();
}