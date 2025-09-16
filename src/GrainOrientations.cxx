/*!
 * \file   GrainOrientations.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   16/09/2025
 */

#include <cmath>
#include <fstream>
#include "MM_OPERA_HPC/GrainOrientations.hxx"

namespace mm_opera_hpc{

[[nodiscard]] static mfem_mgis::real
norm(const std::array<mfem_mgis::real, 3u> &u) {
  mfem_mgis::real ret = 0.;
  for (auto val : u) {
    ret += val * val;
  }
  return (std::sqrt(ret));
}

std::vector<std::array<mfem_mgis::real, 3u>>
readVectorsFromFile(std::string_view filename) {
  // TODO: check if vectors size is equal to nMat
  std::vector<std::array<mfem_mgis::real, 3u>> vectors;
  std::ifstream inputFile(std::string{filename});

  auto checknorm = [](std::array<mfem_mgis::real, 3u> v) {
    mfem_mgis::real normv = norm(v);
    if (normv < 1.e-15) {
      throw std::invalid_argument("checknorm : vectors must not be null");
    }
    if (abs(normv - 1.) > 1.e-15) {
      std::cout << "checknorm : normalizing vector..." << std::endl;
      for (auto &comp : v) {
        comp *= 1. / normv;
      }
    }
    return v;
  };

  if (!inputFile.is_open()) {
    std::cout << "Erreur lors de l'ouverture du fichier " << filename
              << std::endl;
    return vectors;
  }

  std::string line;
  while (std::getline(inputFile, line)) {
    std::istringstream iss(line);
    std::array<mfem_mgis::real, 3u> vector;
    if (!(iss >> vector[0] >> vector[1] >> vector[2])) {
      std::cout << "Erreur de lecture du vecteur : " << line << std::endl;
      continue;
    }
    vectors.push_back(checknorm(vector));
  }

  inputFile.close();
  return vectors;
}

} // end of namespace mm_opera_hpc
