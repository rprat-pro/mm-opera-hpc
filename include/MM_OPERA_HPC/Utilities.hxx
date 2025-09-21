/*!
 * \file   Utilities.hxx
 * \brief
 * \author Thomas Helfer
 * \date   21/09/2025
 */

#ifndef LIB_MM_OPERA_HPC_UTILITIES_HXX
#define LIB_MM_OPERA_HPC_UTILITIES_HXX

#include <string_view>
#include "MM_OPERA_HPC/Config.hxx"

namespace mfem_mgis {

  struct PeriodicNonLinearEvolutionProblem;

}  // end of namespace mfem_mgis

namespace mm_opera_hpc {

  MM_OPERA_HPC_EXPORT void printMeshInformation(
      mfem_mgis::PeriodicNonLinearEvolutionProblem &);

  MM_OPERA_HPC_EXPORT void printMemoryFootprint(std::string_view) noexcept;

}  // namespace mm_opera_hpc

#endif /* LIB_MM_OPERA_HPC_UTILITIES_HXX */
