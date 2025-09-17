/*!
 * \file   GrainOrientations.hxx
 * \brief    
 * \author Thomas Helfer
 * \date   16/09/2025
 */

#ifndef LIB_MM_OPERA_HPC_GRAINORIENTATIONS_HXX
#define LIB_MM_OPERA_HPC_GRAINORIENTATIONS_HXX

#include <array>
#include <vector>
#include <string_view>
#include "MM_OPERA_HPC/Config.hxx"

namespace mm_opera_hpc {

[[nodiscard]] MM_OPERA_HPC_EXPORT std::vector<std::array<mfem_mgis::real, 3u>>
    readVectorsFromFile(std::string_view);

} // end of namespace mm_opera_hpc

#endif /* LIB_MM_OPERA_HPC_GRAINORIENTATIONS_HXX */
