/*!
 * \file   OperaHPC/Utilities.hxx
 * \brief
 * \author Thomas Helfer
 * \date   29/09/2024
 */

#ifndef LIB_OPERA_HPC_UTILITIES_HXX
#define LIB_OPERA_HPC_UTILITIES_HXX

#include "MFEMMGIS/Config.hxx"
#include <array>

namespace mfem_mgis {

struct Material;

} // end of namespace mfem_mgis

namespace opera_hpc {

struct FirstPrincipalStressValueAndLocation {
  mfem_mgis::real value;
  std::array<mfem_mgis::real, 3u> location;
};

FirstPrincipalStressValueAndLocation
findFirstPrincipalStressValueAndLocation(const mfem_mgis::Material &);

std::vector<std::array<mfem_mgis::real, 3u>>
getPointsAboveStressThreshold(mfem_mgis::Material &, const mfem_mgis::real);

std::vector<std::pair<mfem_mgis::real, std::array<mfem_mgis::real, 3u>>>
getPointsandStressAboveStressThreshold(mfem_mgis::Material &m,
                                       const mfem_mgis::real v);

} // end of namespace opera_hpc

#endif /* LIB_OPERA_HPC_UTILITIES_HXX */
