/*!
 * \file   OperaHPC/Utilities.hxx
 * \brief
 * \author Thomas Helfer
 * \date   29/09/2024
 */

#ifndef LIB_OPERA_HPC_UTILITIES_HXX
#define LIB_OPERA_HPC_UTILITIES_HXX

#include <array>
#include "MFEMMGIS/Config.hxx"

namespace mfem_mgis {

  struct Material;

}  // end of namespace mfem_mgis

namespace opera_hpc {

  struct FirstPrincipalStressValueAndLocation {
    mfem_mgis::real value;
    std::array<mfem_mgis::real, 3u> location;
  };

  FirstPrincipalStressValueAndLocation findFirstPrincipalStressValueAndLocation(
      const mfem_mgis::Material &);

  std::vector<std::array<mfem_mgis::real, 3u>>
  getIntegrationPointLocationsWithFirstPrincipalStressGreaterThanThresold(
      mfem_mgis::Material &, const mfem_mgis::real);

}  // end of namespace opera_hpc

#endif /* LIB_OPERA_HPC_UTILITIES_HXX */
