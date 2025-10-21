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
#include <vector>
#include <limits>
#include <algorithm>

namespace mfem_mgis {

  /*!
   * \brief Forward declaration of the Material structure from mfem_mgis.
   */

  struct Material;

}  // end of namespace mfem_mgis

namespace opera_hpc {

  /*!
   * \struct FirstPrincipalStressValueAndLocation
   * \brief A simple structure to store a stress value and its
   * 3D spatial location.
   */ 
  struct FirstPrincipalStressValueAndLocation {
    mfem_mgis::real value;
    std::array<mfem_mgis::real, 3u> location;
  };

  /*!
   * \struct StressValueAndLocation
   * \brief Alias for FirstPrincipalStressValueAndLocation.
   *
   * This provides a more generic name for storing a stress value
   * (which might not be the first principal stress) and its location,
   * while reusing the same data layout.
   */
  struct StressValueAndLocation : FirstPrincipalStressValueAndLocation{};

  /*!
   * \brief Finds the maximum first principal stress value and its location
   * within a material.
   * \param[in] m The material object to search.
   * \return A structure containing the maximum first principal stress
   * value and the coordinates where it occurs.
   */

  FirstPrincipalStressValueAndLocation findFirstPrincipalStressValueAndLocation(
      const mfem_mgis::Material &);
  
  /*!
  * \brief Retrieves the locations of all points within the material
  * where the first principal stress exceeds a given threshold.
  * \param[in] m The material object to query.
  * \param[in] threshold The stress threshold value.
  * \return A vector of 3D coordinates [x, y, z] for all points
  * exceeding the threshold.
  */
  std::vector<std::array<mfem_mgis::real, 3u>> getPointsAboveStressThreshold(
      const mfem_mgis::Material &, const mfem_mgis::real);
/*!
   * \brief Retrieves both the stress value and location for all points
   * within the material where the first principal stress exceeds
   * a given threshold.
   * \param[in] m The material object to query.
   * \param[in] v The stress threshold value.
   * \return A vector of StressValueAndLocation structures, each
   * containing the stress value and coordinates for a point
   * exceeding the threshold.
   */
  std::vector<StressValueAndLocation>
  getPointsandStressAboveStressThreshold(const mfem_mgis::Material &m,
                                         const mfem_mgis::real v);

}  // end of namespace opera_hpc

#endif /* LIB_OPERA_HPC_UTILITIES_HXX */
