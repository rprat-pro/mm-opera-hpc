/*!
 * \file   UniaxialMacroscopicStressPeriodicSimulation.hxx
 * \brief
 * \author Thomas Helfer
 * \date   16/09/2025
 */

#ifndef LIB_MM_OPERA_HPC_UNIAXIALMACROSCOPICSTRESSPERIODICSIMULATION_HXX
#define LIB_MM_OPERA_HPC_UNIAXIALMACROSCOPICSTRESSPERIODICSIMULATION_HXX

#include <array>
#include <vector>
#include <optional>
#include <functional>
#include <string_view>
#include "MM_OPERA_HPC/Config.hxx"
#include "MM_OPERA_HPC/MacroscropicElasticMaterialProperties.hxx"

namespace mfem_mgis {

// forward declaration
struct PeriodicNonLinearEvolutionProblem;

} // end of namespace mfem_mgis

namespace mm_opera_hpc {

MM_OPERA_HPC_EXPORT void print_memory_footprint(std::string_view) noexcept;

struct MM_OPERA_HPC_EXPORT  MacroscopicVariables {
  // diagonal part of the deformation gradient
  std::array<mfem_mgis::real, 3u> F = std::array<mfem_mgis::real, 3u>{1, 1, 1};
  // increment of the diagonal part of the deformation gradient during the
  // previous time step
  std::array<mfem_mgis::real, 3u> dF = std::array<mfem_mgis::real, 3u>{};
  // diagonal components of the macroscropic Cauchy stress
  std::array<mfem_mgis::real, 3u> S = std::array<mfem_mgis::real, 3u>{};
  // previous succesful time step
  std::optional<mfem_mgis::real> previous_time_increment;
};

/*!
 * \brief the role of this class is to perform a periodic simulation where:
 *
 * - the macroscopic stress state is uniaxial along the zz axis (up to a given
 *   tolerance)
 * - the axial componenent Fzz of the macroscopic deformation gradient is
 *   prescribed
 * - the off-diagonal components of the macroscopic deformation gradient are
 *   set to zero
 *
 * The loading is discretized by the user into temporal sequences. Each
 * temporal sequence can be automatically divided into smaller time steps in
 * case of divergence of the resolution.
 *
 * Post-processings are performed at the end of the temporal sequences
 */
struct MM_OPERA_HPC_EXPORT UniaxialMacroscopicStressPeriodicSimulation {
  /*!
   * \brief constructor
   * \param[in] p: problem to be solved
   */
  UniaxialMacroscopicStressPeriodicSimulation(
      mfem_mgis::PeriodicNonLinearEvolutionProblem &,
      const std::function<mfem_mgis::real(mfem_mgis::real)>,
      const MacroscropicElasticMaterialProperties &, const bool) noexcept;
  /*!
   * \brief run the simulation over the given temporal sequences
   * \param[in] temporal_sequences: list of times defining the temporal
   * sequences.
   * \pre the times must be sorted in ascending order
   * \return true on success, false otherwise
   */
  [[nodiscard]] bool run(const std::vector<mfem_mgis::real> &) noexcept;

private:
  //! \brief underlying problem
  mfem_mgis::PeriodicNonLinearEvolutionProblem &problem;
  //
  std::function<mfem_mgis::real(mfem_mgis::real)>
      imposed_axial_deformation_gradient_value;
  //
  MacroscopicVariables macroscopic_variables;
  //! \brief macroscopic material parameters
  const MacroscropicElasticMaterialProperties macroscopic_elastic_material_properties;
  //
  const bool post_processing;
};

} // namespace mm_opera_hpc

#endif /* LIB_MM_OPERA_HPC_UNIAXIALMACROSCOPICSTRESSPERIODICSIMULATION_HXX */
