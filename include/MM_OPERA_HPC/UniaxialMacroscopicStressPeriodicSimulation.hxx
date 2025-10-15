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

}  // end of namespace mfem_mgis

namespace mm_opera_hpc {

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
     * \brief numerical parameters used to find both the microscopic
     * micromorphic satisfying the mechanical equilibrium and the macroscopic
     * deformation gradient to find a macroscopic uniaxial tensile stress state
     */
    struct MM_OPERA_HPC_EXPORT NumericalParameters {
      const mfem_mgis::size_type microscopic_equilibrium_verbosity_level =
          mfem_mgis::size_type{1};
      const mfem_mgis::size_type
          microscopic_equilibrium_maximum_number_of_iterations{20};
      const mfem_mgis::real microscopic_equilibrium_relative_tolerance =
          mfem_mgis::real{1e-4};
      const mfem_mgis::real microscopic_equilibrium_absolute_tolerance =
          mfem_mgis::real{0};
      const bool switch_microscopic_equilibrium_convergence_criterion = true;
      //
      //! \brief macroscopic material parameters
      const MacroscropicElasticMaterialProperties
          macroscopic_elastic_material_properties;
      const mfem_mgis::real macroscopic_stress_absolute_tolerance =
          mfem_mgis::real{1e4};
      const size_t macroscopic_stress_maximum_number_of_iterations = 20;
    };
    //! \brief an internal auxiliary data structure
    struct MM_OPERA_HPC_EXPORT MacroscopicUnknowns {
      // diagonal part of the deformation gradient
      std::array<mfem_mgis::real, 3u> F =
          std::array<mfem_mgis::real, 3u>{1, 1, 1};
      // increment of the diagonal part of the deformation gradient during the
      // previous time step
      std::array<mfem_mgis::real, 3u> dF = std::array<mfem_mgis::real, 3u>{};
      // diagonal components of the macroscropic Cauchy stress
      std::array<mfem_mgis::real, 3u> S = std::array<mfem_mgis::real, 3u>{};
      // previous succesful time step
      std::optional<mfem_mgis::real> previous_time_increment;
    };
    /*!
     * \brief constructor
     * \param[in] p: problem to be solved
     */
    UniaxialMacroscopicStressPeriodicSimulation(
        mfem_mgis::PeriodicNonLinearEvolutionProblem &,
        const std::function<mfem_mgis::real(mfem_mgis::real)>,
        const NumericalParameters &,
        const bool) noexcept;
    /*!
     * \brief run the simulation over the given temporal sequences
     * \param[in] out: output files for the macroscopic results (diagonal
     * components of the deformation gradient and diagonal components of the
     * Cauchy stress).
     * \param[in] temporal_sequences: list of times defining the
     * temporal sequences.
     * \pre the times must be sorted in ascending order
     * \return true on success, false otherwise
     */
    [[nodiscard]] bool run(std::ostream &,
                           const std::vector<mfem_mgis::real> &) noexcept;

   private:
    //! \brief underlying problem
    mfem_mgis::PeriodicNonLinearEvolutionProblem &problem;
    //
    std::function<mfem_mgis::real(mfem_mgis::real)>
        imposed_axial_deformation_gradient_value;
    //
    MacroscopicUnknowns macroscopic_unknowns;
    //
    const NumericalParameters numerical_parameters;
    //
    const bool post_processing;
  };

}  // namespace mm_opera_hpc

#endif /* LIB_MM_OPERA_HPC_UNIAXIALMACROSCOPICSTRESSPERIODICSIMULATION_HXX */
