/*!
 * \file   UniaxialMacroscopicStressPeriodicSimulation.cxx
 * \brief
 * \author Thomas Helfer
 * \date   16/09/2025
 */

#include <cmath>
#include <sys/resource.h>
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MM_OPERA_HPC/UniaxialMacroscopicStressPeriodicSimulation.hxx"

namespace mm_opera_hpc {

// display information
static long get_memory_checkpoint() {
  rusage obj;
  int who = 0;
  [[maybe_unused]] auto test = getrusage(who, &obj);
  assert((test = -1) && "error: getrusage has failed");
  long res;
  MPI_Reduce(&(obj.ru_maxrss), &(res), 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  return res;
};

void print_memory_footprint(std::string_view msg) noexcept {
  long mem = get_memory_checkpoint();
  double m = double(mem) * 1e-6; // conversion kb to Gb
  mfem_mgis::Profiler::Utils::Message(msg, " memory footprint: ", m, " GB");
}

[[nodiscard]] static const std::array<mfem_mgis::real, 3u>
computeMacroscopicCauchyStress(
    mfem_mgis::PeriodicNonLinearEvolutionProblem &problem,
    const std::array<mfem_mgis::real, 3u> &F) {
  // integrals of the diagonal components of the First Piola-Kirchhoff stress
  auto pk1_integral = std::array<mfem_mgis::real, 3u>{};
  auto volume = mfem_mgis::real{};
  // get pk1 integral and volume on all materials //
  auto [pk1_integrals, volumes] =
      mfem_mgis::computeMeanThermodynamicForcesValues<true>(
          problem.template getImplementation<true>());
  // sum volume of all materials
  // sum pk1 integrals
  for (const auto &m : problem.getAssignedMaterialsIdentifiers()) {
    volume += volumes[m];
    pk1_integral[0] += pk1_integrals[m][0];
    pk1_integral[1] += pk1_integrals[m][1];
    pk1_integral[2] += pk1_integrals[m][2];
  }
  assert(volume > 0.0);
  // -- sum volumes of all processess
  MPI_Allreduce(MPI_IN_PLACE, &volume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  // -- sum integrals from all processes
  MPI_Allreduce(MPI_IN_PLACE, pk1_integral.data(), 3, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  // macroscopic values of the First Piola-Kirchhoff stress
  const auto pk1 = std::array<mfem_mgis::real, 3u>{pk1_integral[0] / volume,
                                                   pk1_integral[1] / volume,
                                                   pk1_integral[2] / volume};
  // macroscopic change of volume
  const auto det = F[0] * F[1] * F[2];
  // computation of the diagonal components of the Cauchy stress
  return {pk1[0] * F[0] / det, pk1[1] * F[1] / det, pk1[2] * F[2] / det};
} // end of computeMacroscopicCauchyStress

[[nodiscard]] static bool simulateOverATimeStep(
    mfem_mgis::PeriodicNonLinearEvolutionProblem &problem,
    MacroscopicVariables &macroscopic_variables, //
    const std::function<mfem_mgis::real(mfem_mgis::real)>& imposed_axial_deformation_gradient_value,
    const MacroscropicElasticMaterialProperties &mp, const bool post_processing, 
    const mfem_mgis::real bts, //
    const mfem_mgis::real ets) {
  using namespace mfem_mgis::Profiler::Utils;
  Message("Solving time step from ", bts, " to ", ets);
  // --- fixed-point param
  const double tolFP = 1.e4;
  const int maxitFP = 50;
  // Traction
  const auto dt = ets - bts;

  auto &F = macroscopic_variables.F;
  auto &S = macroscopic_variables.S;
  // deformation gradient at the beginning of the time step
  auto S0 = S;
  auto F0 = F;

  if (macroscopic_variables.previous_time_increment.has_value()) {
    auto &dF = macroscopic_variables.dF;
    const auto pdt = *(macroscopic_variables.previous_time_increment);
    // -- update F from input data, here def = 5e-4. -- //
    F[0] += dF[0] * (dt / pdt);
    F[1] += dF[1] * (dt / pdt);
  }
  F[2] = imposed_axial_deformation_gradient_value(ets);
  // --- fixed point iteration number
  int itFP = 0;
  // --- fixed point loop --- //
  bool converged = false;
  while ((!converged) && (itFP <= maxitFP)) {
    // -- compute correction -- //
    // -- note that S[0] and S[1] are updated at the end of this loop section,
    // default is 0 -- //
    F[0] += ((mp.nueff * mp.nueff - 1) * S[0]) / mp.Eeff +
            (mp.nueff * (mp.nueff + 1) * S[1]) / mp.Eeff;
    F[1] += ((mp.nueff * mp.nueff - 1) * S[1]) / mp.Eeff +
            (mp.nueff * (mp.nueff + 1) * S[0]) / mp.Eeff;
    // -- note that MacroscopicGradientsEvolution taktes into account of new
    // F[0] and F[1] -- //
    auto statistics = problem.solve(bts, dt);
    if (!statistics.status) {
      Message("error: the resolution failed.");
      S = S0;
      F = F0;
      problem.revert();
      return false;
    }

    S = computeMacroscopicCauchyStress(problem, F);
    const auto r = std::sqrt(S[0] * S[0] + S[1] * S[1]);
    //
    Message("Fixed Point iteration", itFP, ": |res| =", r);
    converged = std::abs(r) < tolFP;
    itFP++;
  }
  //
  Message("Solution at time", ets, ":", F[2], S[0], S[1], S[2]);
  //
  if (itFP >= maxitFP) {
    Message("error: maximum number of iterations for the "
            "fixed-point algorithm reached");
    S = S0;
    F = F0;
    problem.revert();
    return false;
  }
  //
  if (post_processing) {
    // add postprocessing
    CatchTimeSection("common::post_processing_step");
    problem.executePostProcessings(bts, dt);
  }
  // update state variable for the next time step
  problem.update();
  // update information for time extrapolation
  auto &dF = macroscopic_variables.dF;
  macroscopic_variables.previous_time_increment = dt;
  for (std::size_t i = 0; i != 3; ++i) {
    dF[i] = F[i] - F0[i];
  }
  return true;
}

[[nodiscard]] static bool simulateOverATemporalSequence(
    mfem_mgis::PeriodicNonLinearEvolutionProblem &problem,
    MacroscopicVariables &macroscopic_variables, //
    const std::function<mfem_mgis::real(mfem_mgis::real)> &imposed_axial_deformation_gradient_value,
    const bool post_processing, const MacroscropicElasticMaterialProperties &mp,
    const mfem_mgis::real bts, //
    const mfem_mgis::real ets) {
  using namespace mfem_mgis::Profiler::Utils;
  constexpr auto maximum_of_substeps = std::size_t{10};
  auto nsubsteps = std::size_t{};
  auto nstep = std::size_t{1};
  auto t = bts;
  auto dt = ets - bts;
  while (nstep != 0) {
    const auto do_post_processing = post_processing && (nstep == 1);
    const auto success =
        simulateOverATimeStep(problem, macroscopic_variables,
                              imposed_axial_deformation_gradient_value, mp,
                              do_post_processing, t, t + dt);
    if (success) {
      t += dt;
      --nstep;
    } else {
      ++nsubsteps;
      if (nsubsteps == maximum_of_substeps) {
        Message("Maximum number of substeps reached for temporal sequence from",
                bts, "to", ets);
        return false;
      } else {
        Message("dividing time step by 2");
      }
      dt *= mfem_mgis::real{0.5};
      nstep *= 2;
    }
  }
  //
  print_memory_footprint(
      "[At end of temporal sequence: " + std::to_string(ets) + "]");
  return true;
} // end of simulateOverATemporalSequence

UniaxialMacroscopicStressPeriodicSimulation::
    UniaxialMacroscopicStressPeriodicSimulation(
        mfem_mgis::PeriodicNonLinearEvolutionProblem &p,
        const std::function<mfem_mgis::real(mfem_mgis::real)> Fzz,
        const MacroscropicElasticMaterialProperties &mp, const bool b) noexcept
    : problem(p), imposed_axial_deformation_gradient_value(Fzz),
      macroscopic_elastic_material_properties(mp), post_processing(b) {
  // warning Fxx, F[1], and Fzz are defined during the following loop
  problem.setMacroscopicGradientsEvolution([this](const double) {
    const auto &F = this->macroscopic_variables.F;
    auto ret = std::vector<mfem_mgis::real>(9, mfem_mgis::real{});
    std::copy(F.begin(), F.end(), ret.begin());
    return ret;
  });
} // end of UniaxialMacroscopicStressPeriodicSimulation

bool UniaxialMacroscopicStressPeriodicSimulation::run(
    const std::vector<mfem_mgis::real> &temporal_sequences) noexcept {
  for (std::size_t i = 0; i + 1 != temporal_sequences.size(); i++) {
    if (!simulateOverATemporalSequence(
            this->problem, this->macroscopic_variables,
            this->imposed_axial_deformation_gradient_value,
            this->post_processing,
            this->macroscopic_elastic_material_properties,
            temporal_sequences[i], temporal_sequences[i + 1])) {
      return false;
    }
  }
  return true;
}

} // end of namespace mm_opera_hpc