/*!
 * \file   main.cxx
 * \brief    
 * \author Thomas Helfer
 * \date   29/09/2024
 */

#include <cstdlib>
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/UniformImposedPressureBoundaryCondition.hxx"
#include "OperaHPC/BubbleDescription.hxx"
#include "OperaHPC/Utilities.hxx"

namespace opera_hpc {

  struct Bubble : BubbleDescription {
    bool broken = false;
  };

  bool areAllBroken(const std::vector<Bubble>& bubbles) {
    for (const auto& b : bubbles) {
      if (!b.broken) {
        return false;
      }
    }
    return true;
  }
}  // end of namespace opera_hpc

// add postprocessing
template <typename Problem>
void post_process(Problem& p, double start, double end) {
  CatchTimeSection("common::post_processing_step");
  p.executePostProcessings(start, end);
}

int main(int argc, char** argv) {
  using namespace mfem_mgis::Profiler::Utils; // Use Message
  // options treatment
  mfem_mgis::initialize(argc, argv);
  mfem_mgis::Profiler::timers::init_timers();
  // reference pressure
  constexpr auto pref = mfem_mgis::real{1e6};
  constexpr auto maximumNumberSteps = mfem_mgis::size_type{100};
  // distance pour déterminer si une bulle casse
  constexpr auto dmin = mfem_mgis::real{1};
  //
  auto bubbles = [] {
    auto r = std::vector<opera_hpc::Bubble>{};
    for (const auto& d : opera_hpc::BubbleDescription::read("bubbles.txt")) {
      auto b = opera_hpc::Bubble{};
      static_cast<opera_hpc::BubbleDescription&>(b) = d;
      r.push_back(b);
    }
    return r;
  }();
  // finite element space
  auto fed = std::make_shared<mfem_mgis::FiniteElementDiscretization>(
      mfem_mgis::Parameters{
          {"MeshFileName", "mesh/mesh_sphere.msh"},
          {"FiniteElementFamily", "H1"},
          {"FiniteElementOrder", 1},
          {"UnknownsSize", 3},
          {"NumberOfUniformRefinements", 0}, // faster for testing
          //{"NumberOfUniformRefinements", parameters.parallel ? 1 : 0},
          //          {"Hypothesis", "Tridimensional"},
          {"Parallel", true}});
  // definition of the nonlinear problem
  auto problem = mfem_mgis::PeriodicNonLinearEvolutionProblem{fed};
  // macroscopic strain
  std::vector<mfem_mgis::real> e(6, mfem_mgis::real{0});
  problem.setMacroscopicGradientsEvolution(
      [e](const mfem_mgis::real) { return e; });
  // imposing pressure on the bubble boundaries
  for (const auto &b : bubbles) {
    problem.addBoundaryCondition(
        std::make_unique<mfem_mgis::UniformImposedPressureBoundaryCondition>(
            problem.getFiniteElementDiscretizationPointer(),
            b.boundary_identifier, [&b](const mfem_mgis::real) {
              // pressure evolution in the bubble
              // if sound, return the reference pressure
              // if broken, return 0
              return b.broken ? 0 : pref;
            }));
    }
    // choix du solver linéaire +
    const int verbosity = 0;
    const int post_processing = 0;
    auto solverParameters = mfem_mgis::Parameters{};
    solverParameters.insert(
        mfem_mgis::Parameters{{"VerbosityLevel", verbosity}});
    //  solverParameters.insert(mfem_mgis::Parameters{{"MaximumNumberOfIterations",
    //  defaultMaxNumOfIt}});
    //  solverParameters.insert(mfem_mgis::Parameters{{"Tolerance", Tol}});
    //
    auto options = mfem_mgis::Parameters{{"VerbosityLevel", verbosity},
                                         {"Strategy", "Elasticity"}};
    auto preconditionner =
        mfem_mgis::Parameters{{"Name", "HypreBoomerAMG"}, {"Options", options}};
    solverParameters.insert(
        mfem_mgis::Parameters{{"Preconditioner", preconditionner}});
    problem.setLinearSolver("HyprePCG", solverParameters);
    // matrice élastique
    problem.addBehaviourIntegrator("Mechanics", 1, "src/libBehaviour.so",
                                   "Elasticity");
    auto &m = problem.getMaterial(1);
    for (auto &s : {&m.s0, &m.s1}) {
      mgis::behaviour::setMaterialProperty(*s, "YoungModulus", 150e9);
      mgis::behaviour::setMaterialProperty(*s, "PoissonRatio", 0.3);
      mgis::behaviour::setExternalStateVariable(*s, "Temperature", 293.15);
      mgis::behaviour::setExternalStateVariable(*s, "Temperature", 293.15);
    }
    if (post_processing){
    problem.addPostProcessing("ParaviewExportResults",
                              {{"OutputFileName", "TestCaseOneBubble"}});
    }
    //
    auto nstep = mfem_mgis::size_type{};
    while ((nstep != maximumNumberSteps) &&
           (!opera_hpc::areAllBroken(bubbles))) {
      problem.solve(0, 1);
      const auto r = opera_hpc::findFirstPrincipalStressValueAndLocation(
          problem.getMaterial(1));
      auto nbroken = mfem_mgis::size_type{};
      auto all_broken_bubbles_identifiers = std::vector<mfem_mgis::size_type>{};
      for (auto &b : bubbles) {
        const auto d = opera_hpc::distance(b, r.location);
        if (d < dmin) {
          b.broken = true;
          ++nbroken;
        }
        if (b.broken) {
          all_broken_bubbles_identifiers.push_back(b.boundary_identifier);
        }
      }
      if (nbroken == 0) {
        Message("no bubble broke at step ", nstep,"\n");
        break;
      }


      Message("step:",nstep);
      Message("-", nbroken, "bubbles broke at this step.");
      Message("-", all_broken_bubbles_identifiers.size(), "bubbles are broken");
      Message("- value of the first principal stress:", r.value, "at coordinate (", r.location[0], ",", r.location[1], "," , r.location[2], ")");
      if (post_processing)
        post_process(problem, 0+double(nstep), 1+double(nstep));
      
      ++nstep;
    }
    mfem_mgis::Profiler::timers::print_and_write_timers();
    return EXIT_SUCCESS;
}
