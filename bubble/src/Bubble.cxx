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
#include "Common/Print.hxx"
#include "Common/Memory.hxx"

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

struct TestParameters {
  const char* mesh_file = "mesh/single_sphere.msh";
  const char* behaviour = "Elasticity";
  const char* library = "src/libBehaviour.so";
  const char* bubble_file = "mesh/single_bubble.txt";
  const char* testcase_name = "TestCaseBubble";
  int parallel_mesh = 0;
  int order = 1;
  int refinement = 0;
  int post_processing = 1;  // default value : activated
  int verbosity_level = 0;  // default value : lower level
  mfem_mgis::real scale_factor_vp = 0.9;
};

void fill_parameters(mfem::OptionsParser& args, TestParameters& p) {
  args.AddOption(&p.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&p.parallel_mesh, "-pm", "--parallel-mesh",
                 "Parallel mesh format or not");
  args.AddOption(&p.library, "-l", "--library", "Material library.");
  args.AddOption(&p.bubble_file, "-f", "--bubble-file",
                 "File containing the bubbles.");
  args.AddOption(&p.order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&p.refinement, "-r", "--refinement",
                 "refinement level of the mesh, default = 0");
  args.AddOption(&p.post_processing, "-p", "--post-processing",
                 "run post processing step");
  args.AddOption(&p.verbosity_level, "-v", "--verbosity-level",
                 "choose the verbosity level");
  args.AddOption(&p.scale_factor_vp, "-sf", "--scale-factor-vp",
                 "Scaling factor for the principal stress");
  args.AddOption(&p.testcase_name, "-n", "--name-case", "Name of the testcase.");

  args.Parse();

  if (!args.Good()) {
    if (mfem_mgis::getMPIrank() == 0) args.PrintUsage(std::cout);
    mfem_mgis::finalize();
    exit(0);
  }
  if (p.mesh_file == nullptr) {
    if (mfem_mgis::getMPIrank() == 0)
      std::cout << "ERROR: Mesh file missing" << std::endl;
    args.PrintUsage(std::cout);
    mfem_mgis::abort(EXIT_FAILURE);
  }
  if (mfem_mgis::getMPIrank() == 0) args.PrintOptions(std::cout);
}

void write_message(std::ofstream& file_to_write, const auto&... args) {
  ((file_to_write << args << " "), ...) << std::endl;
}

int main(int argc, char** argv) {
  using namespace mfem_mgis::Profiler::Utils;  // Use Message
  // options treatment
  mfem_mgis::initialize(argc, argv);
  mfem_mgis::Profiler::timers::init_timers();

  // file collecting the output
  std::string file_name = "bubbles_and_stresses.txt";
  std::ofstream output_file(file_name);
  if (!output_file.is_open()) {
    std::cerr << "Failed to open the file: " << file_name << std::endl;
    return EXIT_FAILURE;  // Exit the program with an error code
  }

  // get parameters
  TestParameters p;
  mfem::OptionsParser args(argc, argv);
  fill_parameters(args, p);


  print_memory_footprint("[Start]");

  // reference pressure
  constexpr mfem_mgis::real pref = 1.0;
  constexpr mfem_mgis::size_type maximumNumberSteps =1000;
  // distance to determine if a bubble breaks
  constexpr mfem_mgis::real dmin = 0.600;
  //
  auto bubbles = [p] {
    auto r = std::vector<opera_hpc::Bubble>{};
    for (const auto& d : opera_hpc::BubbleDescription::read(p.bubble_file)) {
      auto b = opera_hpc::Bubble{};
      static_cast<opera_hpc::BubbleDescription&>(b) = d;
      r.push_back(b);
    }
    return r;
  }();
  // finite element space
  print_memory_footprint("[Building problem ...]");
  auto fed = std::make_shared<mfem_mgis::FiniteElementDiscretization>(
      mfem_mgis::Parameters{
          {"MeshFileName", p.mesh_file},
          {"MeshReadMode", p.parallel_mesh ? "Restart" : "FromScratch"},
          {"FiniteElementFamily", "H1"},
          {"FiniteElementOrder", p.order},
          {"UnknownsSize", 3},
          {"NumberOfUniformRefinements", p.refinement},  // faster for testing
          //{"NumberOfUniformRefinements", parameters.parallel ? 1 : 0},
          {"Parallel", true}});
  // definition of the nonlinear problem
  auto problem = mfem_mgis::PeriodicNonLinearEvolutionProblem{fed};
  print_memory_footprint("[Building problem done]");

  // get problem information
  print_mesh_information(problem.getImplementation<true>());

  // macroscopic strain
  std::vector<mfem_mgis::real> e(6, mfem_mgis::real{0});
  problem.setMacroscopicGradientsEvolution(
      [e](const mfem_mgis::real) { return e; });
  // imposing pressure on the bubble boundaries
  for (const auto& b : bubbles) {
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
  int verbosity = p.verbosity_level;
  int post_processing = p.post_processing;
  mfem_mgis::Parameters solverParameters;
  solverParameters.insert(mfem_mgis::Parameters{{"VerbosityLevel", verbosity}});
  //  solverParameters.insert(mfem_mgis::Parameters{{"MaximumNumberOfIterations",
  //  defaultMaxNumOfIt}});
  //  solverParameters.insert(mfem_mgis::Parameters{{"Tolerance", Tol}});
  //
/*
  auto options = mfem_mgis::Parameters{{"VerbosityLevel", verbosity},
                                       {"Strategy", "Elasticity"}};
  auto preconditioner =
      mfem_mgis::Parameters{{"Name", "HypreBoomerAMG"}, {"Options", options}};
*/
  auto preconditioner =
      mfem_mgis::Parameters{{"Name", "HypreDiagScale"}};

  solverParameters.insert(mfem_mgis::Parameters{
      {"Preconditioner", preconditioner}, {"Tolerance", 1e-9}});
  problem.setLinearSolver("HyprePCG", solverParameters);
  problem.setSolverParameters({{"VerbosityLevel", 1},
                               {"RelativeTolerance", 1e-10},
                               {"AbsoluteTolerance", 0.},
                               {"MaximumNumberOfIterations", 1}});
  // matrix is considered elastic
  problem.addBehaviourIntegrator("Mechanics", 1, p.library, p.behaviour);
  auto& m = problem.getMaterial(1);
  for (auto& s : {&m.s0, &m.s1}) {
    mgis::behaviour::setMaterialProperty(*s, "YoungModulus", 150e9);
    mgis::behaviour::setMaterialProperty(*s, "PoissonRatio", 0.3);
    mgis::behaviour::setExternalStateVariable(*s, "Temperature", 293.15);
    mgis::behaviour::setExternalStateVariable(*s, "Temperature", 293.15);
  }
  if (post_processing) {
    auto results = std::vector<mfem_mgis::Parameter>{"Stress"};
    problem.addPostProcessing(
        "ParaviewExportIntegrationPointResultsAtNodes",
        {{"OutputFileName", p.testcase_name}, {"Results", results}});
  }
  //
  mfem_mgis::size_type nstep{0};
  while ((nstep != maximumNumberSteps) && (!opera_hpc::areAllBroken(bubbles))) {
    CatchTimeSection("BubbleWhileLoop");
    print_memory_footprint("[Step " + std::to_string(nstep) + "]");
    problem.solve(0, 1);
    const auto r = opera_hpc::findFirstPrincipalStressValueAndLocation(
        problem.getMaterial(1));
    const auto max_vp_scaled = p.scale_factor_vp * r.value;
    const auto all_locations_above_threshold = opera_hpc::
    getPointsAboveStressThreshold(
            problem.getMaterial(1), max_vp_scaled);
    // for (auto &location : all_locations_above_threshold){
    //   Message("Here we are");
    //   Message(location[0], "\t",location[1], "\t",location[2], "\t");
    // }
    mfem_mgis::size_type nbroken{0};
    std::vector<mfem_mgis::size_type> all_broken_bubbles_identifiers;
    for (auto &b : bubbles) {
      for (auto &location : all_locations_above_threshold) {
        const auto d = opera_hpc::distance(b, location);
        if (d < dmin) {
          if (!b.broken) {
            ++nbroken;
          }
          b.broken = true;
        }
      }
      if (b.broken) {
        all_broken_bubbles_identifiers.push_back(b.boundary_identifier);
      }
    }

    if (nbroken == 0) {
      Message("no bubble broke at step ", nstep, "\n");
      break;
    }

    Message("step:", nstep);
    Message("-", nbroken, "bubbles broke at this step.");
    Message("-", all_broken_bubbles_identifiers.size(), "bubbles are broken");
    Message("- value of the first principal stress:", r.value,
            "at coordinate (", r.location[0], ",", r.location[1], ",",
            r.location[2], ")");
    
    if (mfem_mgis::getMPIrank() == 0) {
      write_message(output_file, nbroken, "bubbles broke at this step.");
      write_message(output_file, all_broken_bubbles_identifiers.size(),
                    "bubbles are broken");
      write_message(output_file,
                    "Value of the first principal stress:", r.value,
                    "at coordinate (", r.location[0], ",", r.location[1], ",",
                    r.location[2], ")");
    }

    if (post_processing)
      post_process(problem, 0 + double(nstep), 1 + double(nstep));

    ++nstep;
  }
  output_file.close();
  print_memory_footprint("[End]");
  mfem_mgis::Profiler::timers::print_and_write_timers();
  return EXIT_SUCCESS;
}
