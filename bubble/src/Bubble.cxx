/*!
 * \file   main.cxx
 * \brief  Main program for bubble fracture simulation in porous materials
 * \date   29/09/2024
 * 
 * This program simulates the mechanical behavior of a material containing bubbles
 * under pressure. It identifies which bubbles are likely to break based on the
 * principal stress field around them.
 */

#include "MFEMMGIS/Profiler.hxx"
#include "Common/Memory.hxx"
#include "Common/Print.hxx"

#include "MGIS/Context.hxx"

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/MechanicalPostProcessings.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/ParaviewExportIntegrationPointResultsAtNodes.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/UniformImposedPressureBoundaryCondition.hxx"
#include "OperaHPC/BubbleDescription.hxx"
#include "OperaHPC/Utilities.hxx"

#include <cstdlib>
#include <fstream>
#include <tuple>

namespace opera_hpc {

  /**
   * \brief Represents a bubble in the material with its state (broken or intact)
   * 
   * Extends BubbleDescription to add a boolean flag indicating whether
   * the bubble has broken during the simulation
   */
  struct Bubble : BubbleDescription {
    bool broken = false;  ///< True if the bubble has fractured
  };

  /**
   * \brief Check if all bubbles in the system have broken
   * 
   * \param bubbles Vector of all bubbles to check
   * \return true if all bubbles are broken, false otherwise
   */
  bool areAllBroken(const std::vector<Bubble>& bubbles) {
    for (const auto& b : bubbles) {
      if (!b.broken) {
        return false;
      }
    }
    return true;
  }
}  // end of namespace opera_hpc

/**
 * \brief Record storing bubble information for output
 * 
 * Contains the boundary ID, location, and maximum stress value
 * found near a bubble during the analysis
 */
struct BubbleInfoRecord {
  mfem_mgis::size_type boundary_id;  ///< Boundary identifier of the bubble
  mfem_mgis::real location[3];       ///< 3D coordinates of max stress location
  mfem_mgis::real stress_value;      ///< Maximum principal stress value

  /**
   * \brief Constructor
   * 
   * \param bid Boundary identifier
   * \param loc 3D location array
   * \param s_val Stress value at that location
   */
  BubbleInfoRecord(mfem_mgis::size_type bid,
                   const std::array<mfem_mgis::real, 3u>& loc,
                   mfem_mgis::real s_val)
      : boundary_id(bid), stress_value(s_val) {
    location[0] = loc[0];
    location[1] = loc[1];
    location[2] = loc[2];
  }
};

/**
 * \brief Execute post-processing operations on the problem
 * 
 * \tparam Problem Type of the mechanical problem
 * \param p The problem instance
 * \param start Start time for post-processing
 * \param end End time for post-processing
 */
template <typename Problem>
void post_process(Problem& p, double start, double end) {
  CatchTimeSection("common::post_processing_step");
  p.executePostProcessings(start, end);
}

/**
 * \brief Structure holding all simulation parameters
 * 
 * Contains default values for mesh files, material properties,
 * numerical parameters, and output options
 */
struct TestParameters {
  const char* mesh_file = "mesh/single_sphere.msh";        ///< Path to mesh file
  const char* behaviour = "Elasticity";                    ///< Material behavior law
  const char* library = "src/libBehaviour.so";             ///< Material library path
  const char* bubble_file = "mesh/single_bubble.txt";      ///< Bubble definitions file
  const char* testcase_name = "TestCaseBubble";            ///< Name for output files
  int parallel_mesh = 0;                                   ///< Flag for parallel mesh format
  int order = 1;                                           ///< Finite element order
  int refinement = 0;                                      ///< Mesh refinement level
  int post_processing = 1;                                 ///< Enable post-processing (1=yes)
  int verbosity_level = 0;                                 ///< Output verbosity level
  mfem_mgis::real scale_factor_vp = 0.9;                   ///< Scale factor for principal stress threshold
};

/**
 * \brief Parse command-line arguments and fill parameter structure
 * 
 * \param args Command-line argument parser
 * \param p Parameter structure to fill
 */
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
  args.AddOption(&p.testcase_name, "-n", "--name-case",
                 "Name of the testcase.");

  args.Parse();

  // Validate arguments
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

/**
 * \brief Write comma-separated values to a file stream
 * 
 * Variadic template function to write any number of arguments
 * separated by commas, ending with a newline
 * 
 * \param file_to_write Output file stream
 * \param args Variadic arguments to write
 */
void write_message(std::ofstream& file_to_write, const auto&... args) {
  ((file_to_write << args << ","), ...) << std::endl;
}

/**
 * \brief Write bubble information records to a CSV file
 * 
 * \param file_to_write Output file stream
 * \param bubble_infos Vector of bubble information records to write
 */
void write_bubble_infos(std::ofstream& file_to_write,
                        const std::vector<BubbleInfoRecord>& bubble_infos) {
  for (auto& entry : bubble_infos) {
    file_to_write << entry.boundary_id << ",";
    for (auto& el : entry.location) file_to_write << el << ",";
    file_to_write << entry.stress_value << "\n";
  }
}

/**
 * \brief Calculate the average hydrostatic pressure over a material domain
 * 
 * Integrates the hydrostatic pressure field over all elements and
 * computes the volume-averaged value using MPI reduction if needed
 * 
 * \tparam parallel Whether running in parallel mode
 * \param os Output stream for results
 * \param prob The nonlinear evolution problem
 * \param f Quadrature function containing the pressure field
 */
template <bool parallel>
static void calculateAverageHydroPressure(
    std::ostream& os,
    mfem_mgis::NonLinearEvolutionProblemImplementation<parallel>& prob,
    const mfem_mgis::ImmutablePartialQuadratureFunctionView& f) {
  const auto& s = f.getPartialQuadratureSpace();
  const auto& fed = s.getFiniteElementDiscretization();
  const auto& fespace = fed.getFiniteElementSpace<parallel>();
  const auto m = s.getId();
  
  // Storage for local contributions
  std::vector<mfem_mgis::real> integrals;
  std::vector<mfem_mgis::real> volumes;
  mfem_mgis::real integr, vol;
  
  // Loop over all elements in the mesh
  for (mfem_mgis::size_type i = 0; i != fespace.GetNE(); ++i) {
    // Skip elements not in this material
    if (fespace.GetAttribute(i) != m) {
      continue;
    }
    const auto& fe = *(fespace.GetFE(i));
    auto& tr = *(fespace.GetElementTransformation(i));
    const auto& ir = s.getIntegrationRule(fe, tr);
    integr = 0.;
    vol = 0.;
    
    // Integrate over quadrature points
    for (mfem_mgis::size_type g = 0; g != ir.GetNPoints(); ++g) {
      mfem::Vector p;
      const auto& ip = ir.IntPoint(g);
      tr.SetIntPoint(&ip);
      const auto w =
          prob.getBehaviourIntegrator(m).getIntegrationPointWeight(tr, ip);
      tr.Transform(tr.GetIntPoint(), p);

      integr += f.getIntegrationPointValue(i, g) * w;
      vol += w;
    }
    integrals.push_back(integr);
    volumes.push_back(vol);
  }

  // MPI reduction to get global integral and volume
  mfem_mgis::real v = 0.;
  mfem_mgis::real integral = 0.;
  MPI_Reduce(integrals.data(), &integral, 1, MPI_DOUBLE, MPI_SUM, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(volumes.data(), &v, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  
  // Write results on rank 0
  if (mfem_mgis::getMPIrank() == 0) {
    os << v << "," << integral << "," << integral / v << "\n";
  }
}

/**
 * \brief Main program
 * 
 * Workflow:
 * 1. Initialize stuff and parse command-line arguments
 * 2. Load mesh and bubble definitions
 * 3. Set up the finite element problem
 * 4. Apply pressure boundary conditions on bubbles
 * 5. Solve the mechanical equilibrium
 * 6. Find maximum principal stress locations near each bubble
 * 7. Identify bubbles based on stress threshold
 * 8. Export results and post-process data
 */

int main(int argc, char** argv) {
  using namespace mfem_mgis::Profiler::Utils;
  
  // Initialize MPI and profiling
  mfem_mgis::initialize(argc, argv);
  mfem_mgis::Profiler::timers::init_timers();
  auto ctx = mgis::Context();

  // Open output file for bubble stress results
  std::string bubbles_selected = "bubbles_and_stresses_selected.txt";
  std::ofstream output_file_locations(bubbles_selected);

  if (!output_file_locations.is_open()) {
    std::cerr << "Failed to open the file: " << bubbles_selected << std::endl;
    return EXIT_FAILURE;
  }

  // Parse command-line parameters
  TestParameters p;
  mfem::OptionsParser args(argc, argv);
  fill_parameters(args, p);

  print_memory_footprint("[Start]");

  // Physical parameters
  constexpr mfem_mgis::real pref = 1.0;  ///< Reference pressure in bubbles
  constexpr mfem_mgis::size_type maximumNumberSteps = 1000;  ///< Max time steps
  constexpr mfem_mgis::real dmin = 0.650;  ///< Distance threshold to associate stress with bubble
  
  // Load bubble descriptions from file
  auto bubbles = [p] {
    auto r = std::vector<opera_hpc::Bubble>{};
    for (const auto& d : opera_hpc::BubbleDescription::read(p.bubble_file)) {
      auto b = opera_hpc::Bubble{};
      static_cast<opera_hpc::BubbleDescription&>(b) = d;
      r.push_back(b);
    }
    return r;
  }();

  // Create finite element discretization
  print_memory_footprint("[Building problem ...]");
  auto fed = std::make_shared<mfem_mgis::FiniteElementDiscretization>(
      mfem_mgis::Parameters{
          {"MeshFileName", p.mesh_file},
          {"MeshReadMode", p.parallel_mesh ? "Restart" : "FromScratch"},
          {"FiniteElementFamily", "H1"},  // Continuous Lagrange elements
          {"FiniteElementOrder", p.order},
          {"UnknownsSize", 3},  // 3D displacement field
          {"NumberOfUniformRefinements", p.refinement},
          {"Parallel", true}});

  // Define the nonlinear evolution problem
  auto problem = mfem_mgis::PeriodicNonLinearEvolutionProblem{fed};
  print_memory_footprint("[Building problem done]");

  print_mesh_information(problem.getImplementation<true>());

  // Set macroscopic strain to zero
  std::vector<mfem_mgis::real> e(6, mfem_mgis::real{0});
  problem.setMacroscopicGradientsEvolution(
      [e](const mfem_mgis::real) { return e; });
  
  // Apply pressure boundary conditions on each bubble
  for (const auto& b : bubbles) {
    problem.addBoundaryCondition(
        std::make_unique<mfem_mgis::UniformImposedPressureBoundaryCondition>(
            problem.getFiniteElementDiscretizationPointer(),
            b.boundary_identifier, [&b](const mfem_mgis::real) {
              // If bubble is broken, pressure = 0; otherwise pressure = pref
              return b.broken ? 0 : pref;
            }));
  }

  // Configure linear solver
  int verbosity = p.verbosity_level;
  int post_processing = p.post_processing;
  mfem_mgis::Parameters solverParameters;
  solverParameters.insert(mfem_mgis::Parameters{{"VerbosityLevel", verbosity}});
  
  // Use diagonal scaling preconditioner
  auto preconditioner = mfem_mgis::Parameters{{"Name", "HypreDiagScale"}};
  solverParameters.insert(mfem_mgis::Parameters{
      {"Preconditioner", preconditioner}, {"Tolerance", 1e-10}});

  // Set up conjugate gradient solver
  problem.setLinearSolver("HyprePCG", solverParameters);
  problem.setSolverParameters({{"VerbosityLevel", 1},
                               {"RelativeTolerance", 1e-11},
                               {"AbsoluteTolerance", 0.},
                               {"MaximumNumberOfIterations", 1}});

  // Add elastic material behavior
  problem.addBehaviourIntegrator("Mechanics", 1, p.library, p.behaviour);
  auto& m = problem.getMaterial(1);
  
  // Set material properties (elastic constants and temperature)
  for (auto& s : {&m.s0, &m.s1}) {
    mgis::behaviour::setMaterialProperty(*s, "YoungModulus", 150e-3);  // Young's modulus
    mgis::behaviour::setMaterialProperty(*s, "PoissonRatio", 0.3);     // Poisson's ratio
    mgis::behaviour::setExternalStateVariable(*s, "Temperature", 293.15);  // Room temperature
  }

  // Configure post-processing output
  if (post_processing) {
    auto results =
        std::vector<mfem_mgis::Parameter>{"Stress", "HydrostaticPressure"};
    problem.addPostProcessing(
        "ParaviewExportIntegrationPointResultsAtNodes",
        {{"OutputFileName", p.testcase_name}, {"Results", results}});
  }

  // Storage for internal variables
  auto vim = mfem_mgis::PartialQuadratureFunction{
      m.getPartialQuadratureSpacePointer(), 1};

  mfem_mgis::size_type nstep{1};

  // Solve the mechanical equilibrium problem
  problem.solve(0, 1);

  // Find the maximum principal stress in the domain
  const auto r = opera_hpc::findFirstPrincipalStressValueAndLocation(
      problem.getMaterial(1));

  std::vector<BubbleInfoRecord> bubbles_information;

  // Define stress threshold (90% of maximum principal stress by default)
  const auto max_vp_scaled = p.scale_factor_vp * r.value;
  
  // Get all locations where stress exceeds threshold
  const auto all_locations_and_stresses_above_threshold =
      opera_hpc::getPointsandStressAboveStressThreshold(problem.getMaterial(1),
                                                        max_vp_scaled);

  // For each bubble, find the maximum stress within distance dmin
  for (auto& b : bubbles) {
    CatchTimeSection("BubbleLoop");
    Message("Treating bubble #", b.boundary_identifier);
    mfem_mgis::real tmp_stress = 0.;
    std::array<mfem_mgis::real, 3u> tmp_loc{0., 0., 0.};
    
    // Check all high-stress locations
    for (auto& location_and_stress :
         all_locations_and_stresses_above_threshold) {
      const auto d = opera_hpc::distance(b, location_and_stress.location);
      
      // If location is close enough and stress is higher than current max
      if ((d < dmin) && (tmp_stress < location_and_stress.value)) {
        tmp_stress = location_and_stress.value;
        tmp_loc = location_and_stress.location;
      }
    }
    // Store the result for this bubble
    bubbles_information.emplace_back(b.boundary_identifier, tmp_loc,
                                     tmp_stress);
  }

  // Write CSV header and bubble stress data
  write_message(output_file_locations, "Bubble ID", "Location[0]",
                "Location[1]", "Location[2]", "Stress");

  for (auto& el : bubbles_information) {
    Message("Bubble = ", el.boundary_id);
    Message("Stress = ", el.stress_value);
    for (auto& el1 : el.location) {
      Message("Location=", el1);
    }
  }

  write_bubble_infos(output_file_locations, bubbles_information);
  
  // Open file for hydrostatic stress output
  std::ofstream outfile("sig_hydro.txt");

  if (mfem_mgis::getMPIrank() == 0)
    write_message(outfile, "Volume", "Integral", "Average");

  if (post_processing) {
    // Execute standard post-processing (Paraview export)
    post_process(problem, 0 + double(nstep), 1 + double(nstep));

    // Compute and export principal stress field
    auto end = mfem_mgis::Material::END_OF_TIME_STEP;
    mfem_mgis::computeFirstEigenStress(ctx, vim, m, end);
    mfem_mgis::ParaviewExportIntegrationPointResultsAtNodesImplementation<true>
        export_stress(problem.getImplementation<true>(),
                      {{.name = "Stress", .functions = {vim}}}, 
                      std::string(p.testcase_name) + "-eig");
    
    export_stress.execute(problem.getImplementation<true>(), 0 + double(nstep),
                          1 + double(nstep));

    // Extract and compute average hydrostatic pressure
    auto pr = getInternalStateVariable(
        static_cast<const mfem_mgis::Material&>(m), "HydrostaticPressure");

    calculateAverageHydroPressure<true>(outfile,
                                        problem.getImplementation<true>(), pr);
  }

  // Clean up, close and win
  outfile.close();
  output_file_locations.close();
  print_memory_footprint("[End]");
  mfem_mgis::Profiler::timers::print_and_write_timers();
  return EXIT_SUCCESS;
}