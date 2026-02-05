/*!
 * \file   main.cxx
 * \brief  Main program for bubble fracture simulation in porous materials
 * \date   29/09/2024
 *
 * This program simulates the mechanical behavior of a material containing
 * bubbles under pressure. It identifies which bubbles are likely to break based
 * on the principal stress field around them.
 */


#include <cstdlib>
#include <fstream>
#include <tuple>

#include "MGIS/Function/TFEL/Mechanics.hxx"

#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/MechanicalPostProcessings.hxx"
#include "MFEMMGIS/ParaviewExportIntegrationPointResultsAtNodes.hxx"
#include "MFEMMGIS/PartialQuadratureFunction.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/UniformImposedPressureBoundaryCondition.hxx"

#include "Common/Memory.hxx"
#include "Common/Print.hxx"
#include "OperaHPC/BubbleDescription.hxx"
#include "OperaHPC/Utilities.hxx"

#include "mfem/general/forall.hpp"

namespace mfem_mgis {
  template <unsigned short N>
  static bool computeHydrostaticPressure_implementation(
      Context& ctx,
      PartialQuadratureFunction& hyp,
      const Material& m,
      const Material::StateSelection s) {
    using namespace mgis::function;
    const auto sig = getThermodynamicForce(m, "Stress", s);
    const auto ok = sig | as_stensor<N> | hydrostatic_stress | hyp;
    if (!ok) {
      return ctx.registerErrorMessage(
          "computeHydrostaticPressure: computation of the hydrostatic "
          "stress failed");
    }
    return true;
  }  // end of computeHydrostaticPressure_implementation

  static bool computeHydrostaticPressure(Context& ctx,
                                         PartialQuadratureFunction& hyp,
                                         const Material& m,
                                         const Material::StateSelection s) {
    if ((m.b.hypothesis == Hypothesis::PLANESTRESS) ||
        (m.b.hypothesis == Hypothesis::PLANESTRAIN)) {
      return computeHydrostaticPressure_implementation<2>(ctx, hyp, m, s);
    } else if (m.b.hypothesis == Hypothesis::TRIDIMENSIONAL) {
      return computeHydrostaticPressure_implementation<3>(ctx, hyp, m, s);
    }
    return ctx.registerErrorMessage(
        "computeHydrostaticPressure: unsupported modelling hypothesis");
  }  // end of computeHydrostaticPressure

}  // namespace mfem_mgis

namespace opera_hpc {

  /**
   * \brief Represents a bubble in the material with its state (broken or
   * intact)
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
  BubbleInfoRecord(const mfem_mgis::size_type bid,
                   const std::array<mfem_mgis::real, 3u>& loc,
                   const mfem_mgis::real s_val)
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
  p.executePostProcessings(start, end);
}

/**
 * \brief Physical parameters for bubble fracture simulation
 *
 * Contains the key physical parameters that control the mechanical
 * behavior of bubbles and the stress analysis criteria used to
 * determine bubble fracture.
 */
struct PhysicalParameters {
  mfem_mgis::real bubble_pressure = 1.0;  ///< Reference pressure inside bubbles
  mfem_mgis::real d_min =
      0.950;  ///< Distance threshold for associating stress with bubbles
};

/**
 * \brief Structure holding all simulation parameters
 *
 * Contains default values for mesh files, material properties,
 * numerical parameters, and output options
 */
struct TestParameters : PhysicalParameters {
  const char* mesh_file = "mesh/single_sphere.msh";  ///< Path to mesh file
  const char* behaviour = "Elasticity";              ///< Material behavior law
  const char* library = "src/libBehaviour.so";       ///< Material library path
  const char* bubble_file =
      "mesh/single_bubble.txt";                  ///< Bubble definitions file
  const char* testcase_name = "TestCaseBubble";  ///< Name for output files
  int parallel_mesh = 0;    ///< Flag for parallel mesh format
  int order = 2;            ///< Finite element order
  int refinement = 0;       ///< Mesh refinement level
  int post_processing = 1;  ///< Enable post-processing (1=yes)
  int verbosity_level = 0;  ///< Output verbosity level
  mfem_mgis::real scale_factor_vp =
      0.9;  ///< Scale factor for principal stress threshold
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
  args.AddOption(&p.d_min, "-dm", "--dmin",
                 "Distance to evaluate bubble rupture, in µm, default = 0.950.");
  args.AddOption(&p.bubble_pressure, "-pr", "--pref",
                 "Internal bubble pressure, in Pa, default = 1");

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
 * \brief Calculate the average hydrostatic stress over a material domain
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
static void calculateAverageHydrostaticStress(
    std::ostream& os,
    const mfem_mgis::NonLinearEvolutionProblemImplementation<parallel>& prob,
    const mfem_mgis::ImmutablePartialQuadratureFunctionView& f) {
  const auto& s = f.getPartialQuadratureSpace();
  const auto& fed = s.getFiniteElementDiscretization();
  const auto& fespace = fed.getFiniteElementSpace<parallel>();
  const auto m = s.getId();

  mfem_mgis::real integral = 0.;
  mfem_mgis::real volume = 0.;

  // Loop over all elements in the mesh
  for (mfem_mgis::size_type i = 0; i != fespace.GetNE(); ++i) {
    // Skip elements not in this material
    if (fespace.GetAttribute(i) != m) {
      continue;
    }
    const auto& fe = *(fespace.GetFE(i));
    auto& tr = *(fespace.GetElementTransformation(i));
    const auto& ir = s.getIntegrationRule(fe, tr);

    // Integrate over quadrature points
    for (mfem_mgis::size_type g = 0; g != ir.GetNPoints(); ++g) {
      const auto& ip = ir.IntPoint(g);
      tr.SetIntPoint(&ip);
      const auto w =
          prob.getBehaviourIntegrator(m).getIntegrationPointWeight(tr, ip);

      integral += f.getIntegrationPointValue(i, g) * w;
      volume += w;
    }
  }

  // MPI reduction to get global integral and volume
  mfem_mgis::real global_volume = 0.;
  mfem_mgis::real global_integral = 0.;
  MPI_Reduce(&integral, &global_integral, 1, MPI_DOUBLE, MPI_SUM, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(&volume, &global_volume, 1, MPI_DOUBLE, MPI_SUM, 0,
             MPI_COMM_WORLD);

  // Write results on rank 0
  if (mfem_mgis::getMPIrank() == 0) {
    os << global_volume << "," << global_integral << ","
       << global_integral / global_volume << "\n";
  }
}  // end of calculateAverageHydrostaticStress

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
  auto ctx = mfem_mgis::Context();

  // Open output file for bubble stress results
  std::string bubbles_selected = "bubbles_and_stresses_selected.txt";
  std::string hydrostatic_stress_f = "sig_hydro.txt";
  std::ofstream outfile_sighydr(hydrostatic_stress_f);
  std::ofstream output_file_locations(bubbles_selected);

  if (!output_file_locations.is_open()) {
    std::cerr << "Failed to open the file: " << bubbles_selected << std::endl;
    return EXIT_FAILURE;
  }

  if (!outfile_sighydr.is_open()) {
    std::cerr << "Failed to open the file: " << hydrostatic_stress_f
              << std::endl;
    return EXIT_FAILURE;
  }

  // Parse command-line parameters
  TestParameters p;
  mfem::OptionsParser args(argc, argv);
  fill_parameters(args, p);

  print_memory_footprint("[Start]");

  // Physical parameters
  mfem_mgis::real dmin = p.d_min;

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
            b.boundary_identifier,
            [&b, pref = p.bubble_pressure](const mfem_mgis::real) {
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
  // The sphere is not meshed, thus no material properties need to be
  // considered for it. We set the elastic constants for the material
  // which are representative for UO2.
  // NB we are considering an isotropic purely elastic behaviour,
  // imposing only the two elastic moduli. Please note that since
  // we are scaling everything to µm, we consider the Elastic
  // modulus in N/µm² to have a consistent displacement field.
  for (auto& s : {&m.s0, &m.s1}) {
    mgis::behaviour::setMaterialProperty(*s, "YoungModulus",
                                         150e-3);  // Young's modulus
    mgis::behaviour::setMaterialProperty(*s, "PoissonRatio",
                                         0.3);  // Poisson's ratio
    mgis::behaviour::setExternalStateVariable(*s, "Temperature",
                                              293.15);  // Room temperature
  }

  // Configure post-processing output
  if (post_processing) {
    auto results = std::vector<mfem_mgis::Parameter>{"Stress"};
    problem.addPostProcessing(
        "ParaviewExportIntegrationPointResultsAtNodes",
        {{"OutputFileName", p.testcase_name}, {"Results", results}});
  }

  // Partial Quadrature functions for post-processing needs
  // eig for the first eigenstress
  auto eig = mfem_mgis::PartialQuadratureFunction{
      m.getPartialQuadratureSpacePointer(), 1};
  // hyp for the hydrostatic stress
  auto hyp = mfem_mgis::PartialQuadratureFunction{
      m.getPartialQuadratureSpacePointer(), 1};

  mfem_mgis::size_type nstep{1};

  // Solve the mechanical equilibrium problem
  problem.solve(0, 1);

  // Find the maximum principal stress in the domain
  const auto r = opera_hpc::findFirstPrincipalStressValueAndLocation(
      problem.getMaterial(1));

  std::vector<BubbleInfoRecord> bubbles_information;

  // Define stress threshold
  const auto max_vp_scaled = p.scale_factor_vp * r.value;

  // Get all locations where stress exceeds threshold
  auto all_locations_and_stresses_above_threshold =
      opera_hpc::getPointsandStressAboveStressThreshold(problem.getMaterial(1),
                                                        max_vp_scaled);
  // Resize the bubble information vector to match the number of bubbles,
  // initializing each entry with default values (0 identifier, {0,0,0}
  // location, 0.0 stress)
  bubbles_information.resize(bubbles.size(),
                             BubbleInfoRecord(0, {0, 0, 0}, 0.0));

  {
    CatchTimeSection("BubbleLoop");

    // Sort all  points by their X-coordinate for efficient spatial
    // searching
    std::sort(all_locations_and_stresses_above_threshold.begin(),
              all_locations_and_stresses_above_threshold.end(),
              [](const auto& a, const auto& b) {
                return a.location[0] < b.location[0];
              });

    // Pointers to avoid repeated vector access
    const auto* bubbles_data = bubbles.data();
    const auto* points_data = all_locations_and_stresses_above_threshold.data();
    auto* results_data = bubbles_information.data();

    const auto num_bubbles = bubbles.size();
    const auto num_points = all_locations_and_stresses_above_threshold.size();

    // Parallel loop over all bubbles
    // NB it can run OMP on CPU parallelization if hypre was compiled with this
    // support otherwise it degenerates to a simple for. it can do something
    // exotic with gpu but idk.
    mfem::forall(num_bubbles, [=] MFEM_HOST_DEVICE (int i) {
      // 1. Get current bubble to process
      const auto& b = bubbles_data[i];

      // Initialize tracking variables for maximum stress point near this bubble
      auto local_max_stress = mfem_mgis::real{0.0};
      auto local_max_loc = std::array<mfem_mgis::real, 3u>{0.0, 0.0, 0.0};

      // bubble center coordinates
      const auto bx = b.center[0];
      const auto by = b.center[1];
      const auto bz = b.center[2];

      // X-axis slice:
      // Only search points within distance 'dmin' along the X-axis
      const auto min_x_search = bx - dmin;
      const auto max_x_search = bx + dmin;

      // Find the starting point for our search using binary search on sorted
      // X-coordinates
      const auto* start_ptr = points_data;
      const auto* end_ptr = points_data + num_points;

      // Find first point with X-coordinate >= min_x_search
      const auto* it = std::lower_bound(start_ptr, end_ptr, min_x_search,
                                        [](const auto& point, double val) {
                                          return point.location[0] < val;
                                        });

      // Scan through candidate points, using early exits and filtering
      // Stop when we leave the X-slice
      for (; it != end_ptr; ++it) {
        const auto& pt = *it;

        // Early exit: If X-coordinate exceeds our search window, stop
        // immediately (points are sorted by X, so all remaining points are also
        // too far --> no need to do the rest of the stuff, we
        // dont care about the stress there!)
        if (pt.location[0] > max_x_search) break;
        // Now we carry on with a matrioska checking
        // Second check is on the y: check if point is within bounding box on
        // Y-axis
        if (std::abs(pt.location[1] - by) > dmin) continue;

        // Third check is on the z: check if point is within bounding box on
        // 2-axis
        if (std::abs(pt.location[2] - bz) > dmin) continue;
        // if we made here, the point is worthwile of consideration!
        const auto dx = pt.location[0] - bx;
        const auto dy = pt.location[1] - by;
        const auto dz = pt.location[2] - bz;

        const auto dist = [&dx, &dy, &dz] {
          return std::sqrt(tfel::math::power<2>(dx) + tfel::math::power<2>(dy) +
                           tfel::math::power<2>(dz));
        }();

        // If point is within spherical radius 'dmin' and has higher stress than
        // current max, update the maximum stress and its location
        if ((dist < dmin) && (local_max_stress < pt.value)) {
          local_max_stress = pt.value;
          local_max_loc = pt.location;
        }
      }

      // Store the maximum stress information found for this bubble
      results_data[i] = BubbleInfoRecord(b.boundary_identifier, local_max_loc,
                                         local_max_stress);
    });
  }

  if (mfem_mgis::getMPIrank() == 0) {
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
  }

  if (mfem_mgis::getMPIrank() == 0)
    write_message(outfile_sighydr, "Volume", "Integral", "Average");

  if (post_processing) {
    CatchTimeSection("common::post_processing_step");
    // Execute standard post-processing (Paraview export)
    post_process(problem, 0, 1);

    // Compute and export principal stress and hydrostatic stress
    // field
    auto end = mfem_mgis::Material::END_OF_TIME_STEP;
    mfem_mgis::computeFirstEigenStress(ctx, eig, m, end);
    mfem_mgis::computeHydrostaticPressure(ctx, hyp, m, end);

    // calculate the average hydrostatic stress on the rev
    calculateAverageHydrostaticStress<true>(
        outfile_sighydr, problem.getImplementation<true>(), hyp);

    mfem_mgis::ParaviewExportIntegrationPointResultsAtNodesImplementation<true>
      export_stress(ctx, problem.getImplementation<true>(),
                      {{.name = "FirstEigenStress", .functions = {eig}},
                       {.name = "HydrostaticPressure", .functions = {hyp}}},
                      std::string(p.testcase_name) + "-pp");

    export_stress.execute(problem.getImplementation<true>(), 0, 1);
  }

  // Clean up, close, and win
  if (mfem_mgis::getMPIrank() == 0) {
    outfile_sighydr.close();
    output_file_locations.close();
  }

  print_memory_footprint("[End]");
  mfem_mgis::Profiler::timers::print_and_write_timers();
  return EXIT_SUCCESS;
}
