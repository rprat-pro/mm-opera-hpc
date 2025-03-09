/*!
 * \file   Bubble.cxx
 * \brief
 * \author T. Helfer, T. Barani, R. Prat
 */

 #include <cstdlib>
 #include <random>
 #include <cmath>
 #include <iostream>
 #include <memory> 
 
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

std::array<mfem_mgis::real, 3u> generate_random_vector() {
  std::random_device rand_dummy;
  std::mt19937 gen(rand_dummy());
  std::uniform_real_distribution<> distribution(-1.0, 1.0);
  return {distribution(gen), distribution(gen), distribution(gen)};
}

std::array<mfem_mgis::real, 3u> generate_plane_by_random_normal() {
  auto v = opera_hpc::generate_random_vector();
  auto square = [](const mfem_mgis::real x) { return x * x; };
  mfem_mgis::real norm = std::sqrt(square(v[0]) + square(v[1]) + square(v[2]));
  if (norm < std::numeric_limits<mfem_mgis::real>::epsilon())
    return {1.0, 0.0, 0.0};
  return {v[0] / norm, v[1] / norm, v[2] / norm};
}

// Function to calculate the distance from a point to a plane defined by a point
// and normal
mfem_mgis::real distance_point_to_plane(
    const std::array<mfem_mgis::real, 3u>& point,
    const std::array<mfem_mgis::real, 3u>& plane_normal,
    const std::array<mfem_mgis::real, 3u>& plane_point) {
  mfem_mgis::real dist =
      std::abs(plane_normal[0] * (point[0] - plane_point[0]) +
               plane_normal[1] * (point[1] - plane_point[1]) +
               plane_normal[2] * (point[2] - plane_point[2]));
  return dist;
}

// Function to check if a bubble (sphere) intersects a plane
bool is_plane_intersecting_bubble(const Bubble& bubble, 
                                  const std::array<mfem_mgis::real, 3u>& plane_normal,
                                  const std::array<mfem_mgis::real, 3u>& plane_point) {
  mfem_mgis::real dist = distance_point_to_plane(bubble.center, plane_normal, plane_point);
  return dist <= bubble.radius;
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
  int post_processing = 1;
  int verbosity_level = 0;  
  mfem_mgis::real scale_factor_vp = 0.9;
  const char* breaking_criterion = "distance"; 
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
  args.AddOption(&p.testcase_name, "-n", "--name-case",
                 "Name of the testcase.");
  args.AddOption(&p.breaking_criterion, "-bc", "--breaking-criterion",
                 "Bubble breaking criterion ('distance' or 'plane').");

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

void write_bubble_positions(
    std::ofstream& file_to_write,
    const std::array<mfem_mgis::real, 3u> bubble_positions) {
  for (auto& el : bubble_positions) file_to_write << el << "\t";
  file_to_write << std::endl;
}

int main(int argc, char** argv) {
  using namespace mfem_mgis::Profiler::Utils;  // Use Message
  // options treatment
  mfem_mgis::initialize(argc, argv);
  mfem_mgis::Profiler::timers::init_timers();

  // file collecting the output
  std::string file_name = "bubbles_and_stresses.txt";
  std::string bubbles_locations_file = "bubbles_positions.txt";
  std::ofstream output_file(file_name);
  std::ofstream output_file_bubbles(bubbles_locations_file);
  if (!output_file.is_open()) {
    std::cerr << "Failed to open the file: " << file_name << std::endl;
    return EXIT_FAILURE;  // Exit the program with an error code
  }

  if (!output_file_bubbles.is_open()) {
    std::cerr << "Failed to open the file: " << bubbles_locations_file
              << std::endl;
    return EXIT_FAILURE;  // Exit the program with an error code
  }

  // get parameters
  TestParameters p;
  mfem::OptionsParser args(argc, argv);
  fill_parameters(args, p);

  // reference pressure
  constexpr mfem_mgis::real pref = 1.0;
  constexpr mfem_mgis::size_type maximumNumberSteps = 1000;
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

  // get problem information
  print_mesh_information(problem.getImplementation<true>());
  print_memory_footprint("[Building problem]");

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
  auto options = mfem_mgis::Parameters{{"VerbosityLevel", verbosity},
                                       {"Strategy", "Elasticity"}};
  auto preconditioner =
      mfem_mgis::Parameters{{"Name", "HypreBoomerAMG"}, {"Options", options}};
  solverParameters.insert(mfem_mgis::Parameters{
      {"Preconditioner", preconditioner}, {"Tolerance", 1e-8}});
  problem.setLinearSolver("HyprePCG", solverParameters);
  problem.setSolverParameters({{"VerbosityLevel", 1},
                               {"RelativeTolerance", 1e-6},
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
  mfem_mgis::size_type nstep{1};
  while ((nstep != maximumNumberSteps) && (!opera_hpc::areAllBroken(bubbles))) {
    problem.solve(0, 1);
    const auto r = opera_hpc::findFirstPrincipalStressValueAndLocation(
        problem.getMaterial(1));
    const auto max_vp_scaled = p.scale_factor_vp * r.value;
    const auto all_locations_above_threshold =
        opera_hpc::getPointsAboveStressThreshold(problem.getMaterial(1),
                                                 max_vp_scaled);
    
    if (p.breaking_criterion == "plane" &&
        !all_locations_above_threshold.empty()) {
      CatchTimeSection("common::plane_creation_and_bubble_ids_step");
      std::vector<opera_hpc::Bubble*> intersected_bubbles;
      mfem_mgis::size_type nbroken_plane{0};
      std::vector<mfem_mgis::size_type> all_broken_bubbles_identifiers_plane;
      std::vector<std::array<mfem_mgis::real, 3u>> all_broken_bubbles_locations_plane;

      //std::array<mfem_mgis::real, 3u> plane_point = all_locations_above_threshold[0]; // we suppose it's the first one but it can change
      std::array<mfem_mgis::real, 3u> plane_point = {r.location[0], r.location[1], r.location[2]}; // we suppose it's the first one but it can change
      std::array<mfem_mgis::real, 3u> plane_normal = opera_hpc::generate_plane_by_random_normal();
      
      //
      Message("Generated random plane at step ", nstep, 
              " at point (", plane_point[0], ",", plane_point[1], ",", plane_point[2],
              ") with normal (", plane_normal[0], ",", plane_normal[1], ",", plane_normal[2], ")");
      //

      for (auto& b : bubbles) {
        if (opera_hpc::is_plane_intersecting_bubble(b, plane_normal, plane_point)) {
          Message("Bubble", b.boundary_identifier, "intersects the plane.");
          intersected_bubbles.push_back(&b);
        }
      }

      for (auto* b_ptr : intersected_bubbles) {
        auto& b = *b_ptr;
        if (!b.broken) {
          ++nbroken_plane;
          all_broken_bubbles_identifiers_plane.push_back(b.boundary_identifier);
          all_broken_bubbles_locations_plane.push_back(b.center);
          b.broken = true;
        }
      }

      if (nbroken_plane > 0) {
        Message("step:", nstep);
        Message("-", nbroken_plane,
                "bubbles broke (plane criterion) at this step.");
        Message("-", all_broken_bubbles_identifiers_plane.size(),
                " bubbles are broken in total.");
        Message("- value of the first principal stress:", r.value,
                "at coordinate (", r.location[0], ",", r.location[1], ",",
                r.location[2], ")");

        if (mfem_mgis::getMPIrank() == 0) {
          write_message(output_file, nbroken_plane,
                        "bubbles broke (plane criterion) at this step.");
          write_message(
              output_file, all_broken_bubbles_identifiers_plane.size(),
              "bubbles are broken in total");
          write_message(output_file,
                        "Value of the first principal stress:", r.value,
                        "at coordinate (", r.location[0], ",", r.location[1],
                        ",", r.location[2], ")");
          for (auto& el : all_broken_bubbles_locations_plane)
            write_bubble_positions(output_file_bubbles, el);
        }
      }

    } else {  // Distance (old) criterion
      mfem_mgis::size_type nbroken_distance{0};
      std::vector<mfem_mgis::size_type> all_broken_bubbles_identifiers_distance;
      std::vector<std::array<mfem_mgis::real, 3u>> all_broken_bubbles_locations_distance;

      for (auto& b : bubbles) {
        for (auto& location : all_locations_above_threshold) {
          const auto d = opera_hpc::distance(b, location);
          if (d < dmin) {
            if (!b.broken) {
              ++nbroken_distance;
              all_broken_bubbles_locations_distance.push_back(location);
            }
            b.broken = true;
          }
        }
        if (b.broken) {
          all_broken_bubbles_identifiers_distance.push_back(
              b.boundary_identifier);
        }
      }

      if (nbroken_distance > 0) {
        Message("step:", nstep);
        Message("-", nbroken_distance,
                "bubbles broke (distance criterion) at this step.");
        Message("-", all_broken_bubbles_identifiers_distance.size(),
                "bubbles are broken");
        Message("- value of the first principal stress:", r.value,
                "at coordinate (", r.location[0], ",", r.location[1], ",",
                r.location[2], ")");

        if (mfem_mgis::getMPIrank() == 0) {
          write_message(output_file, nbroken_distance,
                        "bubbles broke (distance criterion) at this step.");
          write_message(output_file,
                        all_broken_bubbles_identifiers_distance.size(),
                        "bubbles are broken");
          write_message(output_file,
                        "Value of the first principal stress:", r.value,
                        "at coordinate (", r.location[0], ",", r.location[1],
                        ",", r.location[2], ")");
          for (auto& el : all_broken_bubbles_locations_distance)
            write_bubble_positions(output_file_bubbles, el);
        }
      } else if (p.breaking_criterion != "plane") {
        Message("no bubble broke at step ", nstep, " (distance criterion)\n");
      }
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