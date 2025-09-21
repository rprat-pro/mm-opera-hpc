/*!
 * \file   Mono_UO2_CosH_Jaco.cxx
 * \brief
 * This example is modelling ...
 * \author Raphael Prat, Maxence Wangermez
 * \date   06/2023
 */

#include <memory>
#include <csignal>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <functional>
#include <string_view>

#include "MFEMMGIS/AnalyticalTests.hxx"
#include "MFEMMGIS/Config.hxx"
#include "MFEMMGIS/ImposedDirichletBoundaryConditionAtClosestNode.hxx"
#include "MFEMMGIS/MFEMForward.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/Parameters.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/MechanicalPostProcessings.hxx"
#include "MFEMMGIS/ParaviewExportIntegrationPointResultsAtNodes.hxx"
#include "MFEMMGIS/PartialQuadratureFunctionsSet.hxx"
#include "MM_OPERA_HPC/GrainOrientations.hxx"
#include "MM_OPERA_HPC/MacroscropicElasticMaterialProperties.hxx"
#include "MM_OPERA_HPC/UniaxialMacroscopicStressPeriodicSimulation.hxx"

/*
                Problem :
                Parameters :
*/

// We need this class for test case sources
struct TestParameters {
  const char *output_file = "uniaxial-polycrystal.res";
  const char *mesh_file = "mesh/5cristals.msh";
  const char *vect_file = "mesh/vectors_5cristals.txt";
  //  const char *behaviour = "MonoCristal_UO2";
  const char *behaviour = "Mono_UO2_Cosh_Jaco3";
  const char *library = "src/libBehaviour.so";
  int order = 1;
  bool parallel = true;
  bool export_von_Mises_stress = false;
  bool export_first_eigen_stress = false;
  int refinement = 0;
  int post_processing = 1; // default value : disabled
  int verbosity_level = 0; // default value : lower level
  mfem_mgis::real duration = 200;
  int nstep = 200;
};

static void common_parameters(mfem::OptionsParser &args, TestParameters &p) {
  args.AddOption(&p.output_file, "", "--macroscopic-stress-output-file",
                 "main output file containing the evolution of the diagonal "
                 "components of the deformation gradient and the  diagonal "
                 "components of the Cauchy stress.");
  args.AddOption(&p.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&p.vect_file, "-f", "--vect", "Vector file to use.");
  args.AddOption(&p.library, "-l", "--library", "Material library.");
  args.AddOption(&p.order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&p.refinement, "-r", "--refinement",
                 "refinement level of the mesh, default = 0");
  args.AddOption(&p.export_von_Mises_stress, "",
                 "--enable-export-von-Mises-stress", "",
                 "--disable-export-von-Mises-stress",
                 "export the von Mises stress", false);
  args.AddOption(&p.export_first_eigen_stress, "",
                 "--enable-export-first_eigen_stress", "",
                 "--disable-export-first_eigen_stress",
                 "export first eigen stress", false);
  args.AddOption(&p.post_processing, "-p", "--post-processing",
                 "run post processing step");
  args.AddOption(&p.verbosity_level, "-v", "--verbosity-level",
                 "choose the verbosity level");
  args.AddOption(&p.duration, "-d", "--duration",
                 "choose the duration (default = 5)");
  args.AddOption(&p.nstep, "-n", "--nstep",
                 "choose the number of steps (default = 40)");

  args.Parse();

  if (!args.Good()) {
    if (mfem_mgis::getMPIrank() == 0)
      args.PrintUsage(std::cout);
    mfem_mgis::finalize();
    exit(0);
  }
  if (p.mesh_file == nullptr) {
    if (mfem_mgis::getMPIrank() == 0)
      std::cout << "ERROR: Mesh file missing" << std::endl;
    args.PrintUsage(std::cout);
    mfem_mgis::abort(EXIT_FAILURE);
  }
  if (mfem_mgis::getMPIrank() == 0)
    args.PrintOptions(std::cout);
}

static void
print_mesh_information(mfem_mgis::PeriodicNonLinearEvolutionProblem &problem) {
  CatchTimeSection("common::print_mesh_information");

  using mfem_mgis::Profiler::Utils::Message;
  using mfem_mgis::Profiler::Utils::sum;

  auto &fespace = problem.getImplementation<true>().getFiniteElementSpace();
  // getMesh
  auto mesh = fespace.GetMesh();

  // get the number of vertices
  int64_t numbers_of_vertices_local = mesh->GetNV();
  int64_t numbers_of_vertices = sum(numbers_of_vertices_local);

  // get the number of elements
  int64_t numbers_of_elements_local = mesh->GetNE();
  int64_t numbers_of_elements = sum(numbers_of_elements_local);

  // get the element size
  const mfem_mgis::real h = mesh->GetElementSize(0);

  // get n dofs
  int64_t unknowns_local = fespace.GetTrueVSize();
  int64_t unknowns = sum(unknowns_local);

  Message("Info problem: number of vertices -> ", numbers_of_vertices);
  Message("Info problem: number of elements -> ", numbers_of_elements);
  Message("Info prolbem: element size -> ", h);
  Message("Info porblem: Number of finite element unknowns: ", unknowns);
}

static void add_post_processings(
    mfem_mgis::PeriodicNonLinearEvolutionProblem &p,
    const TestParameters &params,
    const std::string &msg) {
  p.addPostProcessing("ParaviewExportResults", {{"OutputFileName", msg}});
  p.addPostProcessing("MeanThermodynamicForces",
                      {{"OutputFileName", "avgStressPolycristal"}});
#ifdef MGIS_FUNCTION_SUPPORT
  if (params.export_von_Mises_stress) {
    p.getImplementation<true>().addPostProcessing(
        std::make_unique<
            mfem_mgis::
                ParaviewExportIntegrationPointPostProcessingsResultsAtNodes<
                    true>>(
            p.getImplementation<true>(), "vonMisesStress",
            p.getAssignedMaterialsIdentifiers(), 1,
            [&p](mfem_mgis::Context &ctx,
                 mfem_mgis::PartialQuadratureFunction &f) {
              const auto mid = f.getPartialQuadratureSpace().getId();
              const auto &m = p.getBehaviourIntegrator(mid).getMaterial();
              return mfem_mgis::computeVonMisesEquivalentStress(
                  ctx, f, m, mfem_mgis::Material::END_OF_TIME_STEP);
            },
            "vonMisesStressOutput"));
  }
  if (params.export_first_eigen_stress) {
    p.getImplementation<true>().addPostProcessing(
        std::make_unique<
            mfem_mgis::
                ParaviewExportIntegrationPointPostProcessingsResultsAtNodes<
                    true>>(
            p.getImplementation<true>(), "FirstEigenStress",
            p.getAssignedMaterialsIdentifiers(), 1,
            [&p](mfem_mgis::Context &ctx,
                 mfem_mgis::PartialQuadratureFunction &f) {
              const auto mid = f.getPartialQuadratureSpace().getId();
              const auto &m = p.getBehaviourIntegrator(mid).getMaterial();
              return mfem_mgis::computeFirstEigenStress(
                  ctx, f, m, mfem_mgis::Material::END_OF_TIME_STEP);
            },
            "FirstEigenStressOutput"));
  }
#endif
  // p.addPostProcessing(
  // 		"ParaviewExportIntegrationPointResultsAtNodes",
  // 		{{"OutputFileName", msg + "IntegrationPointOutputPKI"},
  // 		 {"Results", {"FirstPiolaKirchhoffStress"}}});
  // p.addPostProcessing(
  // 		"ParaviewExportIntegrationPointResultsAtNodes",
  // 		{{"OutputFileName", msg + "IntegrationPointOutputDG"},
  // 		 {"Results", {"DeformationGradient"}}});
}  // end timer add_postprocessing_and_outputs

static void
setup_properties(mfem_mgis::PeriodicNonLinearEvolutionProblem &problem,
                 mm_opera_hpc::MacroscropicElasticMaterialProperties &mp,
                 const TestParameters &p) {
  using namespace mgis::behaviour;
  using real = mfem_mgis::real;

  CatchTimeSection("set_mgis_stuff");

  // const int nMat = 8;
  const int nMat =
      getMaterialsAttributes(*(problem.getFiniteElementDiscretizationPointer()))
          .Max();
  mfem_mgis::Profiler::Utils::Message("Nombre de matériaux : ", nMat);

  for (int i = 0; i < nMat; i++) {
    problem.addBehaviourIntegrator("Mechanics", i + 1, p.library, p.behaviour);
  }

  // cubic symmetry elasticity
  const mfem_mgis::real young1 = 222.e9;
  const mfem_mgis::real young2 = young1;
  const mfem_mgis::real young3 = young1;
  const mfem_mgis::real poisson12 = 0.27;
  const mfem_mgis::real poisson23 = poisson12;
  const mfem_mgis::real poisson13 = poisson12;
  const mfem_mgis::real shear12 = 54.e9;
  const mfem_mgis::real shear23 = shear12;
  const mfem_mgis::real shear13 = shear12;

  //
  mp.update(young1, poisson12, shear12);
  // materials
  auto set_properties =
      [](auto &m, const mfem_mgis::real yo1, const mfem_mgis::real yo2,
         const mfem_mgis::real yo3, const mfem_mgis::real po12,
         const mfem_mgis::real po23, const mfem_mgis::real po13,
         const mfem_mgis::real sm12, const mfem_mgis::real sm23,
         const mfem_mgis::real sm13) {
        setMaterialProperty(m.s0, "YoungModulus1", yo1);
        setMaterialProperty(m.s0, "YoungModulus2", yo2);
        setMaterialProperty(m.s0, "YoungModulus3", yo3);
        setMaterialProperty(m.s0, "PoissonRatio12", po12);
        setMaterialProperty(m.s0, "PoissonRatio23", po23);
        setMaterialProperty(m.s0, "PoissonRatio13", po13);
        setMaterialProperty(m.s0, "ShearModulus12", sm12);
        setMaterialProperty(m.s0, "ShearModulus23", sm23);
        setMaterialProperty(m.s0, "ShearModulus13", sm13);

        setMaterialProperty(m.s1, "YoungModulus1", yo1);
        setMaterialProperty(m.s1, "YoungModulus2", yo2);
        setMaterialProperty(m.s1, "YoungModulus3", yo3);
        setMaterialProperty(m.s1, "PoissonRatio12", po12);
        setMaterialProperty(m.s1, "PoissonRatio23", po23);
        setMaterialProperty(m.s1, "PoissonRatio13", po13);
        setMaterialProperty(m.s1, "ShearModulus12", sm12);
        setMaterialProperty(m.s1, "ShearModulus23", sm23);
        setMaterialProperty(m.s1, "ShearModulus13", sm13);
      };

  auto set_temperature = [](auto &m) {
    setExternalStateVariable(m.s0, "Temperature", 1600.);
    setExternalStateVariable(m.s1, "Temperature", 1600.);
  };

  for (int i = 0; i < nMat; i++) {
    auto &mat = problem.getMaterial(i + 1);
    set_properties(mat, young1, young2, young3,     // young modulus
                   poisson12, poisson23, poisson13, // poisson ration
                   shear12, shear23, shear13        // shear modulus
    );
    set_temperature(mat);
  }

  //
  const auto vectors = mm_opera_hpc::readVectorsFromFile(p.vect_file);
  if (vectors.size() != 2 * nMat) {
    throw std::invalid_argument(
        "setup_properties : incorrect number of vectors in vector file");
  }

  std::array<mfem_mgis::MaterialAxis3D, 2u> r;
  for (int i = 0; i < nMat; i++) {
    auto &mat = problem.getMaterial(i + 1);
    if (mat.b.symmetry == mgis::behaviour::Behaviour::ORTHOTROPIC) {
      r[0] = vectors[2 * i];
      r[1] = vectors[2 * i + 1];
      mat.setRotationMatrix(mfem_mgis::RotationMatrix3D{r});
    }
  }
}

static void setLinearSolver(mfem_mgis::PeriodicNonLinearEvolutionProblem &p,
                            const int verbosity = 0,
                            const mfem_mgis::real Tol = 1e-12) {
  CatchTimeSection("set_linear_solver");
  // pilote
  constexpr int defaultMaxNumOfIt = 5000; // MaximumNumberOfIterations
  auto solverParameters = mfem_mgis::Parameters{};
  solverParameters.insert(mfem_mgis::Parameters{{"VerbosityLevel", verbosity}});
  solverParameters.insert(
      mfem_mgis::Parameters{{"MaximumNumberOfIterations", defaultMaxNumOfIt}});
  // solverParameters.insert(mfem_mgis::Parameters{{"AbsoluteTolerance", Tol}});
  // solverParameters.insert(mfem_mgis::Parameters{{"RelativeTolerance", Tol}});
  solverParameters.insert(mfem_mgis::Parameters{{"Tolerance", Tol}});

  // preconditionner hypreBoomerAMG
  auto options = mfem_mgis::Parameters{{"VerbosityLevel", verbosity}};
  // auto preconditionner = mfem_mgis::Parameters{{"Name","HypreDiagScale"},
  // {"Options",options}};
  auto preconditionner =
      mfem_mgis::Parameters{{"Name", "HypreBoomerAMG"}, {"Options", options}};
  solverParameters.insert(
      mfem_mgis::Parameters{{"Preconditioner", preconditionner}});
  // solver HyprePCG
  p.setLinearSolver("HyprePCG", solverParameters);
  // p.setLinearSolver("MUMPSSolver", solverParameters);
  // p.setLinearSolver("CGSolver", solverParameters);
}

int main(int argc, char *argv[]) {
  // mpi initialization here
  mfem_mgis::initialize(argc, argv);

  // init timers
  mfem_mgis::Profiler::timers::init_timers();

  // get parameters
  TestParameters p;
  mfem::OptionsParser args(argc, argv);
  common_parameters(args, p);

  // add post processing
  const bool use_post_processing = (p.post_processing == 1);

  // 3D
  constexpr const auto dim = mfem_mgis::size_type{3};

  // creating the finite element workspace
  auto fed = std::make_shared<mfem_mgis::FiniteElementDiscretization>(
      mfem_mgis::Parameters{
          {"MeshFileName", p.mesh_file},
          {"FiniteElementFamily", "H1"},
          {"FiniteElementOrder", p.order},
          {"UnknownsSize", dim},
          {"NumberOfUniformRefinements", p.parallel ? p.refinement : 0},
          {"Parallel", p.parallel}});
  mfem_mgis::PeriodicNonLinearEvolutionProblem problem(fed);
  print_mesh_information(problem);

  // set problem
  mm_opera_hpc::MacroscropicElasticMaterialProperties mp;
  setup_properties(problem, mp, p);
  setLinearSolver(problem, p.verbosity_level);

  // add post processings
  if (use_post_processing)
    add_post_processings(problem, p, "OutputFile-Uniaxial-polycristal");

  // main function here
  const auto te = p.duration;
  const auto nsteps = p.nstep;
  const auto temporal_sequences = [&te, &nsteps] {
    auto times = std::vector<mfem_mgis::real>{};
    const double dt = te / nsteps;
    times.reserve(nsteps + 1);
    for (std::size_t i = 0; i != nsteps + 1; ++i) {
      times.push_back(i * dt);
    }
    return times;
  }();
  //
  const auto Fzz = [](const mfem_mgis::real ets) {
    constexpr auto def = mfem_mgis::real{5e-4};
    return 1 + def * ets;
  };
  //
  std::ofstream out(p.output_file);
  if (!out) {
    std::cerr << "can't open output file '" << p.output_file << "'\n";
  }
  out.precision(14);
  //
  const auto np = mm_opera_hpc::UniaxialMacroscopicStressPeriodicSimulation::
      NumericalParameters{.macroscopic_elastic_material_properties = mp};
  mm_opera_hpc::UniaxialMacroscopicStressPeriodicSimulation s(
      problem, Fzz, np, use_post_processing);
  const auto success = s.run(out, temporal_sequences);

  // print and write timetable
  mfem_mgis::Profiler::timers::print_and_write_timers();
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
