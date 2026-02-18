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

#include "mfem/fem/datacollection.hpp"
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
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/L2Projection.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include "MM_OPERA_HPC/Utilities.hxx"
#include "MM_OPERA_HPC/GrainOrientations.hxx"
#include "MM_OPERA_HPC/MacroscropicElasticMaterialProperties.hxx"
#include "MM_OPERA_HPC/UniaxialMacroscopicStressPeriodicSimulation.hxx"

/*!
 * \brief parameters related to the mesh
 */
struct MeshParameters {
  //! \brief default value for the mesh file
  const char *mesh_file = "mesh/5crystals.msh";
  //! \brief number of uniform refinement of the mesh
  int refinement = 0;
};

/*!
 * \brief parameters related to the numerical aspects of the simulation
 */
struct NumericalParameters {
  //! \brief polynomial order of the finite element used
  int order = 1;
  //! \brief name of the linear solver
  const char *linear_solver = "HyprePCG";
  //! \brief name of the preconditioner for the linear solver
  const char *linear_solver_preconditioner = "HypreBoomerAMG";
  //! \brief default tolerance for the convergence of the linear solver
  mfem_mgis::real linear_solver_tolerance = mfem_mgis::real{1e-12};
};

/*!
 * \brief parameters describing the elastic material properties
 */
struct ElasticMaterialPropertiesParameters {
  //! \brief Young's modulus of uranium dioxide
  mfem_mgis::real young_modulus = mfem_mgis::real{222.e9};
  //! \brief Poisson's ratio of uranium dioxide
  mfem_mgis::real poisson_ratio = mfem_mgis::real{0.27};
  //! \brief shear modulus of uranium dioxide
  mfem_mgis::real shear_modulus = mfem_mgis::real{54.e9};
};

struct MaterialParameters : ElasticMaterialPropertiesParameters {
  //! \brief default material library
  const char *library = "src/libBehaviour.so";
  //! \brief default mechanical behaviour
  const char *behaviour = "Mono_UO2_Cosh_Jaco3";
  //! \brief default file name describing orientation vectors
  const char *vect_file = "mesh/vectors_5crystals.txt";
};

struct LoadingParameters {
  //! \brief default value of the temperature
  mfem_mgis::real temperature = mfem_mgis::real{1600};
  //! \brief default value of the imposed linear strain rate
  mfem_mgis::real linear_strain_rate = mfem_mgis::real{5e-4};
};

struct TimeDiscretizationParameters {
  //! \brief total duration of the simulated test
  mfem_mgis::real duration = mfem_mgis::real{200};
  //! \brief number of step describing the temporal sequences
  int nstep = 200;
};

struct PostProcessingParameters {
  //! \brief boolean state if post-processings shall be executed
  bool post_processings = true;
  //! \brief boolean stating if the von Mises stress of the Cauchy stress is
  //! exported to VTK
  bool export_von_Mises_stress = false;
  //! \brief boolean stating if first eigen stress of the Cauchy stress is
  //! exported to VTK
  bool export_first_eigen_stress = false;
};

// main parameters
struct TestParameters : MeshParameters,
    NumericalParameters,
    MaterialParameters,
    LoadingParameters,
    TimeDiscretizationParameters,
    PostProcessingParameters {
      //! \brief default output file for macroscopic results
      const char *output_file = "uniaxial-polycrystal.res";
      //! \brief default verbositiy level
      int verbosity_level = 0;  // default value : lower level
    };

// a few utility functions defined after the main functions

static void parseCommandLineArguments(mfem::OptionsParser &, TestParameters &);
static void addPostProcessings(mfem_mgis::Context&,
                               mfem_mgis::PeriodicNonLinearEvolutionProblem &,
                               const std::string &);
static void addOptionalPostProcessings(mfem_mgis::Context&,
                                       mfem_mgis::PeriodicNonLinearEvolutionProblem &,
                                       const PostProcessingParameters &,
                                       mfem::ParaViewDataCollection*&);
static void setupMaterials(
    mfem_mgis::PeriodicNonLinearEvolutionProblem &,
    mm_opera_hpc::MacroscropicElasticMaterialProperties &,
    const TestParameters &);
static void setLinearSolver(mfem_mgis::PeriodicNonLinearEvolutionProblem &,
                            const TestParameters &);

int main(int argc, char *argv[]) {
  auto ctx = mfem_mgis::Context{};
  // mpi initialization here
  mfem_mgis::initialize(argc, argv);

  // init timers
  mfem_mgis::Profiler::timers::init_timers();

  // get parameters
  TestParameters p;
  mfem::OptionsParser args(argc, argv);
  parseCommandLineArguments(args, p);

  // creating the finite element workspace
  auto fed = std::make_shared<mfem_mgis::FiniteElementDiscretization>(
      mfem_mgis::Parameters{{"MeshFileName", p.mesh_file},
      {"FiniteElementFamily", "H1"},
      {"FiniteElementOrder", p.order},
      {"UnknownsSize", mfem_mgis::size_type{3}},
      {"NumberOfUniformRefinements", p.refinement},
      {"Parallel", true}});
  mfem_mgis::PeriodicNonLinearEvolutionProblem problem(fed);
  mm_opera_hpc::printMeshInformation(problem);

  // set problem
  mm_opera_hpc::MacroscropicElasticMaterialProperties mp;
  setupMaterials(problem, mp, p);
  setLinearSolver(problem, p);
  // add post processings
  mfem::ParaViewDataCollection* paraview_exporter;
  if (p.post_processings) {
    addPostProcessings(ctx, problem, "OutputFile-Uniaxial-polycrystal");
    addOptionalPostProcessings(ctx, problem, p, paraview_exporter);
  }
  // definition of the temporal sequences
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
  // evolution of the axial component of the deformation gradient
  const auto Fzz = [&p](const mfem_mgis::real ets) {
    return 1 + p.linear_strain_rate * ets;
  };
  // main output file
  std::ofstream out(p.output_file);
  if (!out) {
    std::cerr << "can't open output file '" << p.output_file << "'\n";
  }
  out.precision(14);
  // setting the simulation and running it
  const auto np = mm_opera_hpc::UniaxialMacroscopicStressPeriodicSimulation::
      NumericalParameters{.macroscopic_elastic_material_properties = mp};
  mm_opera_hpc::UniaxialMacroscopicStressPeriodicSimulation s(
      problem, Fzz, np, p.post_processings);
  const auto success = s.run(out, temporal_sequences);

  delete paraview_exporter;
  // print and write timetable
  mfem_mgis::Profiler::timers::print_and_write_timers();
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}

static void parseCommandLineArguments(mfem::OptionsParser &args,
                                      TestParameters &p) {
  args.AddOption(&p.output_file, "", "--macroscopic-stress-output-file",
                 "main output file containing the evolution of the diagonal "
                 "components of the deformation gradient and the  diagonal "
                 "components of the Cauchy stress.");
  args.AddOption(&p.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&p.vect_file, "-f", "--vect", "Vector file to use.");
  args.AddOption(&p.library, "-l", "--library", "Material library.");
  args.AddOption(&p.behaviour, "-b", "--behaviour", "Mechanical behaviour.");
  args.AddOption(&p.linear_solver, "", "--linear-solver", "linear solver to be used");
  args.AddOption(&p.linear_solver_preconditioner, "",
                 "--linear-solver-preconditioner", "preconditioner of the linear solver to be used");
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
  args.AddOption(&p.post_processings, "", "--enable-post-processing", "",
                 "--disable-post-processing", "run post processing steps",
                 false);
  args.AddOption(&p.verbosity_level, "-v", "--verbosity-level",
                 "choose the verbosity level");
  args.AddOption(&p.duration, "-d", "--duration",
                 "choose the duration (default = 5)");
  args.AddOption(&p.nstep, "-n", "--nstep",
                 "choose the number of steps (default = 40)");

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

static void setupMaterials(
    mfem_mgis::PeriodicNonLinearEvolutionProblem &problem,
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
  const mfem_mgis::real young1 = p.young_modulus;
  const mfem_mgis::real young2 = young1;
  const mfem_mgis::real young3 = young1;
  const mfem_mgis::real poisson12 = p.poisson_ratio;
  const mfem_mgis::real poisson23 = poisson12;
  const mfem_mgis::real poisson13 = poisson12;
  const mfem_mgis::real shear12 = p.shear_modulus;
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
        for (auto s : {&m.s0, &m.s1}) {
          setMaterialProperty(*s, "YoungModulus1", yo1);
          setMaterialProperty(*s, "YoungModulus2", yo2);
          setMaterialProperty(*s, "YoungModulus3", yo3);
          setMaterialProperty(*s, "PoissonRatio12", po12);
          setMaterialProperty(*s, "PoissonRatio23", po23);
          setMaterialProperty(*s, "PoissonRatio13", po13);
          setMaterialProperty(*s, "ShearModulus12", sm12);
          setMaterialProperty(*s, "ShearModulus23", sm23);
          setMaterialProperty(*s, "ShearModulus13", sm13);
        }
      };

  auto set_temperature = [&p](auto &m) {
    setExternalStateVariable(m.s0, "Temperature", p.temperature);
    setExternalStateVariable(m.s1, "Temperature", p.temperature);
  };

  for (int i = 0; i < nMat; i++) {
    auto &mat = problem.getMaterial(i + 1);
    set_properties(mat, young1, young2, young3,      // young modulus
                   poisson12, poisson23, poisson13,  // poisson ration
                   shear12, shear23, shear13         // shear modulus
                  );
    set_temperature(mat);
  }

  //
  const auto vectors = mm_opera_hpc::readVectorsFromFile(p.vect_file);
  if (vectors.size() != 2 * nMat) {
    throw std::invalid_argument(
        "setupMaterials : incorrect number of vectors in vector file");
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
                            const TestParameters &params) {
  CatchTimeSection("set_linear_solver");
  // pilote
  constexpr int defaultMaxNumOfIt = 5000;  // MaximumNumberOfIterations
  auto solverParameters = mfem_mgis::Parameters{};
  solverParameters.insert(
      mfem_mgis::Parameters{{"VerbosityLevel", params.verbosity_level}});
  solverParameters.insert(
      mfem_mgis::Parameters{{"MaximumNumberOfIterations", defaultMaxNumOfIt}});
  // solverParameters.insert(mfem_mgis::Parameters{{"AbsoluteTolerance", Tol}});
  // solverParameters.insert(mfem_mgis::Parameters{{"RelativeTolerance", Tol}});
  solverParameters.insert(
      mfem_mgis::Parameters{{"Tolerance", params.linear_solver_tolerance}});

  // preconditioner
  auto options =
      mfem_mgis::Parameters{{"VerbosityLevel", params.verbosity_level}};
  // auto preconditioner = mfem_mgis::Parameters{{"Name","HypreDiagScale"},
  // {"Options",options}};
  if(std::string_view{params.linear_solver_preconditioner} != "none"){
    auto preconditioner = mfem_mgis::Parameters{
      {"Name", params.linear_solver_preconditioner}, {"Options", options}};
    solverParameters.insert(mfem_mgis::Parameters{{"Preconditioner", preconditioner}});
  }
  p.setLinearSolver(params.linear_solver, solverParameters);
}


static void addOptionalPostProcessings(mfem_mgis::Context& ctx,
                                       mfem_mgis::PeriodicNonLinearEvolutionProblem &p,
                                       const PostProcessingParameters &params,
                                       mfem::ParaViewDataCollection *&exporter) {
  // setup paraview
  exporter = new mfem::ParaViewDataCollection("PolycrystalParaviewOutput");
  exporter->SetDataFormat(mfem::VTKFormat::BINARY);
  exporter->SetMesh(&(p.getFiniteElementDiscretization().getMesh<true>()));

  if (params.export_von_Mises_stress || params.export_first_eigen_stress) {
    p.getImplementation<true>().addPostProcessing(
        [&ctx, &p, &params, &exporter] (mfem_mgis::real t, mfem_mgis::real dt) {
        static int count = 0;
        exporter->SetCycle(count++);
        exporter->SetTime(t);
        const int nMat =
        getMaterialsAttributes(*(p.getFiniteElementDiscretizationPointer()))
        .Max();
        // fes = FirstEigenStress
        // vmes = VonMisesEquivalentStress
        std::vector<mfem_mgis::ImmutablePartialQuadratureFunctionView> fesView{};
        std::vector<mfem_mgis::ImmutablePartialQuadratureFunctionView> vmesView{};
        auto ets = mfem_mgis::Material::END_OF_TIME_STEP;
        // For each material
        for (int i = 1; i <= nMat; i++) {
        auto& m = p.getMaterial(i);
        if (params.export_first_eigen_stress) {
        auto fes_pqf = mfem_mgis::PartialQuadratureFunction{
        m.getPartialQuadratureSpacePointer(), 1};
        mfem_mgis::computeFirstEigenStress(ctx, fes_pqf, m, ets);
        fesView.push_back(fes_pqf);
        }
        if (params.export_von_Mises_stress) {
          auto vmes_pqf = mfem_mgis::PartialQuadratureFunction{
            m.getPartialQuadratureSpacePointer(), 1};
          bool success = mfem_mgis::computeVonMisesEquivalentStress(ctx, vmes_pqf, m, ets);
          if (!success) {
            std::cout << "ERROR: computeVonMisesEquivalentStress returns"
                << " false in addOptionalPostProcessings" << std::endl;
          }
          vmesView.push_back(vmes_pqf);
        }
        }

        // Used by computeL2Projection
        auto& lsf = mfem_mgis::LinearSolverFactory<true>::getFactory();
        auto& fespace =
            p.getFiniteElementDiscretization().getFiniteElementSpace<true>();
        auto linear_solver = lsf.generate(
            ctx, "HyprePCG", fespace,
            mfem_mgis::Parameters{
            {"Preconditioner", mfem_mgis::Parameters{{"Name", "HypreILU"}}}});

        // project solutions
        std::optional<mfem_mgis::L2ProjectionResult<true>> fes_projection;
        std::optional<mfem_mgis::L2ProjectionResult<true>> vmes_projection;
        if (params.export_first_eigen_stress) {
          fes_projection =
              mfem_mgis::computeL2Projection<true>(ctx, linear_solver, fesView);
          exporter->RegisterField("FirstEigenStress",
                                  fes_projection->result.get());
        }
        if (params.export_von_Mises_stress) {
          vmes_projection =
              mfem_mgis::computeL2Projection<true>(ctx, linear_solver, vmesView);
          exporter->RegisterField("vonMisesStress",
                                  vmes_projection->result.get());
        }

        // save paraview files
        exporter->Save();
        });
  }
}

static void addPostProcessings(mfem_mgis::Context& ctx,
                               mfem_mgis::PeriodicNonLinearEvolutionProblem &p,
                               const std::string &msg) {
  p.addPostProcessing("ParaviewExportResults", {{"OutputFileName", msg}});
  p.addPostProcessing("MeanThermodynamicForces",
                      {{"OutputFileName", "avgStressPolycrystal"}});
}  // end timer add_postprocessing_and_outputs
