#include <fstream>
#include <cstdlib>
#include <optional>
#include <functional>

#include "mfem/fem/datacollection.hpp"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/ParaviewExportIntegrationPointResultsAtNodes.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/MechanicalPostProcessings.hxx"
#include "MFEMMGIS/LinearSolverFactory.hxx"
#include "MFEMMGIS/L2Projection.hxx"
#include "MM_OPERA_HPC/Utilities.hxx"
#include "MM_OPERA_HPC/GrainOrientations.hxx"
#include "MM_OPERA_HPC/UniaxialMacroscopicStressPeriodicSimulation.hxx"

/*
Problem:

This simulation consists of applying a tensile load on a cermet RVE. The cermet
is a polycrystal where each grain has a material ID (from 2 to Nmat) and a
different orientation. A metallic interface is present between the different
polycrystals (material ID 1)

Definition:

- [Crystal] Behavior law: Mono_UO2_Cosh_Jaco3
-- Parameters:
--- Young Modulus = 222.e9
--- Poisson ratio = 0.27
--- Shear Modulus = 54.e9
- [Interface] Behavior law: Norton
-- Parameters (@see law/NortonCr.mfront):
--- Young Modulus = 276e+09
--- Poisson ratio = 0.3
--- A             = 2.5e+11
--- n1            = 4.75
--- Q             = 306.27e+03
--- D0            = 1.55e-5
--- b             = 2.5e-10
- FEM Parameters
-- Element: H1
-- Order: 1
- T_final: 200s
- Number of timesteps: 500

*/

struct TestParameters {
  const char *output_file = "cermet.res";
  const char *mesh_file = "mesh/5grains.msh";
  //  const char *behaviourGrain = "MonoCristal_UO2";
  const char *behaviourGrain = "Mono_UO2_Cosh_Jaco3";
  const char *behaviourMetal = "NortonCr";
  const char *libraryGrain = "src/libBehaviour.so";
  const char *libraryMetal = "src/libBehaviour.so";
  const char *vector_file = "vectors_5grains.txt";
  int order = 1;
  int refinement = 0;
  int post_processing = 1;  // default value : activated
  //! \brief boolean stating if the von Mises stress of the Cauchy stress is
  //! exported to VTK
  bool export_von_Mises_stress = false;
  //! \brief boolean stating if first eigen stress of the Cauchy stress is
  //! exported to VTK
  bool export_first_eigen_stress = false;
  int verbosity_level = 0;  // default value : lower level
  int nsteps = 500;
  double duration = 200.;
};

static void parseCommandLineOptions(mfem::OptionsParser &args,
                                    TestParameters &p) {
  args.AddOption(&p.output_file, "", "--macroscopic-stress-output-file",
                 "main output file containing the evolution of the diagonal "
                 "components of the deformation gradient and the  diagonal "
                 "components of the Cauchy stress.");
  args.AddOption(&p.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&p.order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&p.refinement, "-r", "--refinement",
                 "refinement level of the mesh");
  args.AddOption(&p.post_processing, "-p", "--postprocessing",
                 "run post processing step");
  args.AddOption(&p.export_von_Mises_stress, "",
                 "--enable-export-von-Mises-stress", "",
                 "--disable-export-von-Mises-stress",
                 "export the von Mises stress", false);
  args.AddOption(&p.export_first_eigen_stress, "",
                 "--enable-export-first_eigen_stress", "",
                 "--disable-export-first_eigen_stress",
                 "export first eigen stress", false);
  args.AddOption(&p.verbosity_level, "-v", "--verbosity-level",
                 "choose the verbosity level");
  args.AddOption(&p.duration, "-d", "--duration",
                 "choose the duration");
  args.AddOption(&p.nsteps, "-n", "--nstep",
                 "choose the number of steps");
  args.AddOption(&p.vector_file, "-f", "--file", "Vector file to use");
  args.Parse();
  if (!args.Good()) {
    if (mfem_mgis::isMPIMainProcessForDefaultCommunicator()) args.PrintUsage(std::cout);
    mfem_mgis::finalize();
    exit(args.Help() ? EXIT_SUCCESS : EXIT_FAILURE);
  }
  if (p.mesh_file == nullptr) {
    if (mfem_mgis::isMPIMainProcessForDefaultCommunicator())
      std::cout << "ERROR: Mesh file missing" << std::endl;
    args.PrintUsage(std::cout);
    mfem_mgis::abort(EXIT_FAILURE);
  }
  if (mfem_mgis::isMPIMainProcessForDefaultCommunicator()) args.PrintOptions(std::cout);
}

static void setup_material_properties(
    mfem_mgis::attributes::MayAbort,
    mfem_mgis::Context &ctx,
    mfem_mgis::PeriodicNonLinearEvolutionProblem &problem,
    TestParameters &p,
    mm_opera_hpc::MacroscopicElasticMaterialProperties &mp) {
  auto or_die = ctx.getFatalFailureHandler();
  // cubic symmetry elasticity
  const double young1 = 222.e9;
  const double young2 = young1;
  const double young3 = young1;
  const double poisson12 = 0.27;
  const double poisson23 = poisson12;
  const double poisson13 = poisson12;
  const double shear12 = 54.e9;
  const double shear23 = shear12;
  const double shear13 = shear12;

  mp.update(young1, poisson12, shear12);  // used by the fixed-point algorithm

  // ceramic
  auto set_temperature = [&ctx, &or_die](auto &m) {
    setExternalStateVariable(ctx, m.s0, "Temperature", 1600) | or_die;
    setExternalStateVariable(ctx, m.s1, "Temperature", 1600) | or_die;
  };

  auto set_properties_mono =
      [&ctx, &or_die](auto &m, const double yo1, const double yo2,
                      const double yo3, const double po12, const double po23,
                      const double po13, const double sm12, const double sm23,
                      const double sm13) {
        setMaterialProperty(ctx, m.s0, "YoungModulus1", yo1) | or_die;
        setMaterialProperty(ctx, m.s0, "YoungModulus2", yo2) | or_die;
        setMaterialProperty(ctx, m.s0, "YoungModulus3", yo3) | or_die;
        setMaterialProperty(ctx, m.s0, "PoissonRatio12", po12) | or_die;
        setMaterialProperty(ctx, m.s0, "PoissonRatio23", po23) | or_die;
        setMaterialProperty(ctx, m.s0, "PoissonRatio13", po13) | or_die;
        setMaterialProperty(ctx, m.s0, "ShearModulus12", sm12) | or_die;
        setMaterialProperty(ctx, m.s0, "ShearModulus23", sm23) | or_die;
        setMaterialProperty(ctx, m.s0, "ShearModulus13", sm13) | or_die;

        setMaterialProperty(ctx, m.s1, "YoungModulus1", yo1) | or_die;
        setMaterialProperty(ctx, m.s1, "YoungModulus2", yo2) | or_die;
        setMaterialProperty(ctx, m.s1, "YoungModulus3", yo3) | or_die;
        setMaterialProperty(ctx, m.s1, "PoissonRatio12", po12) | or_die;
        setMaterialProperty(ctx, m.s1, "PoissonRatio23", po23) | or_die;
        setMaterialProperty(ctx, m.s1, "PoissonRatio13", po13) | or_die;
        setMaterialProperty(ctx, m.s1, "ShearModulus12", sm12) | or_die;
        setMaterialProperty(ctx, m.s1, "ShearModulus23", sm23) | or_die;
        setMaterialProperty(ctx, m.s1, "ShearModulus13", sm13) | or_die;
      };

  const int nMat =
      getMaterialsAttributes(problem.getFiniteElementDiscretization())
          .Max();

  const auto vectors = mm_opera_hpc::readVectorsFromFile(p.vector_file);
  if (vectors.size() != 2 * (nMat - 1)) {
    std::cout << vectors.size() << "!=" << 2 * (nMat - 1) << std::endl;
    throw std::invalid_argument(
        "setup_properties : incorrect number of vectors in vector file.");
  }

  // TEST metal
  {
    problem.addBehaviourIntegrator(ctx, "Mechanics", 1, p.libraryMetal,
                                   p.behaviourMetal) |
        or_die;
    auto &metal = problem.getMaterial(ctx, 1, 0) | or_die;
    set_temperature(metal);
  }

  for (int grainID = 2; grainID <= nMat; grainID++) {
    problem.addBehaviourIntegrator(ctx, "Mechanics", grainID, p.libraryGrain,
                                   p.behaviourGrain) |
        or_die;
    auto &grain = problem.getMaterial(ctx, grainID, 0) | or_die;
    set_properties_mono(grain, young1, young2, young3, poisson12, poisson23,
                        poisson13, shear12, shear23, shear13);
    set_temperature(grain);

    std::array<mfem_mgis::MaterialAxis3D, 2u> r;
    if (grain.b.symmetry == mgis::behaviour::Behaviour::ORTHOTROPIC) {
      r[0] = vectors.at(2 * (grainID - 2));
      r[1] = vectors.at(2 * (grainID - 2) + 1);
      grain.setRotationMatrix(mfem_mgis::RotationMatrix3D{r});
    }
  }
}

static void addOptionalPostProcessings(
    mfem_mgis::attributes::MayAbort,
    mfem_mgis::Context &ctx,
    mfem_mgis::PeriodicNonLinearEvolutionProblem &p,
    const TestParameters &params,
    mfem::ParaViewDataCollection *&exporter) {
  auto or_die = ctx.getFatalFailureHandler();
  // setup paraview
  exporter = new mfem::ParaViewDataCollection("CermetParaviewOutput");
  exporter->SetDataFormat(mfem::VTKFormat::BINARY);
  exporter->SetMesh(&(p.getFiniteElementDiscretization().getMesh<true>()));

  if (params.export_von_Mises_stress || params.export_first_eigen_stress) {
    p.getImplementation<true>().addPostProcessing(ctx, [&ctx, &p, &params,
                                                        &exporter, or_die](
                                                           mfem_mgis::real t,
                                                           mfem_mgis::real dt) {
      static int count = 0;
      exporter->SetCycle(count++);
      // the stresses are computed at the end of the time step
      exporter->SetTime(t + dt);
      const int nMat =
          getMaterialsAttributes(p.getFiniteElementDiscretization())
              .Max();
      // fes = FirstEigenStress
      // vmes = VonMisesEquivalentStress
      std::vector<mfem_mgis::ImmutablePartialQuadratureFunctionView> fesView{};
      std::vector<mfem_mgis::ImmutablePartialQuadratureFunctionView> vmesView{};
      // the functions must outlive the views used by the projections
      std::vector<mfem_mgis::PartialQuadratureFunction> fes_pqfs;
      std::vector<mfem_mgis::PartialQuadratureFunction> vmes_pqfs;
      fes_pqfs.reserve(nMat);
      vmes_pqfs.reserve(nMat);
      auto ets = mfem_mgis::Material::END_OF_TIME_STEP;
      // For each material: the metallic interface and the grains
      for (int i = 1; i <= nMat; i++) {
        auto &m = p.getMaterial(ctx, i, 0) | or_die;
        if (params.export_first_eigen_stress) {
          auto &fes_pqf =
              fes_pqfs.emplace_back(m.getPartialQuadratureSpace(), 1);
          mfem_mgis::computeFirstEigenStress(ctx, fes_pqf, m, ets);
          fesView.push_back(fes_pqf);
        }
        if (params.export_von_Mises_stress) {
          auto &vmes_pqf =
              vmes_pqfs.emplace_back(m.getPartialQuadratureSpace(), 1);
          bool success =
              mfem_mgis::computeVonMisesEquivalentStress(ctx, vmes_pqf, m, ets);
          if (!success) {
            std::cout << "ERROR: computeVonMisesEquivalentStress returns"
                      << " false in addOptionalPostProcessings" << std::endl;
          }
          vmesView.push_back(vmes_pqf);
        }
      }

      // Used by computeL2Projection
      auto &lsf = mfem_mgis::LinearSolverFactory<true>::getFactory();
      auto &fespace =
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
    }) | or_die;
  }
}

int main(int argc, char **argv) {
  using namespace mfem_mgis::Profiler::Utils;  // Use Message
                                               // options treatment
  mfem_mgis::initialize(argc, argv);
  auto ctx = mfem_mgis::Context{};
  auto or_die = ctx.getFatalFailureHandler();
  ctx.enableProfiling(true);

  // get parameters
  TestParameters p;
  mfem::OptionsParser args(argc, argv);
  parseCommandLineOptions(args, p);
  // definition of the nonlinear problem
  auto fed =
      mfem_mgis::construct<mfem_mgis::FiniteElementDiscretization>(
          ctx, mfem_mgis::Parameters{{"MeshFileName", p.mesh_file},
                                     {"FiniteElementFamily", "H1"},
                                     {"FiniteElementOrder", p.order},
                                     {"UnknownsSize", 3},
                                     {"NumberOfUniformRefinements",
                                      p.refinement},  // faster for testing
                                     {"MeshReadMode", "FromScratch"},
                                     {"Parallel", true}}) |
      or_die;
  auto problem =
      mfem_mgis::construct<mfem_mgis::PeriodicNonLinearEvolutionProblem>(ctx,
                                                                         fed) |
      or_die;

  // get problem information
  mm_opera_hpc::printMeshInformation(ctx, problem);
  mm_opera_hpc::printMemoryFootprint("[Building problem]");

  // choix du solver linéaire +
  int verbosity = p.verbosity_level;
  int post_processing = p.post_processing;
  double Tol = 1e-12;
  int defaultMaxNumOfIt = 2000;
  auto solverParameters = mfem_mgis::Parameters{};
  solverParameters.insert(mfem_mgis::may_throw,
                          mfem_mgis::Parameters{{"VerbosityLevel", verbosity}});
  solverParameters.insert(
      mfem_mgis::may_throw,
      mfem_mgis::Parameters{{"MaximumNumberOfIterations", defaultMaxNumOfIt}});
  solverParameters.insert(mfem_mgis::may_throw,
                          mfem_mgis::Parameters{{"Tolerance", Tol}});

  //
  auto preconditionner =
      //    mfem_mgis::Parameters{{"Name", "HypreDiagScale"}};
      mfem_mgis::Parameters{
          {"Name", "HypreBoomerAMG"},
          {"Options", mfem_mgis::Parameters{{"VerbosityLevel", verbosity}}}};
  solverParameters.insert(
      mfem_mgis::may_throw,
      mfem_mgis::Parameters{{"Preconditioner", preconditionner}});
  problem.setLinearSolver(ctx, "HyprePCG", solverParameters) | or_die;
  // problem.setLinearSolver("HypreGMRES", solverParameters);
  // problem.setLinearSolver("MUMPSSolver", {});

  //
  const auto Fzz = [](const mfem_mgis::real ets) {
    constexpr auto def = mfem_mgis::real{5e-4};
    return 1 + def * ets;
  };
  //
  mm_opera_hpc::MacroscopicElasticMaterialProperties mp;
  setup_material_properties(mfem_mgis::may_abort, ctx, problem, p, mp);

  mfem::ParaViewDataCollection *paraview_exporter = nullptr;
  if (post_processing) {
    mfem_mgis::Profiler::Utils::Message("Define post processings");
    problem.addPostProcessing(ctx, "ParaviewExportResults",
                              {{"OutputFileName", "Displacement"}}) |
        or_die;
    problem.addPostProcessing(ctx, "MeanThermodynamicForces",
                              {{"OutputFileName", "avgStress"}}) |
        or_die;
    /* BUG
       problem.addPostProcessing(
       "ParaviewExportIntegrationPointResultsAtNodes",
       {{"OutputFileName", "IntegrationPointOutput"}});
     */
    addOptionalPostProcessings(mfem_mgis::may_abort, ctx, problem, p,
                               paraview_exporter);
  }
  //
  const auto te = p.duration;
  const auto nsteps = p.nsteps;
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
  std::ofstream out(p.output_file);
  if (!out) {
    std::cerr << "can't open output file '" << p.output_file << "'\n";
  }
  out.precision(14);
  //
  const auto np = mm_opera_hpc::UniaxialMacroscopicStressPeriodicSimulation::
      NumericalParameters{.macroscopic_elastic_material_properties = mp};
  mm_opera_hpc::UniaxialMacroscopicStressPeriodicSimulation s(problem, Fzz, np,
                                                              post_processing);
  const auto success = s.run(ctx, out, temporal_sequences);

  if (paraview_exporter != nullptr) {
    delete paraview_exporter;
  }
  //
  mm_opera_hpc::printMemoryFootprint("[End]");
  mfem_mgis::Profiler::OutputManager::printTimeTable(ctx);
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
