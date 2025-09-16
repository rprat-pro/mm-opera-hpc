#include <cstdlib>
#include <optional>
#include <functional>

#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/ParaviewExportIntegrationPointResultsAtNodes.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MM_OPERA_HPC/GrainOrientations.hxx"
#include "MM_OPERA_HPC/UniaxialMacroscopicStressPeriodicSimulation.hxx"

/*
Problem:

This simulation consists of applying a tensile load on a cermet RVE. The cermet
is a polycrystal where each grain has a material ID (from 2 to Nmat – 1) and a
different orientation. A metallic interface is present between the different
polycrystals (material ID 1)

Definition:

- contact law

MonoCristal_UO2 for all materials

- Parameters:

material: [grain, metal]

- FEM

element: H1
order: 1

*/

struct TestParameters {
  const char *mesh_file = "mesh/5grains.msh";
  const char *behaviourGrain = "MonoCristal_UO2";
  //  const char *behaviourGrain = "Mono_UO2_Cosh_Jaco2";
  const char *behaviourMetal = "NortonCr";
  const char *libraryGrain = "src/libBehaviour.so";
  const char *libraryMetal = "src/libBehaviour.so";
  const char *vector_file = "vectors_5grains.txt";
  int order = 1;
  int refinement = 0;
  int post_processing = 1; // default value : activated
  int verbosity_level = 0; // default value : lower level
  int nsteps = 500;
  double duration = 200.;
};

static void fill_parameters(mfem::OptionsParser &args, TestParameters &p) {
  args.AddOption(&p.mesh_file, "-m", "--mesh", "Mesh file to use.");
  args.AddOption(&p.order, "-o", "--order",
                 "Finite element order (polynomial degree).");
  args.AddOption(&p.refinement, "-r", "--refinement",
                 "refinement level of the mesh, default = 1");
  args.AddOption(&p.post_processing, "-p", "--post-processing",
                 "run post processing step");
  args.AddOption(&p.verbosity_level, "-v", "--verbosity-level",
                 "choose the verbosity level");
  args.AddOption(&p.duration, "-d", "--duration",
                 "choose the duration (default = 5)");
  args.AddOption(&p.nsteps, "-n", "--nsteps",
                 "choose the number of steps (default = 40)");
  args.AddOption(&p.vector_file, "-f", "--file", "Vector file to use");
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


static void print_mesh_information(mfem_mgis::Mesh<true> &mesh) {
  using mfem_mgis::Profiler::Utils::Message;
  using mfem_mgis::Profiler::Utils::sum;

  // get the number of vertices
  int64_t numbers_of_vertices_local = mesh.GetNV();
  int64_t numbers_of_vertices = sum(numbers_of_vertices_local);

  // get the number of elements
  int64_t numbers_of_elements_local = mesh.GetNE();
  int64_t numbers_of_elements = sum(numbers_of_elements_local);

  Message("INFO: number of vertices -> ", numbers_of_vertices);
  Message("INFO: number of elements -> ", numbers_of_elements);
}

static void
setup_material_properties(mfem_mgis::PeriodicNonLinearEvolutionProblem &problem,
                          TestParameters &p,
                          mm_opera_hpc::MacroscropicElasticMaterialProperties &mp) {

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

  mp.update(young1, poisson12, shear12); // used in solve_null_strain

  // ceramic
  auto set_temperature = [](auto &m) {
    setExternalStateVariable(m.s0, "Temperature", 1600);
    setExternalStateVariable(m.s1, "Temperature", 1600);
  };

  auto set_properties_mono =
      [](auto &m, const double yo1, const double yo2, const double yo3,
         const double po12, const double po23, const double po13,
         const double sm12, const double sm23, const double sm13) {
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

  const int nMat =
      getMaterialsAttributes(*(problem.getFiniteElementDiscretizationPointer()))
          .Max();

  const auto vectors = mm_opera_hpc::readVectorsFromFile(p.vector_file);
  if (vectors.size() != 2 * (nMat - 1)) {
    std::cout << vectors.size() << "!=" << 2 * (nMat - 1) << std::endl;
    throw std::invalid_argument(
        "setup_properties : incorrect number of vectors in vector file.");
  }

  // TEST metal
  {
    problem.addBehaviourIntegrator("Mechanics", 1, p.libraryMetal,
                                   p.behaviourMetal);
    auto &metal = problem.getMaterial(1);
    set_temperature(metal);
  }

  for (int grainID = 2; grainID <= nMat; grainID++) {
    problem.addBehaviourIntegrator("Mechanics", grainID, p.libraryGrain,
                                   p.behaviourGrain);
    auto &grain = problem.getMaterial(grainID);
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

int main(int argc, char **argv) {
  using namespace mfem_mgis::Profiler::Utils; // Use Message
                                              // options treatment
  mfem_mgis::initialize(argc, argv);
  mfem_mgis::Profiler::timers::init_timers();

  // get parameters
  TestParameters p;
  mfem::OptionsParser args(argc, argv);
  fill_parameters(args, p);
  // definition of the nonlinear problem
  auto fed = std::make_shared<mfem_mgis::FiniteElementDiscretization>(
      mfem_mgis::Parameters{
          {"MeshFileName", p.mesh_file},
          {"FiniteElementFamily", "H1"},
          {"FiniteElementOrder", p.order},
          {"UnknownsSize", 3},
          {"NumberOfUniformRefinements", p.refinement}, // faster for testing
          {"MeshReadMode", "FromScratch"},
          {"Parallel", true}});
  mfem_mgis::PeriodicNonLinearEvolutionProblem problem(fed);

  // get problem information
  print_mesh_information(fed->getMesh<true>());
  mm_opera_hpc::print_memory_footprint("[Building problem]");

  // choix du solver linéaire +
  int verbosity = p.verbosity_level;
  int post_processing = p.post_processing;
  double Tol = 1e-12;
  int defaultMaxNumOfIt = 2000;
  auto solverParameters = mfem_mgis::Parameters{};
  solverParameters.insert(mfem_mgis::Parameters{{"VerbosityLevel", verbosity}});
  solverParameters.insert(
      mfem_mgis::Parameters{{"MaximumNumberOfIterations", defaultMaxNumOfIt}});
  solverParameters.insert(mfem_mgis::Parameters{{"Tolerance", Tol}});

  //
  auto preconditionner =
      //    mfem_mgis::Parameters{{"Name", "HypreDiagScale"}};
      mfem_mgis::Parameters{
          {"Name", "HypreBoomerAMG"},
          {"Options", mfem_mgis::Parameters{{"VerbosityLevel", verbosity}}}};
  solverParameters.insert(
      mfem_mgis::Parameters{{"Preconditioner", preconditionner}});
  problem.setLinearSolver("HyprePCG", solverParameters);
  // problem.setLinearSolver("HypreGMRES", solverParameters);
  // problem.setLinearSolver("MUMPSSolver", {});

  //
  const auto Fzz = [](const mfem_mgis::real ets) {
    constexpr auto def = mfem_mgis::real{5e-4};
    return 1 + def * ets;
  };
  //
  mm_opera_hpc::MacroscropicElasticMaterialProperties mp;
  setup_material_properties(problem, p, mp);

  if (post_processing) {
    mfem_mgis::Profiler::Utils::Message("Define post processings");
    problem.addPostProcessing("ParaviewExportResults",
                              {{"OutputFileName", "Displacement"}});
    problem.addPostProcessing("MeanThermodynamicForces",
                              {{"OutputFileName", "avgStress"}});
    /* BUG
       problem.addPostProcessing(
       "ParaviewExportIntegrationPointResultsAtNodes",
       {{"OutputFileName", "IntegrationPointOutput"}});
     */
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
  const auto np = mm_opera_hpc::UniaxialMacroscopicStressPeriodicSimulation::
      NumericalParameters{.macroscopic_elastic_material_properties = mp};
  mm_opera_hpc::UniaxialMacroscopicStressPeriodicSimulation s(problem, Fzz, np,
                                                              post_processing);
  const auto success = s.run(temporal_sequences);
  //
  mm_opera_hpc::print_memory_footprint("[End]");
  mfem_mgis::Profiler::timers::print_and_write_timers();
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
