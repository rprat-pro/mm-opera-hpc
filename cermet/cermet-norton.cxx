#include <cstdlib>
#include <sys/resource.h>
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/ParaviewExportIntegrationPointResultsAtNodes.hxx"


/*

Parameters : 

start time = 0
end time = 5s
number of time step = 40

Strain Gradient matrix : val = 0.012
[ val / 2 ,         0 ,   0 ]
[ 0         , val / 2 ,   0 ] 
[ 0         ,          0, -val ]

Solver : HyprePCG
Preconditionner : HypreBoomerAMG

Behavior law parameters : SaintVenantKirchhoffElasticity, ImplicitNortonThreshold
[ parameters       , grain    , metal      ];    
[ Young Modulus    , 8.182e9  , 2*8.182e9  ];
[ Poisson Ratio    , 0.364    , 0.364      ];
[ Stress Threshold , 100.0e6  , 100.0e12   ];
[ Norton Exponent  , 3.333333 , 3.333333   ];
[ Temperature      , 293.15   , 293.15     ];
Element :

Familly H1
Order 1
 */

// display information
long get_memory_checkpoint()
{
  rusage obj;
  int who = 0;
  [[maybe_unused]] auto test = getrusage(who, &obj);
  assert((test = -1) && "error: getrusage has failed");
  long res;
  MPI_Reduce(&(obj.ru_maxrss), &(res), 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  return res;
};
void print_memory_footprint(std::string msg)
{
  long mem = get_memory_checkpoint();
  double m = double(mem) * 1e-6; // conversion kb to Gb
  mfem_mgis::Profiler::Utils::Message(msg, " memory footprint: ", m, " GB");
}

template <typename Mesh>
void print_mesh_information(Mesh& mesh) {
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


// add postprocessing
template <typename Problem>
void post_process(Problem &p, double start, double end) {
  CatchTimeSection("common::post_processing_step");
  p.executePostProcessings(start, end);
}
struct TestParameters {
  //  const char *metal = "MO";
  const char *mesh_file = "mesh/cermet-mini.msh";
  //const char *behaviourGrain = "Elasticity";
  const char *behaviourGrain = "ImplicitNortonThreshold";
  const char *behaviourMetal = "ImplicitNortonThreshold";
  const char *libraryGrain = "src/libBehaviour.so";
  const char *libraryMetal = "src/libBehaviour.so";
  int order = 1;
  int refinement = 1;
  int post_processing = 1; // default value : activated
  int verbosity_level = 0; // default value : lower level
  int nstep = 40;
  double duration = 5.;
};

void fill_parameters(mfem::OptionsParser &args, TestParameters &p) {
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
      mfem_mgis::Parameters{{"MeshFileName", p.mesh_file},
      {"FiniteElementFamily", "H1"},
      {"FiniteElementOrder", p.order},
      {"UnknownsSize", 3},
      {"NumberOfUniformRefinements", p.refinement}, // faster for testing
      {"MeshReadMode", "FromScratch"},
      {"Parallel", true}});
  mfem_mgis::PeriodicNonLinearEvolutionProblem problem(fed);

  // get problem information
  print_mesh_information(fed->getMesh<true>());
  print_memory_footprint("[Building problem]");

  // non linear solver
  problem.setSolverParameters({{"VerbosityLevel", 1},
      {"RelativeTolerance", 1e-6},
      {"AbsoluteTolerance", 0.},
      {"MaximumNumberOfIterations", 6}});
  // choix du solver linéaire +
  int verbosity = p.verbosity_level;
  int post_processing = p.post_processing;
  double Tol = 1e-6;
  int defaultMaxNumOfIt = 2000;
  auto solverParameters = mfem_mgis::Parameters{};
  solverParameters.insert(mfem_mgis::Parameters{{"VerbosityLevel", verbosity}});
  solverParameters.insert(
      mfem_mgis::Parameters{{"MaximumNumberOfIterations", defaultMaxNumOfIt}});
  solverParameters.insert(mfem_mgis::Parameters{{"Tolerance", Tol}});

  //
  auto preconditionner =
    mfem_mgis::Parameters{{"Name", "HypreDiagScale"}};
  //    mfem_mgis::Parameters{{"Name", "HypreBoomerAMG"}};
  solverParameters.insert(
      mfem_mgis::Parameters{{"Preconditioner", preconditionner}});
  problem.setLinearSolver("HyprePCG", solverParameters);

  // Boundary Conditions
  // TO DO
  // ceramic
  problem.addBehaviourIntegrator("Mechanics", 1, p.libraryGrain, p.behaviourGrain);
  problem.addBehaviourIntegrator("Mechanics", 2, p.libraryMetal, p.behaviourMetal);
  auto& grain = problem.getMaterial(1);
  auto& metal = problem.getMaterial(2);


  auto set_properties_elasticity = [](auto& m, const double l, const double mu) {
    mgis::behaviour::setMaterialProperty(m.s0, "FirstLameCoefficient", l);
    mgis::behaviour::setMaterialProperty(m.s0, "ShearModulus", mu);
    mgis::behaviour::setMaterialProperty(m.s1, "FirstLameCoefficient", l);
    mgis::behaviour::setMaterialProperty(m.s1, "ShearModulus", mu);
  };

  auto set_properties_norton = [](auto& m, const double yo, const double po, const double st, const double no) 
  {
    setMaterialProperty(m.s0, "YoungModulus", yo);
    setMaterialProperty(m.s0, "PoissonRatio", po);
    setMaterialProperty(m.s0, "StressThreshold", st);
    setMaterialProperty(m.s0, "NortonExponent", no);

    setMaterialProperty(m.s1, "YoungModulus", yo);
    setMaterialProperty(m.s1, "PoissonRatio", po);
    setMaterialProperty(m.s1, "StressThreshold", st);
    setMaterialProperty(m.s1, "NortonExponent", no);
  };

  auto set_temperature = [](auto& m) {
    setExternalStateVariable(m.s0, "Temperature", 293.15);
    setExternalStateVariable(m.s1, "Temperature", 293.15);
  };

  //  set_properties_elasticity(grain, 100 /*lambda*/, 75 /*mu*/); // TODO
  set_properties_norton(grain, 8.182e9, 0.364, 100.0e6, 3.333333);
  set_properties_norton(metal, 16.182e9, 0.364, 100.0e12, 3.333333); // TODO

  set_temperature(grain);
  set_temperature(metal);




  // macroscopic strain
  // macroscopic strain
  std::vector<mfem_mgis::real> e(6, mfem_mgis::real{0});
  const int xx = 0;
  const int yy = 1;
  const int zz = 2;

  /* bar{E} = e33 *(-1/2 E1 x E1 + (-1/2) * E2 x E2 + E3 x E3)*/
  const double eps = -0.012;
  e[xx] = -0.5*eps;
  e[yy] = -0.5*eps;
  e[zz] = eps;
  problem.setMacroscopicGradientsEvolution([e](const double t) { 
      auto ret = e;
      for(auto& it : ret) it *= t;
      return ret; 
      });


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

  const double dt = p.duration / double(p.nstep);
  for (int i = 0; i < p.nstep; i++) {
    mfem_mgis::Profiler::Utils::Message("Solving: from ", i*dt, " to ", (i+1)*dt);
    auto statistics = problem.solve(i * dt, dt);
    if (!statistics.status) {
      mfem_mgis::Profiler::Utils::Message("INFO: FAILED");
    }
    if(post_processing) problem.executePostProcessings(i * dt, dt);
    problem.update();
  }
  print_memory_footprint("[End]");
  mfem_mgis::Profiler::timers::print_and_write_timers();
  return EXIT_SUCCESS;
}

