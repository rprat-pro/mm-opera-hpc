#include <cstdlib>
#include <sys/resource.h>
#include <functional>

#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/UniformDirichletBoundaryCondition.hxx"
#include "MFEMMGIS/ParaviewExportIntegrationPointResultsAtNodes.hxx"


/*

  Problem:

  Definition:

  - contact law

  - Parameters:

  material: [grain, metal]

  element: H1
  order: 2

 */
struct TestParameters {
  const char *mesh_file = "mesh/5grains.msh";
  const char *behaviourGrain = "MonoCristal_UO2";
  const char *behaviourMetal = "NortonCr";
  const char *libraryGrain = "src/libBehaviour.so";
  const char *libraryMetal = "src/libBehaviour.so";
  const char *vector_file = "vectors_5grains.txt";
  int order = 2;
  int refinement = 0;
  int post_processing = 1; // default value : activated
  int verbosity_level = 0; // default value : lower level
  int nstep = 40;
  double duration = 1.;
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
  args.AddOption(&p.vector_file, "-f", "--file",
      "Vector file to use");
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

auto norm(std::array<mfem_mgis::real, 3u> &u)
{
  mfem_mgis::real ret = 0.;
  for (auto val : u)
  {
    ret += val * val;
  }
  return (std::sqrt(ret));
}

std::vector<std::array<mfem_mgis::real, 3u>> readVectorsFromFile(const std::string &filename)
{
  // TODO: check if vectors size is equal to nMat
  std::vector<std::array<mfem_mgis::real, 3u>> vectors;
  std::ifstream inputFile(filename);

  auto checknorm = [](std::array<mfem_mgis::real, 3u> v)
  {
    mfem_mgis::real normv = norm(v);
    if (normv < 1.e-15)
    {
      throw std::invalid_argument("checknorm : vectors must not be null");
    }
    if (abs(normv - 1.) > 1.e-15)
    {
      std::cout << "checknorm : normalizing vector..." << std::endl;
      for (auto &comp : v)
      {
        comp *= 1. / normv;
      }
    }
    return v;
  };

  if (!inputFile.is_open())
  {
    std::cout << "Erreur lors de l'ouverture du fichier " << filename << std::endl;
    return vectors;
  }

  std::string line;
  while (std::getline(inputFile, line))
  {
    std::istringstream iss(line);
    std::array<mfem_mgis::real, 3u> vector;
    if (!(iss >> vector[0] >> vector[1] >> vector[2]))
    {
      std::cout << "Erreur de lecture du vecteur : " << line << std::endl;
      continue;
    }
    vectors.push_back(checknorm(vector));
  }

  inputFile.close();
  return vectors;
}

// add postprocessing
template <typename Problem>
void post_process(Problem &p, double start, double end) {
  CatchTimeSection("common::post_processing_step");
  p.executePostProcessings(start, end);
}

  template<typename ProblemT>
void setup_material_properties(ProblemT& problem, TestParameters& p)
{

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

  // ceramic
  problem.addBehaviourIntegrator("Mechanics", 1, p.libraryMetal, p.behaviourMetal);
  auto& metal = problem.getMaterial(1);

  auto set_properties_norton = [](auto &m,
      const double young, const double poisson)
  {
    setMaterialProperty(m.s0, "YoungModulus", young);
    setMaterialProperty(m.s0, "PoissonRatio", poisson);
    setMaterialProperty(m.s1, "YoungModulus", young);
    setMaterialProperty(m.s1, "PoissonRatio", poisson);
  };

  auto set_temperature = [](auto& m) {
    setExternalStateVariable(m.s0, "Temperature", 1600);
    setExternalStateVariable(m.s1, "Temperature", 1600);
  };

  set_properties_norton(metal, young1, poisson12);
  set_temperature(metal);


  auto set_properties_mono = [](auto &m,
      const double yo1, const double yo2, const double yo3,
      const double po12, const double po23, const double po13,
      const double sm12, const double sm23, const double sm13)
  {
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

  const int nMat = getMaterialsAttributes(*(problem.getFiniteElementDiscretizationPointer())).Max();

  std::vector<std::array<mfem_mgis::real, 3u>> vectors = readVectorsFromFile(p.vector_file);
  if (vectors.size()!=2*(nMat-1)) {
   std::cout << vectors.size() << "!=" << 2*(nMat-1) << std::endl; 
   throw std::invalid_argument("setup_properties : incorrect number of vectors in vector file");
  }

  for(int grainID = 2 ; grainID <= nMat ; grainID++)
  {
    problem.addBehaviourIntegrator("Mechanics", grainID, p.libraryGrain, p.behaviourGrain);
    auto& grain = problem.getMaterial(grainID);
    set_properties_mono(grain, young1, young2, young3, poisson12, poisson23, poisson13, shear12, shear23, shear13);
    set_temperature(grain);

    std::array<mfem_mgis::MaterialAxis3D, 2u> r;
    if (grain.b.symmetry == mgis::behaviour::Behaviour::ORTHOTROPIC)
    {
      r[0] = vectors[2 * (grainID - 1)];
      r[1] = vectors[2 * (grainID - 1) + 1 ];
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
  solverParameters.insert(
      mfem_mgis::Parameters{{"Preconditioner", preconditionner}});
  problem.setLinearSolver("HyprePCG", solverParameters);

  setup_material_properties(problem, p);

  const double def = -0.01;
  problem.setMacroscopicGradientsEvolution([def](const double t) {
      const int xx = 0;
      const int yy = 1;
      const int zz = 2;
      auto ret = std::vector<mfem_mgis::real>(9, mfem_mgis::real{});
      ret[xx] = 1 - 0.5 * def;
      ret[yy] = 1 - 0.5 * def;
      ret[zz] = 1 + def;
      return ret; });

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
    print_memory_footprint("[At timestep: " + std::to_string(i) + "]");
  }
  print_memory_footprint("[End]");
  mfem_mgis::Profiler::timers::print_and_write_timers();
  return EXIT_SUCCESS;
}
