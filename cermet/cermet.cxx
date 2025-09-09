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

This simulation consists of applying a tensile load on a cermet RVE. The cermet is a polycrystal where each grain has a material ID (from 2 to Nmat – 1) and a different orientation. A metallic interface is present between the different polycrystals (material ID 1)

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
  const char *behaviourMetal = "MonoCristal_UO2";
  const char *libraryGrain = "src/libBehaviour.so";
  const char *libraryMetal = "src/libBehaviour.so";
  const char *vector_file = "vectors_5grains.txt";
  int order = 1;
  int refinement = 0;
  int post_processing = 1; // default value : activated
  int verbosity_level = 0; // default value : lower level
  int nstep = 400;
  double duration = 200.;
  int bcs_type = 0;
};


struct MaterialParameters
{
  double nueff;
  double Eeff;
  double Geff;

  void operator()(double young, double poisson, double shear)
  {
    // ===========================
    // effective moduli estimations
    // ===========================
    const double C11 = young * (poisson - 1) / (2 * poisson * poisson + poisson - 1);
    const double C12 = -(poisson * young) / (2 * poisson * poisson + poisson - 1);
    const double C44 = shear;
    // ---- Bulk Modulus ----
    const double Keff = 1. / 3. * (C11 + 2 * C12);
    // ---- Shear Modulus ----
    // Reuss & Voigt bounds
    const double aniso = 2 * C44 / (C11 - C12);
    const double GReuss = 5 * C44 / (3 + 2 * aniso);
    const double GVoigt = (3 * aniso + 2) / (5 * aniso) * C44;
    // Hashin–Shtrikman bounds
    const double Gv = (C11 - C12) / 2;
    const double beta1 = -3 * (Keff + 2 * Gv) / (5 * Gv * (3 * Keff + 4 * Gv));
    const double beta2 = -3 * (Keff + 2 * C44) / (5 * C44 * (3 * Keff + 4 * C44));
    const double GHSl = Gv + 3 * ((C44 - Gv) - 4 * beta1) / 5;
    const double GHSu = C44 + 3 * ((Gv - C44) - 6 * beta2) / 5;
    // Kroener
    double m0, mus;
    double D1, sD1, D2;
    double GK = 0.;
    do
    {
      m0 = GK;
      D1 = Keff + 2 * GK;
      sD1 = 6 * D1;
      mus = GK * (9 * Keff + 8 * GK) / sD1;
      D2 = 5 * mus + 2 * C44 + 3 * Gv;
      GK = (3 * C44 * mus + 2 * Gv * mus + 5 * Gv * C44) / D2;
    } while (fabs(GK - m0) / Keff > 1e-10);

    // Shear modulus choice
    mfem_mgis::Profiler::Utils::Message("GReuss =", GReuss);
    mfem_mgis::Profiler::Utils::Message("GVoigt =", GVoigt);
    mfem_mgis::Profiler::Utils::Message("GHSl =", GHSl);
    mfem_mgis::Profiler::Utils::Message("GHSu =", GHSu);
    mfem_mgis::Profiler::Utils::Message("GK =", GK);

    Geff = GK; // GHSu
               // Young and Poisson effective moduli
    nueff = (3 * Keff - 2 * Geff) / (2 * (3 * Keff + Geff));
    Eeff = 9 * Keff * Geff / (3 * Keff + Geff);
  }
}; // struct MaterialParameters

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
  args.AddOption(&p.bcs_type, "-bcs", "--bcs_type",
      "Types of boundary conditions. For all BCs, a displacement of def = 5e-4 m/s is imposed. Type 0: zero strain in xx and yy. Type 1: imposed displacement of -0.3 * def in xx and yy.");
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
void setup_material_properties(ProblemT& problem, TestParameters& p, MaterialParameters& mp)
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

  mp(young1, poisson12, shear12); // used in solve_null_strain

  // ceramic
  auto set_temperature = [](auto& m) {
    setExternalStateVariable(m.s0, "Temperature", 1600);
    setExternalStateVariable(m.s1, "Temperature", 1600);
  };

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
    throw std::invalid_argument("setup_properties : incorrect number of vectors in vector file.");
  }

  // TEST metal
  {
    problem.addBehaviourIntegrator("Mechanics", 1, p.libraryMetal, p.behaviourMetal);
    auto& metal = problem.getMaterial(1);
    std::array<mfem_mgis::MaterialAxis3D, 2u> r;
    r[0] = std::array<mfem_mgis::real, 3u>{0.7071067811865475, -0.4086070447619255, -0.5770964243269279};
    r[1] = std::array<mfem_mgis::real, 3u>{0.7071067811865475, 0.4086070447619256, 0.5770964243269281};
    metal.setRotationMatrix(mfem_mgis::RotationMatrix3D{r});
    set_properties_mono(metal, young1, young2, young3, poisson12, poisson23, poisson13, shear12, shear23, shear13);
    set_temperature(metal);
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
      r[0] = vectors.at(2*(grainID - 2));
      r[1] = vectors.at(2*(grainID - 2) + 1);
      grain.setRotationMatrix(mfem_mgis::RotationMatrix3D{r});
    }
  }
}

// Version 1
  template<typename Problem>
void solve_impose_displ(Problem& problem, double dt, int nstep, bool post_processing)
{
  // Traction
  const double def = 5e-4; // 0.01;
  problem.setMacroscopicGradientsEvolution([def](const double t) {
      const int xx = 0;
      const int yy = 1;
      const int zz = 2;
      auto ret = std::vector<mfem_mgis::real>(9, mfem_mgis::real{});
      ret[xx] = 1 - 0.3 * def * t;
      ret[yy] = 1 - 0.3 * def * t;
      ret[zz] = 1 + def * t;
      return ret; });

  for (int i = 0; i < nstep; i++) {
    mfem_mgis::Profiler::Utils::Message("Solving: from ", i*dt, " to ", (i+1)*dt);
    auto statistics = problem.solve(i * dt, dt);
    if (!statistics.status) {
      mfem_mgis::Profiler::Utils::Message("INFO: FAILED");
    }
    if(post_processing) problem.executePostProcessings(i * dt, dt);
    problem.update();
    print_memory_footprint("[At timestep: " + std::to_string(i) + "]");
  }
}

// Version 2
  template<typename Problem>
void solve_null_strain(Problem& problem, double dt, int nstep, bool post_processing, MaterialParameters& mp)
{
  // Traction
  const double def = 5e-4;

  double Fxx, Fyy, Fzz;
  // init macro cauchy stress components
  double Sxx = 0.;
  double Syy = 0.;

  // --- use for MPI Reductions --- //
  double SIG_local[2];
  double SIG_global[2];

  // --- variables used to loop to compute Fyy and Fxx --- //
  // --- fixed-point param
  const double tolFP = 1.e4;
  const int maxitFP = 50;

  // warning Fxx, Fyy, and Fzz are defined during the following loop
  problem.setMacroscopicGradientsEvolution([&Fxx, &Fyy, &Fzz](const double t)
      { 
      auto ret = std::vector<mfem_mgis::real>(9, mfem_mgis::real{});
      ret[0] = Fxx;
      ret[1] = Fyy;
      ret[2] = Fzz;
      return ret; });


  for (int i = 0; i < nstep; i++) 
  {
    mfem_mgis::Profiler::Utils::Message("Solving: from ", i*dt, " to ", (i+1)*dt);

    // -- update F from input data, here def = 5e-4. -- //
    Fzz = 1. + def * (i + 1.) * dt; // / end);


    if (i == 0) 
    {
      Fxx = 1.;// / std::sqrt(Fzz);
      Fyy = 1.;// / std::sqrt(Fzz);
    }
    // try another start for fixed point
    //Fxx =  1.;
    //Fyy =  1.;
    Sxx = 0.0;
    Syy = 0.0;
    // --- reset fixed point variables --- //
    int itFP = 0;
    double newRes = 0.0;
    double oldRes = 1.e6;

    // --- fixed point loop --- // 
    while ((abs(oldRes - newRes) > tolFP) && (itFP <= maxitFP))
    {
      // -- Debug -- //
      // -- mfem_mgis::Profiler::Utils::Message("debug: Sxx -> ", Sxx, " Syy -> ", Syy );
      // -- mfem_mgis::Profiler::Utils::Message("debug:: Eeff -> ", mp.Eeff, " mp.nueff -> ", mp.nueff );

      double Fcorr[2];
      // -- compute correction -- //
      // -- note that Sxx and Syy are computed at the end of this loop section, default is 0 -- //
      Fcorr[0] = ((mp.nueff*mp.nueff-1)*Sxx)/mp.Eeff+(mp.nueff*(mp.nueff+1)*Syy)/mp.Eeff;
      Fcorr[1] = ((mp.nueff*mp.nueff-1)*Syy)/mp.Eeff+(mp.nueff*(mp.nueff+1)*Sxx)/mp.Eeff;

      // -- Increment F with F correction -- //
      Fxx += Fcorr[0];
      Fyy += Fcorr[1];

      // -- Debug -- //
      // -- mfem_mgis::Profiler::Utils::Message("debug: ", Fxx, Fyy, Fcorr[0], Fcorr[1], Fzz);

      // -- note that MacroscopicGradientsEvolution taktes into account of new Fxx and Fyy -- //
      auto statistics = problem.solve(i * dt, dt);
      if (!statistics.status) {
        mfem_mgis::Profiler::Utils::Message("The Linear Solver has failed.");
        //std::exit(EXIT_FAILURE);
      }
      problem.update();

      // -- Get stress and volume for all mateirals -- //
      auto [tf_integrals, volumes] = mfem_mgis::computeMeanThermodynamicForcesValues<true>(problem.template getImplementation<true>());

      // -- compute volume on all Materials Identifiers -- //
      // -- compute partial contribution to Sxx and Syy
      double volume = 0.0;
      Sxx = 0.0;
      Syy = 0.0;
      for (const auto &m : problem.getAssignedMaterialsIdentifiers())
      {
        volume += volumes[m];
        Sxx += tf_integrals[m][0];
        Syy += tf_integrals[m][1];
      }
      assert(volume > 0.0);

      // -- Get volume from other subdomains -- //
      // -- local volume -> global volume
      MPI_Allreduce(MPI_IN_PLACE, &volume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

      // -- compute mean value of Thermodynamic Forces on all Materials Identifiers -- //
      Sxx /= volume;
      Syy /= volume;

      // -- PKI to VM -- //
      Sxx = Sxx*Fxx/(Fxx*Fyy*Fzz);
      Syy = Syy*Fyy/(Fxx*Fyy*Fzz);

      // -- get Sxx and Syy from other subdomains -- //
      SIG_local[0] = Sxx;
      SIG_local[1] = Syy;
      MPI_Allreduce(SIG_local, SIG_global, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      Sxx = SIG_global[0];
      Syy = SIG_global[1];

      oldRes = newRes;
      newRes = std::sqrt(Sxx * Sxx + Syy * Syy);

      mfem_mgis::Profiler::Utils::Message("Fixed Point iteration", itFP, ": |res| =", abs(oldRes - newRes));
      if (abs(oldRes - newRes) > tolFP) { problem.revert(); }      
      itFP++;
    }

    if (itFP >= maxitFP) mfem_mgis::Profiler::Utils::Message("warning: maximum number of iterations for the fixed-point algorithm attained, before the requested tolerance is reached");

    if(post_processing) problem.executePostProcessings(i * dt, dt);
    print_memory_footprint("[At timestep: " + std::to_string(i) + "]");
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
      {"RelativeTolerance", 1e-8},
      {"AbsoluteTolerance", 0.},
      {"MaximumNumberOfIterations", 15}});
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
    mfem_mgis::Parameters{{"Name", "HypreBoomerAMG"}, {"Options",  mfem_mgis::Parameters{{"VerbosityLevel", verbosity}}}};
  solverParameters.insert(
      mfem_mgis::Parameters{{"Preconditioner", preconditionner}});
  problem.setLinearSolver("HyprePCG", solverParameters);
  //problem.setLinearSolver("HypreGMRES", solverParameters);
  //problem.setLinearSolver("MUMPSSolver", {});


  // --- Define material porperties directly into setup_material_properties --- //
  // --- Initialize "MaterialParameters mp" that is required for "solve_null_strain" --- // 
  MaterialParameters mp;
  setup_material_properties(problem, p, mp);

  if (post_processing) 
  {
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
  if(p.bcs_type == 0) solve_null_strain(problem, dt, p.nstep, post_processing, mp);
  if(p.bcs_type == 1) solve_impose_displ(problem, dt, p.nstep, post_processing);
  if(p.bcs_type > 1) std::exit(EXIT_FAILURE); 

  print_memory_footprint("[End]");
  mfem_mgis::Profiler::timers::print_and_write_timers();
  return EXIT_SUCCESS;
}
