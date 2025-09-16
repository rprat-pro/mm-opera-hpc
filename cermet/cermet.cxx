#include <cstdlib>
#include <optional>
#include <functional>
#include <sys/resource.h>

#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/ParaviewExportIntegrationPointResultsAtNodes.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"
#include "MFEMMGIS/Profiler.hxx"

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

struct MaterialParameters {
  double nueff;
  double Eeff;
  double Geff;

  void operator()(double young, double poisson, double shear) {
    // ===========================
    // effective moduli estimations
    // ===========================
    const double C11 =
        young * (poisson - 1) / (2 * poisson * poisson + poisson - 1);
    const double C12 =
        -(poisson * young) / (2 * poisson * poisson + poisson - 1);
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
    const double beta2 =
        -3 * (Keff + 2 * C44) / (5 * C44 * (3 * Keff + 4 * C44));
    const double GHSl = Gv + 3 * ((C44 - Gv) - 4 * beta1) / 5;
    const double GHSu = C44 + 3 * ((Gv - C44) - 6 * beta2) / 5;
    // Kroener
    double m0, mus;
    double D1, sD1, D2;
    double GK = 0.;
    do {
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

// display information
static long get_memory_checkpoint() {
  rusage obj;
  int who = 0;
  [[maybe_unused]] auto test = getrusage(who, &obj);
  assert((test = -1) && "error: getrusage has failed");
  long res;
  MPI_Reduce(&(obj.ru_maxrss), &(res), 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  return res;
};

static void print_memory_footprint(std::string msg) {
  long mem = get_memory_checkpoint();
  double m = double(mem) * 1e-6; // conversion kb to Gb
  mfem_mgis::Profiler::Utils::Message(msg, " memory footprint: ", m, " GB");
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

[[nodiscard]] static mfem_mgis::real
norm(const std::array<mfem_mgis::real, 3u> &u) {
  mfem_mgis::real ret = 0.;
  for (auto val : u) {
    ret += val * val;
  }
  return (std::sqrt(ret));
}

[[nodiscard]] static std::vector<std::array<mfem_mgis::real, 3u>>
readVectorsFromFile(const std::string &filename) {
  // TODO: check if vectors size is equal to nMat
  std::vector<std::array<mfem_mgis::real, 3u>> vectors;
  std::ifstream inputFile(filename);

  auto checknorm = [](std::array<mfem_mgis::real, 3u> v) {
    mfem_mgis::real normv = norm(v);
    if (normv < 1.e-15) {
      throw std::invalid_argument("checknorm : vectors must not be null");
    }
    if (abs(normv - 1.) > 1.e-15) {
      std::cout << "checknorm : normalizing vector..." << std::endl;
      for (auto &comp : v) {
        comp *= 1. / normv;
      }
    }
    return v;
  };

  if (!inputFile.is_open()) {
    std::cout << "Erreur lors de l'ouverture du fichier " << filename
              << std::endl;
    return vectors;
  }

  std::string line;
  while (std::getline(inputFile, line)) {
    std::istringstream iss(line);
    std::array<mfem_mgis::real, 3u> vector;
    if (!(iss >> vector[0] >> vector[1] >> vector[2])) {
      std::cout << "Erreur de lecture du vecteur : " << line << std::endl;
      continue;
    }
    vectors.push_back(checknorm(vector));
  }

  inputFile.close();
  return vectors;
}

// add postprocessing
static void post_process(mfem_mgis::PeriodicNonLinearEvolutionProblem &p,
                         double start, double end) {
  CatchTimeSection("common::post_processing_step");
  p.executePostProcessings(start, end);
}

static void
setup_material_properties(mfem_mgis::PeriodicNonLinearEvolutionProblem &problem,
                          TestParameters &p, MaterialParameters &mp) {

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

  std::vector<std::array<mfem_mgis::real, 3u>> vectors =
      readVectorsFromFile(p.vector_file);
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

[[nodiscard]] static const std::array<mfem_mgis::real, 3u>
computeMacroscopicCauchyStress(
    mfem_mgis::PeriodicNonLinearEvolutionProblem &problem,
    const std::array<mfem_mgis::real, 3u> &F) {
  // integrals of the diagonal components of the First Piola-Kirchhoff stress
  auto pk1_integral = std::array<mfem_mgis::real, 3u>{};
  auto volume = mfem_mgis::real{};
  // get pk1 integral and volume on all materials //
  auto [pk1_integrals, volumes] =
      mfem_mgis::computeMeanThermodynamicForcesValues<true>(
          problem.template getImplementation<true>());
  // sum volume of all materials
  // sum pk1 integrals
  for (const auto &m : problem.getAssignedMaterialsIdentifiers()) {
    volume += volumes[m];
    pk1_integral[0] += pk1_integrals[m][0];
    pk1_integral[1] += pk1_integrals[m][1];
    pk1_integral[2] += pk1_integrals[m][2];
  }
  assert(volume > 0.0);
  // -- sum volumes of all processess
  MPI_Allreduce(MPI_IN_PLACE, &volume, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  // -- sum integrals from all processes
  MPI_Allreduce(MPI_IN_PLACE, pk1_integral.data(), 3, MPI_DOUBLE, MPI_SUM,
                MPI_COMM_WORLD);
  // macroscopic values of the First Piola-Kirchhoff stress
  const auto pk1 = std::array<mfem_mgis::real, 3u>{pk1_integral[0] / volume,
                                                   pk1_integral[1] / volume,
                                                   pk1_integral[2] / volume};
  // macroscopic change of volume
  const auto det = F[0] * F[1] * F[2];
  // computation of the diagonal components of the Cauchy stress
  return {pk1[0] * F[0] / det, pk1[1] * F[1] / det, pk1[2] * F[2] / det};
} // end of computeMacroscopicCauchyStress

struct MacroscopicVariables {
  // diagonal part of the deformation gradient
  std::array<mfem_mgis::real, 3u> F = std::array<mfem_mgis::real, 3u>{1, 1, 1};
  // increment of the diagonal part of the deformation gradient during the
  // previous time step
  std::array<mfem_mgis::real, 3u> dF = std::array<mfem_mgis::real, 3u>{};
  // diagonal components of the macroscropic Cauchy stress
  std::array<mfem_mgis::real, 3u> S = std::array<mfem_mgis::real, 3u>{};
  // previous succesful time step
  std::optional<mfem_mgis::real> previous_time_increment;
};

[[nodiscard]] static bool simulateOverATimeStep(
    mfem_mgis::PeriodicNonLinearEvolutionProblem &problem,
    MacroscopicVariables &macroscopic_variables, //
    const bool post_processing, const MaterialParameters &mp,
    const mfem_mgis::real bts, //
    const mfem_mgis::real ets) {
  using namespace mfem_mgis::Profiler::Utils;
  Message("Solving time step from ", bts, " to ", ets);
  // --- fixed-point param
  const double tolFP = 1.e4;
  const int maxitFP = 50;
  // Traction
  const double def = 5e-4;
  const auto dt = ets - bts;

  auto &F = macroscopic_variables.F;
  auto &S = macroscopic_variables.S;
  // deformation gradient at the beginning of the time step
  auto S0 = S;
  auto F0 = F;

  if (macroscopic_variables.previous_time_increment.has_value()) {
    auto &dF = macroscopic_variables.dF;
    const auto pdt = *(macroscopic_variables.previous_time_increment);
    // -- update F from input data, here def = 5e-4. -- //
    F[0] += dF[0] * (dt / pdt);
    F[1] += dF[1] * (dt / pdt);
  }
  F[2] = 1. + def * ets; // / end);
  // --- fixed point iteration number
  int itFP = 0;
  // --- fixed point loop --- //
  bool converged = false;
  while ((!converged) && (itFP <= maxitFP)) {
    // -- compute correction -- //
    // -- note that S[0] and S[1] are updated at the end of this loop section,
    // default is 0 -- //
    F[0] += ((mp.nueff * mp.nueff - 1) * S[0]) / mp.Eeff +
            (mp.nueff * (mp.nueff + 1) * S[1]) / mp.Eeff;
    F[1] += ((mp.nueff * mp.nueff - 1) * S[1]) / mp.Eeff +
            (mp.nueff * (mp.nueff + 1) * S[0]) / mp.Eeff;
    // -- note that MacroscopicGradientsEvolution taktes into account of new
    // F[0] and F[1] -- //
    auto statistics = problem.solve(bts, dt);
    if (!statistics.status) {
      Message("error: the resolution failed.");
      S = S0;
      F = F0;
      problem.revert();
      return false;
    }

    S = computeMacroscopicCauchyStress(problem, F);
    const auto r = std::sqrt(S[0] * S[0] + S[1] * S[1]);
    //
    Message("Fixed Point iteration", itFP, ": |res| =", r);
    converged = std::abs(r) < tolFP;
    itFP++;
  }
  //
  Message("Solution at time", ets, ":", F[2], S[0], S[1], S[2]);
  //
  if (itFP >= maxitFP) {
    Message("error: maximum number of iterations for the "
            "fixed-point algorithm reached");
    S = S0;
    F = F0;
    problem.revert();
    return false;
  }
  //
  if (post_processing) {
    problem.executePostProcessings(bts, dt);
  }
  // update state variable for the next time step
  problem.update();
  // update information for time extrapolation
  auto &dF = macroscopic_variables.dF;
  macroscopic_variables.previous_time_increment = dt;
  for (std::size_t i = 0; i != 3; ++i) {
    dF[i] = F[i] - F0[i];
  }
  return true;
}

[[nodiscard]] static bool simulateOverATemporalSequence(
    mfem_mgis::PeriodicNonLinearEvolutionProblem &problem,
    MacroscopicVariables &macroscopic_variables, //
    const bool post_processing, const MaterialParameters &mp,
    const mfem_mgis::real bts, //
    const mfem_mgis::real ets) {
  using namespace mfem_mgis::Profiler::Utils;
  constexpr auto maximum_of_substeps = std::size_t{10};
  auto nsubsteps = std::size_t{};
  auto nstep = std::size_t{1};
  auto t = bts;
  auto dt = ets - bts;
  while (nstep != 0) {
    const auto do_post_processing = post_processing && (nstep == 1);
    const auto success = simulateOverATimeStep(
        problem, macroscopic_variables, do_post_processing, mp, t, t + dt);
    if (success) {
      t += dt;
      --nstep;
    } else {
      ++nsubsteps;
      if (nsubsteps == maximum_of_substeps) {
        Message("Maximum number of substeps reached for temporal sequence from",
                bts, "to", ets);
        return false;
      } else {
        Message("dividing time step by 2");
      }
      dt *= mfem_mgis::real{0.5};
      nstep *= 2;
    }
  }
  //
  print_memory_footprint("[At end of temporal sequence: " + std::to_string(ets) + "]");
  return true;
} // end of simulateOverATemporalSequence

[[nodiscard]] static bool
solve_null_strain(mfem_mgis::PeriodicNonLinearEvolutionProblem &problem,
                  const double te, const std::size_t nsteps,
                  const bool post_processing, const MaterialParameters &mp) {
  //
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
  MacroscopicVariables macroscopic_variables;

  // warning Fxx, F[1], and Fzz are defined during the following loop
  problem.setMacroscopicGradientsEvolution(
      [&macroscopic_variables](const double t) {
        auto &F = macroscopic_variables.F;
        auto ret = std::vector<mfem_mgis::real>(9, mfem_mgis::real{});
        std::copy(F.begin(), F.end(), ret.begin());
        return ret;
      });

  for (std::size_t i = 0; i < temporal_sequences.size(); i++) {
    if (!simulateOverATemporalSequence(
            problem, macroscopic_variables, post_processing, mp,
            temporal_sequences[i], temporal_sequences[i + 1])) {
      return false;
    }
  }
  return true;
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
  print_memory_footprint("[Building problem]");

  // non linear solver
  problem.setSolverParameters({{"VerbosityLevel", 1},
                               {"RelativeTolerance", 1e-6},
                               //      {"AbsoluteTolerance", 0.},
                               {"AbsoluteTolerance", 1e-6},
                               {"MaximumNumberOfIterations", 20}});

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

  // --- Define material porperties directly into setup_material_properties ---
  // //
  // --- Initialize "MaterialParameters mp" that is required for
  // "solve_null_strain" --- //
  MaterialParameters mp;
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

  const auto success =
      solve_null_strain(problem, p.duration, p.nsteps, post_processing, mp);

  print_memory_footprint("[End]");
  mfem_mgis::Profiler::timers::print_and_write_timers();
  return success ? EXIT_SUCCESS : EXIT_FAILURE;
}
