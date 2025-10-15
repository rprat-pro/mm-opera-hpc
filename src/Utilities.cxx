/*!
 * \file   Utilities.cxx
 * \brief
 * \author Thomas Helfer
 * \date   21/09/2025
 */

#include <sys/resource.h>
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/PeriodicNonLinearEvolutionProblem.hxx"
#include "MM_OPERA_HPC/Utilities.hxx"

namespace mm_opera_hpc {

  void printMeshInformation(
      mfem_mgis::PeriodicNonLinearEvolutionProblem &problem) {
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
    Message("Info problem: element size -> ", h);
    Message("Info problem: number of finite element unknowns: ", unknowns);
  }

  // display information
  static long get_memory_checkpoint() {
    rusage obj;
    int who = 0;
    [[maybe_unused]] auto test = getrusage(who, &obj);
    assert((test = -1) && "error: getrusage has failed");
    long res;
    MPI_Reduce(&(obj.ru_maxrss), &(res), 1, MPI_LONG, MPI_SUM, 0,
               MPI_COMM_WORLD);
    return res;
  };

  void printMemoryFootprint(std::string_view msg) noexcept {
    long mem = get_memory_checkpoint();
    double m = double(mem) * 1e-6;  // conversion kb to Gb
    mfem_mgis::Profiler::Utils::Message(msg, " memory footprint: ", m, " GB");
  }

}  // end of namespace mm_opera_hpc
