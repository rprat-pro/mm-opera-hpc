#pragma once

#include "MFEMMGIS/NonLinearEvolutionProblemImplementation.hxx"
#include "MFEMMGIS/NonLinearEvolutionProblem.hxx"
#include <MFEMMGIS/Profiler.hxx>


  template<typename Implementation>
void print_mesh_information(Implementation& impl)
{

  using mfem_mgis::Profiler::Utils::sum;
  using mfem_mgis::Profiler::Utils::Message;
  Message("INFO: print_mesh_information");

  //getMesh
  auto mesh = impl.getFiniteElementSpace().GetMesh();

  //get the number of vertices
  int64_t numbers_of_vertices_local = mesh->GetNV();
  int64_t  numbers_of_vertices = sum(numbers_of_vertices_local);

  //get the number of elements
  int64_t numbers_of_elements_local = mesh->GetNE();
  int64_t numbers_of_elements = sum(numbers_of_elements_local);

  //get the element size
  double h = mesh->GetElementSize(0);

  // get n dofs
  auto& fespace = impl.getFiniteElementSpace();
  int64_t unknowns_local = fespace.GetTrueVSize();
  int64_t unknowns = sum(unknowns_local);

  Message("INFO: number of vertices -> ", numbers_of_vertices);
  Message("INFO: number of elements -> ", numbers_of_elements);
  Message("INFO: element size -> ", h);
  Message("INFO: Number of finite element unknowns: " , unknowns);
}


