/*!
 * \file   src/FirstPrincipalStressLocation.cxx
 * \brief
 * \author Thomas Helfer
 * \date   29/09/2024
 */

#include <array>
#include <vector>
#include <limits>
#include <algorithm>
#include <TFEL/Math/stensor.hxx>
#include "mfem/linalg/densemat.hpp"
#include "mfem/fem/fespace.hpp"
#ifdef MFEM_USE_MPI
#include "mfem/fem/pfespace.hpp"
#endif /* MFEM_USE_MPI */
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/Material.hxx"
#include "OperaHPC/Utilities.hxx"

namespace opera_hpc {

  template <bool parallel>
  FirstPrincipalStressValueAndLocation
  findFirstPrincipalStressValueAndLocationImplementation(
      const mfem_mgis::Material &m) {
    constexpr auto stress_size = mfem_mgis::size_type{6};
    const auto &s = m.getPartialQuadratureSpace();
    const auto &fed = s.getFiniteElementDiscretization();
    const auto &fespace = fed.getFiniteElementSpace<parallel>();
    // material identifier
    const auto mid = s.getId();
    // values of the stress at the end of the time step
    const auto *stress_values = m.s1.thermodynamic_forces.data();
    auto max_stress = -std::numeric_limits<mfem_mgis::real>::max();
    auto max_stress_location = std::array<mfem_mgis::real, 3u>{};
    auto tmp = mfem::Vector{};
    // loop over the elements
    for (mfem_mgis::size_type i = 0; i != fespace.GetNE(); ++i) {
      if (fespace.GetAttribute(i) != mid) {
        // the element is not associated with the considered material
        continue;
      }
      const auto &fe = *(fespace.GetFE(i));
      auto &tr = *(fespace.GetElementTransformation(i));
      const auto &ir = s.getIntegrationRule(fe, tr);
      // offset of the element
      const auto eo = s.getOffset(i);
      // loop over the integration points
      for (mfem_mgis::size_type g = 0; g != ir.GetNPoints(); ++g) {
        const auto *const stress = stress_values + (eo + g) * stress_size;
        auto sig = tfel::math::stensor<3u, mfem_mgis::real>{stress};
        const auto sig_vp = sig.computeEigenValues<
            tfel::math::stensor_common::FSESJACOBIEIGENSOLVER>();
        const auto mvp = *(tfel::fsalgo::max_element<3>::exe(sig_vp.begin()));
        if (mvp > max_stress) {
          max_stress = mvp;
          // update the location
          const auto &ip = ir.IntPoint(g);
          tr.SetIntPoint(&ip);
          tr.Transform(tr.GetIntPoint(), tmp);
          max_stress_location = {tmp[0], tmp[1], tmp[2]};
        }
      }
    }
    if constexpr (parallel) {
      int nprocs;
      MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
      auto all_max_stresses = std::vector<double>(nprocs);
      auto all_locations = std::vector<double>(3 * nprocs);
      MPI_Allgather(&max_stress, 1, MPI_DOUBLE, all_max_stresses.data(), 1,
                    MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Allgather(max_stress_location.data(), 3, MPI_DOUBLE,
                    all_locations.data(), 3, MPI_DOUBLE, MPI_COMM_WORLD);
      const auto i =
          std::max_element(all_max_stresses.begin(), all_max_stresses.end()) -
          all_max_stresses.begin();
      max_stress = all_max_stresses[i];
      std::copy(all_locations.begin() + 3 * i,
                all_locations.begin() + 3 * (i + 1),
                max_stress_location.begin());
    }
    return {.value = max_stress, .location = max_stress_location};
  }  // end of findFirstPrincipalStressValueAndLocationImplementation

  FirstPrincipalStressValueAndLocation findFirstPrincipalStressValueAndLocation(
      const mfem_mgis::Material &m) {
    const auto &s = m.getPartialQuadratureSpace();
    const auto &fed = s.getFiniteElementDiscretization();
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      return findFirstPrincipalStressValueAndLocationImplementation<true>(m);
#else  /* MFEM_USE_MPI */
      mfem_mgsi::reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
    }
    return findFirstPrincipalStressValueAndLocationImplementation<false>(m);
  }  // end of findFirstPrincipalStressValueAndLocation

}  // end of namespace opera_hpc
