/*!
 * \file   src/Utilities.cxx
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
#include "MFEMMGIS/Profiler.hxx"
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/Material.hxx"
#include "OperaHPC/Utilities.hxx"

namespace opera_hpc {

template <bool parallel>
static FirstPrincipalStressValueAndLocation
findFirstPrincipalStressValueAndLocationImplementation(
    const mfem_mgis::Material &m) {
  CatchTimeSection("OperaHPC::FindFirstPrincipal");
  constexpr mfem_mgis::size_type stress_size = 6;
  const auto &s = m.getPartialQuadratureSpace();
  const auto &fed = s.getFiniteElementDiscretization();
  const auto &fespace = fed.getFiniteElementSpace<parallel>();
  // material identifier
  const auto mid = s.getId();
  // values of the stress at the end of the time step
  const auto *stress_values = m.s1.thermodynamic_forces.data();
  auto max_stress = -std::numeric_limits<mfem_mgis::real>::max();
  std::array<mfem_mgis::real, 3u> max_stress_location;
  mfem::Vector tmp;
  // loop over the elements
  {
    CatchTimeSection("FindFirstPrincipal::LoopOverElem");
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
  }
  if constexpr (parallel) {
    CatchTimeSection("FindFirstPrincipal::GatherResults");
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
              all_locations.begin() + 3 * (i + 1), max_stress_location.begin());
  }
  return {.value = max_stress, .location = max_stress_location};
} // end of findFirstPrincipalStressValueAndLocationImplementation

template <bool parallel>
static std::vector<std::array<mfem_mgis::real, 3u>>
getPointsAboveStressThresholdImplementation(const mfem_mgis::Material &m,
                                            const mfem_mgis::real threshold) {
  CatchTimeSection("OperaHPC::GetPointsAbove");
  constexpr mfem_mgis::size_type stress_size = 6;
  const auto &s = m.getPartialQuadratureSpace();
  const auto &fed = s.getFiniteElementDiscretization();
  const auto &fespace = fed.getFiniteElementSpace<parallel>();
  // material identifier
  const auto mid = s.getId();
  // values of the stress at the end of the time step
  const auto *stress_values = m.s1.thermodynamic_forces.data();
  std::vector<std::array<mfem_mgis::real, 3u>> local_locations;
  // loop over the elements
  mfem::Vector tmp;
  {
    CatchTimeSection("GetPointsAbove::LoopOverElem");
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
        if (mvp > threshold) {
          const auto &ip = ir.IntPoint(g);
          tr.SetIntPoint(&ip);
          tr.Transform(tr.GetIntPoint(), tmp);
          local_locations.push_back({tmp[0], tmp[1], tmp[2]});
        }
      }
    }
  }
  if constexpr (parallel) {
    CatchTimeSection("GetPointsAbove::GatherResults");
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    const auto lsize = static_cast<int>(local_locations.size());
    std::vector<double> local_locations_tmp;
    local_locations_tmp.resize(3 * lsize);
    for (int idx = 0; idx != lsize; ++idx) {
      local_locations_tmp[3 * idx] = local_locations[idx][0];
      local_locations_tmp[3 * idx + 1] = local_locations[idx][1];
      local_locations_tmp[3 * idx + 2] = local_locations[idx][2];
      // std::cout << local_locations[idx][0] << "\t" << local_locations[idx][1]
      // << "\t" << local_locations[idx][2] << std::endl;
    }
    int gsize;
    MPI_Allreduce(&lsize, &gsize, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
    std::vector<double> all_locations_tmp(3 * gsize);
    std::vector<int> counts_recv(nprocs);
    std::vector<int> displacements(nprocs);

    MPI_Allgather(&lsize, 1, MPI_INT, counts_recv.data(), 1, MPI_INT,
                  MPI_COMM_WORLD);

    displacements[0] = 0;
    for (int i = 1; i < nprocs; ++i) {
      displacements[i] = displacements[i - 1] + counts_recv[i - 1] * 3;
    }

    for (int i = 0; i < nprocs; ++i) {
      counts_recv[i] *= 3;
    }

    MPI_Allgatherv(local_locations_tmp.data(), 3 * lsize, MPI_DOUBLE,
                   all_locations_tmp.data(), counts_recv.data(),
                   displacements.data(), MPI_DOUBLE, MPI_COMM_WORLD);
    std::vector<std::array<mfem_mgis::real, 3u>> all_locations;
    all_locations.resize(gsize);

    for (auto &locs : all_locations) {
      std::fill(locs.begin(), locs.end(), -1e6);
    }

    for (int idx = 0; idx != gsize; ++idx) {
      all_locations[idx][0] = all_locations_tmp[3 * idx];
      all_locations[idx][1] = all_locations_tmp[3 * idx + 1];
      all_locations[idx][2] = all_locations_tmp[3 * idx + 2];
    }

    auto out = std::all_of(all_locations.begin(), all_locations.end(),
                           [&](const auto &el) {
                             for (int i = 0; i < 3; ++i) {
                               if (el[i] <= -1e6)
                                 return false;
                             }
                             return true;
                           });

    if (!out)
      mgis::raise(
          "Error in fetching the positions from the distributed processes.");

    return all_locations;
  } else {
    return local_locations;
  }
} // end of getPointsAboveStressThresholdImplementation

template <bool parallel>
static std::vector<std::pair<mfem_mgis::real, std::array<mfem_mgis::real, 3u>>>
getPointsandStressAboveStressThresholdImplementation(
    const mfem_mgis::Material &m, const mfem_mgis::real threshold) {
  constexpr mfem_mgis::size_type stress_size = 6;
  const auto &s = m.getPartialQuadratureSpace();
  const auto &fed = s.getFiniteElementDiscretization();
  const auto &fespace = fed.getFiniteElementSpace<parallel>();
  // material identifier
  const auto mid = s.getId();
  // values of the stress at the end of the time step
  const auto *stress_values = m.s1.thermodynamic_forces.data();
  // std::vector<std::array<mfem_mgis::real, 3u>> local_locations;
  std::vector<std::pair<mfem_mgis::real, std::array<mfem_mgis::real, 3u>>>
      collected_points_with_svp;
  // loop over the elements
  mfem::Vector tmp;
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
      if (mvp > threshold) {
        const auto &ip = ir.IntPoint(g);
        tr.SetIntPoint(&ip);
        tr.Transform(tr.GetIntPoint(), tmp);
        collected_points_with_svp.push_back({mvp, {tmp[0], tmp[1], tmp[2]}});
      }
    }
  }
  std::vector<std::pair<mfem_mgis::real, std::array<mfem_mgis::real, 3u>>>
      locations_and_svp_values;

  if constexpr (parallel) {
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    const auto lsize = static_cast<int>(collected_points_with_svp.size());
    const int values_per_item = 4;
    std::vector<mfem_mgis::real> local_items_tmp;
    local_items_tmp.resize(lsize * values_per_item);
    for (int idx = 0; idx != lsize; ++idx) {
      local_items_tmp[values_per_item * idx] =
          collected_points_with_svp[idx].first;
      local_items_tmp[values_per_item * idx + 1] =
          collected_points_with_svp[idx].second[0];
      local_items_tmp[values_per_item * idx + 2] =
          collected_points_with_svp[idx].second[1];
      local_items_tmp[values_per_item * idx + 3] =
          collected_points_with_svp[idx].second[2];
      // std::cout << local_locations[idx][0] << "\t" << local_locations[idx][1]
      // << "\t" << local_locations[idx][2] << std::endl;
    }

    std::vector<int> item_counts_per_proc(nprocs);
    MPI_Allgather(&lsize, 1, MPI_INT, item_counts_per_proc.data(), 1, MPI_INT,
                  MPI_COMM_WORLD);

    int gsize_items = 0;
    std::vector<int> displacements_items(nprocs);

    if (nprocs > 0) {
      displacements_items[0] = 0;
      gsize_items += item_counts_per_proc[0];
      for (int i = 1; i < nprocs; ++i) {
        displacements_items[i] =
            displacements_items[i - 1] + item_counts_per_proc[i - 1];
        gsize_items += item_counts_per_proc[i];
      }
    }

    std::vector<int> counts_recv_doubles(nprocs);
    std::vector<int> displacements_doubles(nprocs);
    for (int i = 0; i < nprocs; ++i) {
      counts_recv_doubles[i] = item_counts_per_proc[i] * values_per_item;
      displacements_doubles[i] = displacements_items[i] * values_per_item;
    }

    std::vector<mfem_mgis::real> all_items_flat(gsize_items * values_per_item);
    if (gsize_items > 0 ||
        lsize > 0) { // Only call MPI if there's data to send/recv
      MPI_Allgatherv(local_items_tmp.data(), lsize * values_per_item,
                     MPI_DOUBLE, all_items_flat.data(),
                     counts_recv_doubles.data(), displacements_doubles.data(),
                     MPI_DOUBLE, MPI_COMM_WORLD);
    }

    locations_and_svp_values.resize(gsize_items);
    for (int idx = 0; idx != gsize_items; ++idx) {
      locations_and_svp_values[idx].first =
          all_items_flat[idx * values_per_item + 0];
      locations_and_svp_values[idx].second[0] =
          all_items_flat[idx * values_per_item + 1]; // xyz
      locations_and_svp_values[idx].second[1] =
          all_items_flat[idx * values_per_item + 2];
      locations_and_svp_values[idx].second[2] =
          all_items_flat[idx * values_per_item + 3];
    }
    return locations_and_svp_values;
  } else {
    return collected_points_with_svp;
  }
} // end of getPointsandStressAboveStressThresholdImplementation

FirstPrincipalStressValueAndLocation
findFirstPrincipalStressValueAndLocation(const mfem_mgis::Material &m) {
  const auto &s = m.getPartialQuadratureSpace();
  const auto &fed = s.getFiniteElementDiscretization();
  if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
    return findFirstPrincipalStressValueAndLocationImplementation<true>(m);
#else  /* MFEM_USE_MPI */
    mfem_mgis::reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  }
  return findFirstPrincipalStressValueAndLocationImplementation<false>(m);
} // end of findFirstPrincipalStressValueAndLocation

std::vector<std::array<mfem_mgis::real, 3u>>
getPointsAboveStressThreshold(mfem_mgis::Material &m, const mfem_mgis::real v) {
  const auto &s = m.getPartialQuadratureSpace();
  const auto &fed = s.getFiniteElementDiscretization();
  if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
    return getPointsAboveStressThresholdImplementation<true>(m, v);
#else  /* MFEM_USE_MPI */
    mfem_mgis::reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  }
  return getPointsAboveStressThresholdImplementation<false>(m, v);
}

std::vector<std::pair<mfem_mgis::real, std::array<mfem_mgis::real, 3u>>>
getPointsandStressAboveStressThreshold(mfem_mgis::Material &m,
                                       const mfem_mgis::real v) {
  const auto &s = m.getPartialQuadratureSpace();
  const auto &fed = s.getFiniteElementDiscretization();
  if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
    return getPointsandStressAboveStressThresholdImplementation<true>(m, v);
#else  /* MFEM_USE_MPI */
    mfem_mgis::reportUnsupportedParallelComputations();
#endif /* MFEM_USE_MPI */
  }
  return getPointsandStressAboveStressThresholdImplementation<false>(m, v);
}

} // end of namespace opera_hpc
