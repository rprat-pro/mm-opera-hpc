/*!
 * \file   src/Utilities.cxx
 * \brief  Utility functions for stress analysis in finite element simulations
 * \author Thomas Helfer
 * \date   29/09/2024
 *
 * This file provides utilities to extract and analyze stress fields from
 * finite element simulations.
 */

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

  /**
   * \brief Find the maximum first principal stress and its location in the
   * domain
   *
   * This function computes the eigenvalues of the stress tensor at every
   * integration point in the mesh, identifies the maximum first principal
   * stress, and returns both its value and spatial location.
   *
   * In parallel mode, results from all MPI processes are gathered and the
   * global maximum is determined.
   *
   * \tparam parallel Whether running in parallel (MPI) mode
   * \param m Material containing stress field data
   * \return Structure containing the maximum principal stress value and its 3D
   * location
   */
  template <bool parallel>
  static FirstPrincipalStressValueAndLocation
  findFirstPrincipalStressValueAndLocationImplementation(
      const mfem_mgis::Material &m) {
    CatchTimeSection("OperaHPC::FindFirstPrincipal");

    // Stress tensor has 6 components
    constexpr mfem_mgis::size_type stress_size = 6;

    // Get finite element space information
    const auto &s = m.getPartialQuadratureSpace();
    const auto &fed = s.getFiniteElementDiscretization();
    const auto &fespace = fed.getFiniteElementSpace<parallel>();

    // Material identifier for filtering elements
    const auto mid = s.getId();

    // Pointer to stress values at end of time step
    const auto *stress_values = m.s1.thermodynamic_forces.data();

    // Initialize with very small value to find maximum
    auto max_stress = -std::numeric_limits<mfem_mgis::real>::max();
    std::array<mfem_mgis::real, 3u> max_stress_location;
    mfem::Vector tmp;

    // Loop over all elements in the mesh
    {
      CatchTimeSection("FindFirstPrincipal::LoopOverElem");
      for (mfem_mgis::size_type i = 0; i != fespace.GetNE(); ++i) {
        // Skip elements not belonging to this material
        if (fespace.GetAttribute(i) != mid) {
          continue;
        }

        const auto &fe = *(fespace.GetFE(i));
        auto &tr = *(fespace.GetElementTransformation(i));
        const auto &ir = s.getIntegrationRule(fe, tr);

        // Offset in the global stress array for this element
        const auto eo = s.getOffset(i);

        // Loop over integration (quadrature) points within the element
        for (mfem_mgis::size_type g = 0; g != ir.GetNPoints(); ++g) {
          // Get stress tensor at this integration point
          const auto *const stress = stress_values + (eo + g) * stress_size;

          // Create stress tensor object in the tensorial algebra library
          // tfel::math
          auto sig = tfel::math::stensor<3u, mfem_mgis::real>{stress};

          // Compute eigenvalues
          const auto sig_vp = sig.computeEigenValues<
              tfel::math::stensor_common::FSESJACOBIEIGENSOLVER>();

          // Find maximum eigenvalue (first principal stress)
          const auto mvp = *(tfel::fsalgo::max_element<3>::exe(sig_vp.begin()));

          // Update global maximum if this value is larger
          if (mvp > max_stress) {
            max_stress = mvp;

            // Get spatial location of this integration point
            const auto &ip = ir.IntPoint(g);
            tr.SetIntPoint(&ip);
            tr.Transform(tr.GetIntPoint(), tmp);
            max_stress_location = {tmp[0], tmp[1], tmp[2]};
          }
        }
      }
    }

    // In parallel mode, gather results from all processes
    if constexpr (parallel) {
      CatchTimeSection("FindFirstPrincipal::GatherResults");

      int nprocs;
      MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

      // Collect all maximum stresses and locations from all processes
      auto all_max_stresses = std::vector<double>(nprocs);
      auto all_locations = std::vector<double>(3 * nprocs);

      MPI_Allgather(&max_stress, 1, MPI_DOUBLE, all_max_stresses.data(), 1,
                    MPI_DOUBLE, MPI_COMM_WORLD);
      MPI_Allgather(max_stress_location.data(), 3, MPI_DOUBLE,
                    all_locations.data(), 3, MPI_DOUBLE, MPI_COMM_WORLD);

      // Find the process with the global maximum
      const auto i =
          std::max_element(all_max_stresses.begin(), all_max_stresses.end()) -
          all_max_stresses.begin();

      // Extract the global maximum stress and its location
      max_stress = all_max_stresses[i];
      std::copy(all_locations.begin() + 3 * i,
                all_locations.begin() + 3 * (i + 1),
                max_stress_location.begin());
    }

    return {.value = max_stress, .location = max_stress_location};
  }  // end of findFirstPrincipalStressValueAndLocationImplementation

  /**
   * \brief Get all spatial locations where principal stress exceeds a threshold
   *
   * Scans all integration points in the mesh and collects the 3D coordinates
   * of points where the first principal stress exceeds the specified threshold.
   * This is useful for identifying regions of high stress concentration.
   *
   * In parallel mode, locations from all processes are gathered into a single
   * list.
   *
   * \tparam parallel Whether running in parallel (MPI) mode
   * \param m Material containing stress field data
   * \param threshold Stress value threshold for filtering points
   * \return Vector of 3D coordinates where stress exceeds threshold
   */
  template <bool parallel>
  static std::vector<std::array<mfem_mgis::real, 3u>>
  getPointsAboveStressThresholdImplementation(const mfem_mgis::Material &m,
                                              const mfem_mgis::real threshold) {
    CatchTimeSection("OperaHPC::GetPointsAbove");

    constexpr mfem_mgis::size_type stress_size = 6;
    const auto &s = m.getPartialQuadratureSpace();
    const auto &fed = s.getFiniteElementDiscretization();
    const auto &fespace = fed.getFiniteElementSpace<parallel>();

    const auto mid = s.getId();
    const auto *stress_values = m.s1.thermodynamic_forces.data();

    // Local collection of locations exceeding threshold
    std::vector<std::array<mfem_mgis::real, 3u>> local_locations;
    mfem::Vector tmp;

    // Loop over elements and integration points
    {
      CatchTimeSection("GetPointsAbove::LoopOverElem");
      for (mfem_mgis::size_type i = 0; i != fespace.GetNE(); ++i) {
        if (fespace.GetAttribute(i) != mid) {
          continue;
        }

        const auto &fe = *(fespace.GetFE(i));
        auto &tr = *(fespace.GetElementTransformation(i));
        const auto &ir = s.getIntegrationRule(fe, tr);
        const auto eo = s.getOffset(i);

        for (mfem_mgis::size_type g = 0; g != ir.GetNPoints(); ++g) {
          const auto *const stress = stress_values + (eo + g) * stress_size;
          auto sig = tfel::math::stensor<3u, mfem_mgis::real>{stress};

          // Compute principal stresses
          const auto sig_vp = sig.computeEigenValues<
              tfel::math::stensor_common::FSESJACOBIEIGENSOLVER>();
          const auto mvp = *(tfel::fsalgo::max_element<3>::exe(sig_vp.begin()));

          // If stress exceeds threshold, store location
          if (mvp > threshold) {
            const auto &ip = ir.IntPoint(g);
            tr.SetIntPoint(&ip);
            tr.Transform(tr.GetIntPoint(), tmp);
            local_locations.push_back({tmp[0], tmp[1], tmp[2]});
          }
        }
      }
    }

    // Gather results from all MPI processes
    if constexpr (parallel) {
      CatchTimeSection("GetPointsAbove::GatherResults");

      int nprocs;
      MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

      // Flatten local locations into a single array (x1,y1,z1,x2,y2,z2,...)
      const auto lsize = static_cast<int>(local_locations.size());
      std::vector<double> local_locations_tmp;
      local_locations_tmp.resize(3 * lsize);
      for (int idx = 0; idx != lsize; ++idx) {
        local_locations_tmp[3 * idx] = local_locations[idx][0];
        local_locations_tmp[3 * idx + 1] = local_locations[idx][1];
        local_locations_tmp[3 * idx + 2] = local_locations[idx][2];
      }

      // Calculate total number of locations across all processes
      int gsize;
      MPI_Allreduce(&lsize, &gsize, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

      // Prepare arrays for variable-length gather operation
      std::vector<double> all_locations_tmp(3 * gsize);
      std::vector<int> counts_recv(nprocs);
      std::vector<int> displacements(nprocs);

      // Get count from each process
      MPI_Allgather(&lsize, 1, MPI_INT, counts_recv.data(), 1, MPI_INT,
                    MPI_COMM_WORLD);

      // Calculate displacement offsets for each process
      displacements[0] = 0;
      for (int i = 1; i < nprocs; ++i) {
        displacements[i] = displacements[i - 1] + counts_recv[i - 1] * 3;
      }

      // Multiply counts by 3 (for x, y, z coordinates)
      for (int i = 0; i < nprocs; ++i) {
        counts_recv[i] *= 3;
      }

      // Gather all locations from all processes
      MPI_Allgatherv(local_locations_tmp.data(), 3 * lsize, MPI_DOUBLE,
                     all_locations_tmp.data(), counts_recv.data(),
                     displacements.data(), MPI_DOUBLE, MPI_COMM_WORLD);

      // Reconstruct array of 3D locations
      std::vector<std::array<mfem_mgis::real, 3u>> all_locations;
      all_locations.resize(gsize);

      // Initialize with sentinel values for error checking
      for (auto &locs : all_locations) {
        std::fill(locs.begin(), locs.end(), -1e6);
      }

      // Unflatten the gathered data
      for (int idx = 0; idx != gsize; ++idx) {
        all_locations[idx][0] = all_locations_tmp[3 * idx];
        all_locations[idx][1] = all_locations_tmp[3 * idx + 1];
        all_locations[idx][2] = all_locations_tmp[3 * idx + 2];
      }

      // Verify all locations were properly gathered (no sentinel values remain)
      auto out = std::all_of(all_locations.begin(), all_locations.end(),
                             [&](const auto &el) {
                               for (int i = 0; i < 3; ++i) {
                                 if (el[i] <= -1e6) return false;
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
  }  // end of getPointsAboveStressThresholdImplementation

  /**
   * \brief Get locations AND stress values where principal stress exceeds a
   * threshold
   *
   * Similar to getPointsAboveStressThreshold, but also returns the actual
   * stress value at each location.
   *
   * \tparam parallel Whether running in parallel (MPI) mode
   * \param m Material containing stress field data
   * \param threshold Stress value threshold for filtering points
   * \return Vector of structures containing both stress values and 3D locations
   */
  template <bool parallel>
  static std::vector<StressValueAndLocation>
  getPointsandStressAboveStressThresholdImplementation(
      const mfem_mgis::Material &m, const mfem_mgis::real threshold) {
    constexpr mfem_mgis::size_type stress_size = 6;
    const auto &s = m.getPartialQuadratureSpace();
    const auto &fed = s.getFiniteElementDiscretization();
    const auto &fespace = fed.getFiniteElementSpace<parallel>();

    const auto mid = s.getId();
    const auto *stress_values = m.s1.thermodynamic_forces.data();

    // Collect both stress values and locations
    std::vector<StressValueAndLocation> collected_points_with_svp;
    mfem::Vector tmp;

    // Loop over elements and integration points
    for (mfem_mgis::size_type i = 0; i != fespace.GetNE(); ++i) {
      if (fespace.GetAttribute(i) != mid) {
        continue;
      }

      const auto &fe = *(fespace.GetFE(i));
      auto &tr = *(fespace.GetElementTransformation(i));
      const auto &ir = s.getIntegrationRule(fe, tr);
      const auto eo = s.getOffset(i);

      for (mfem_mgis::size_type g = 0; g != ir.GetNPoints(); ++g) {
        const auto *const stress = stress_values + (eo + g) * stress_size;
        auto sig = tfel::math::stensor<3u, mfem_mgis::real>{stress};

        // Compute principal stresses
        const auto sig_vp = sig.computeEigenValues<
            tfel::math::stensor_common::FSESJACOBIEIGENSOLVER>();
        const auto mvp = *(tfel::fsalgo::max_element<3>::exe(sig_vp.begin()));

        // Store both stress value and location if above threshold
        if (mvp > threshold) {
          const auto &ip = ir.IntPoint(g);
          tr.SetIntPoint(&ip);
          tr.Transform(tr.GetIntPoint(), tmp);
          collected_points_with_svp.push_back({mvp, {tmp[0], tmp[1], tmp[2]}});
        }
      }
    }

    std::vector<StressValueAndLocation> locations_and_svp_values;

    // Gather results in parallel mode
    if constexpr (parallel) {
      int nprocs;
      MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

      const auto lsize = static_cast<int>(collected_points_with_svp.size());

      // Each item has 4 values: stress_value, x, y, z
      const int values_per_item = 4;
      std::vector<mfem_mgis::real> local_items_tmp;
      local_items_tmp.resize(lsize * values_per_item);

      // Flatten data: [stress1, x1, y1, z1, stress2, x2, y2, z2, ...]
      for (int idx = 0; idx != lsize; ++idx) {
        local_items_tmp[values_per_item * idx] =
            collected_points_with_svp[idx].value;
        local_items_tmp[values_per_item * idx + 1] =
            collected_points_with_svp[idx].location[0];
        local_items_tmp[values_per_item * idx + 2] =
            collected_points_with_svp[idx].location[1];
        local_items_tmp[values_per_item * idx + 3] =
            collected_points_with_svp[idx].location[2];
      }

      // Get item counts from all processes
      std::vector<int> item_counts_per_proc(nprocs);
      MPI_Allgather(&lsize, 1, MPI_INT, item_counts_per_proc.data(), 1, MPI_INT,
                    MPI_COMM_WORLD);

      // Calculate total number of items and displacements
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

      // Convert item counts/displacements to double counts/displacements
      std::vector<int> counts_recv_doubles(nprocs);
      std::vector<int> displacements_doubles(nprocs);
      for (int i = 0; i < nprocs; ++i) {
        counts_recv_doubles[i] = item_counts_per_proc[i] * values_per_item;
        displacements_doubles[i] = displacements_items[i] * values_per_item;
      }

      // Gather all flattened data from all processes
      std::vector<mfem_mgis::real> all_items_flat(gsize_items *
                                                  values_per_item);
      if (gsize_items > 0 || lsize > 0) {
        MPI_Allgatherv(local_items_tmp.data(), lsize * values_per_item,
                       MPI_DOUBLE, all_items_flat.data(),
                       counts_recv_doubles.data(), displacements_doubles.data(),
                       MPI_DOUBLE, MPI_COMM_WORLD);
      }

      // Reconstruct structured data from flattened array
      locations_and_svp_values.resize(gsize_items);
      for (int idx = 0; idx != gsize_items; ++idx) {
        locations_and_svp_values[idx].value =
            all_items_flat[idx * values_per_item + 0];
        locations_and_svp_values[idx].location = {
            all_items_flat[idx * values_per_item + 1],   // x
            all_items_flat[idx * values_per_item + 2],   // y
            all_items_flat[idx * values_per_item + 3]};  // z
      }
      return locations_and_svp_values;
    } else {
      return collected_points_with_svp;
    }
  }  // end of getPointsandStressAboveStressThresholdImplementation

  /**
   * \brief Public interface to find maximum principal stress and location
   *
   * Detects whether computation is parallel or serial
   * and calls the appropriate implementation.
   *
   * \param m Material containing stress field data
   * \return Maximum principal stress value and its 3D location
   */
  FirstPrincipalStressValueAndLocation findFirstPrincipalStressValueAndLocation(
      const mfem_mgis::Material &m) {
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
  }  // end of findFirstPrincipalStressValueAndLocation

  /**
   * \brief Public interface to get locations above stress threshold
   *
   * \param m Material containing stress field data
   * \param v Stress threshold value
   * \return Vector of 3D locations where stress exceeds threshold
   */
  std::vector<std::array<mfem_mgis::real, 3u>> getPointsAboveStressThreshold(
      const mfem_mgis::Material &m, const mfem_mgis::real v) {
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

  /**
   * \brief Public interface to get locations and stress values above threshold
   *
   * \param m Material containing stress field data
   * \param v Stress threshold value
   * \return Vector of structures containing stress values and 3D locations
   */
  std::vector<StressValueAndLocation> getPointsandStressAboveStressThreshold(
      const mfem_mgis::Material &m, const mfem_mgis::real v) {
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

}  // end of namespace opera_hpc