

#include "TFEL/Math/power.hxx"
#include "mfem/fem/fespace.hpp"
#include "mfem/linalg/densemat.hpp"
#include <MFEMMGIS/Config.hxx>
#include <TFEL/Math/stensor.hxx>
#include <limits>
#include <vector>
#ifdef MFEM_USE_MPI
#include "mfem/fem/pfespace.hpp"
#endif /* MFEM_USE_MPI */
#include "MFEMMGIS/FiniteElementDiscretization.hxx"
#include "MFEMMGIS/Material.hxx"
#include "MFEMMGIS/PartialQuadratureSpace.hxx"
#include "MFEMMGIS/Profiler.hxx"
#include "OperaHPC/BoundaryElsUtilities.hxx"

namespace opera_hpc {

  template <bool parallel>
  std::vector<StressValueAndId> accessBoundaryClosestQPDataImplementation(
      const mfem_mgis::Material &m,
      std::vector<mfem_mgis::size_type> &target_boundary_attribute) {
    constexpr mfem_mgis::size_type stress_size = 6;

    const auto &pqs = m.getPartialQuadratureSpace();
    const auto &fed = pqs.getFiniteElementDiscretization();
    const auto &fec = fed.getFiniteElementSpace<parallel>();

    // Get the underlying MFEM mesh object from the FES.
    auto *mesh = fec.GetMesh();

    const auto &fespace_mgis = fed.getFiniteElementSpace<parallel>();
    const auto mid = pqs.getId();
    const auto *stress_values_ptr = m.s1.thermodynamic_forces.data();

    // std::cout << "Accessing data for boundary attribute: " << std::endl;
    // for (auto &bid : target_boundary_attribute)
    //   std::cout << bid << std::endl;

    // std::cout << "A total of " << target_boundary_attribute.size()
    //           << " surfaces has been requested.\n";

    auto storing_results = std::vector<StressValueAndId>{};

    // Iterate over each target boundary attribute.
    for (mfem_mgis::size_type &bid : target_boundary_attribute) {
      // Iterate over all boundary elements in the mesh.
      for (int i = 0; i < mesh->GetNBE(); i++) {
        if (mesh->GetBdrAttribute(i) == bid) {
          // Get the element transformation for the current boundary element
          // This transformation maps from reference space to physical space for
          // the boundary face
          mfem::ElementTransformation *bdr_tr =
              mesh->GetBdrElementTransformation(i);

          // Get face-element transformations.
          // This provides information about the domain element adjacent to this
          // boundary face. For boundary faces, Elem1 is the interior element
          mfem::FaceElementTransformations *face_elem_tr =
              mesh->GetFaceElementTransformations(
                  i, true);  // NB true selects the boundary elements

          // If face_elem_tr is null, something is wrong
          if (!face_elem_tr) continue;

          int domain_elem_idx = face_elem_tr->Elem1No;

          // Check if the parent domain element belongs to the current material.
          // If not, its quadrature points are not managed by this 'pqs' and
          // 'm'.
          if (fespace_mgis.GetAttribute(domain_elem_idx) != mid) {
            continue;
          }

          // Get information for the parent domain element to build an
          // integration rule we need the finite element type and the polynomial
          // degree
          const auto &domain_fe_mgis = *(fespace_mgis.GetFE(domain_elem_idx));
          auto &domain_elem_tr =
              *(fespace_mgis.GetElementTransformation(domain_elem_idx));

          // Get the integration rule (set of quadrature points) for the domain
          // element as defined by the partial quadrature space.
          const auto &domain_qps_ir_mgis =
              pqs.getIntegrationRule(domain_fe_mgis, domain_elem_tr);
          // Get the offset for this domain element like thomas
          const auto domain_elem_offset = pqs.getOffset(domain_elem_idx);

          // std::cout << "  Boundary Face " << i << " (Attr "
          //           << mesh->GetBdrAttribute(i) << "), Parent Domain Element
          //           "
          //           << domain_elem_idx << std::endl;

          // Determine the integration rule for the boundary face.

          int bdr_geom = bdr_tr->GetGeometryType();
          // The order of the integration rule for the boundary face:
          // Jacobian order of the transformation and the order of the domain
          // FE.
          int bdr_ir_order = bdr_tr->OrderJ() + domain_fe_mgis.GetOrder();
          const mfem::IntegrationRule &bdr_face_ir =
              mfem::IntRules.Get(bdr_geom, bdr_ir_order);

          // Pre-calculate and store the physical coordinates of all quadrature
          // points within the parent domain element.
          std::vector<mfem::Vector> parent_domain_qp_phys_coords_list(
              domain_qps_ir_mgis.GetNPoints());
          for (int dq_idx = 0; dq_idx < domain_qps_ir_mgis.GetNPoints();
               ++dq_idx) {
            // Get the reference coordinates of the current domain quadrature
            // point.
            const mfem::IntegrationPoint &domain_ip_ref =
                domain_qps_ir_mgis.IntPoint(dq_idx);
            // Allocate space for physical coordinates.
            parent_domain_qp_phys_coords_list[dq_idx].SetSize(
                mesh->SpaceDimension());
            // Set the integration point for the domain element transformation.
            domain_elem_tr.SetIntPoint(&domain_ip_ref);
            // Transform reference coordinates (-1,1) to physical coordinates.
            domain_elem_tr.Transform(domain_ip_ref,
                                     parent_domain_qp_phys_coords_list[dq_idx]);
          }

          for (int bg = 0; bg < bdr_face_ir.GetNPoints(); ++bg) {
            // Get the reference coordinates of the current boundary face
            // quadrature point.
            const mfem::IntegrationPoint &bdr_face_ip_ref =
                bdr_face_ir.IntPoint(bg);

            // Transform the boundary face quadrature point to physical
            // coordinates.
            mfem::Vector bdr_face_qp_phys_coords(mesh->SpaceDimension());
            bdr_tr->SetIntPoint(&bdr_face_ip_ref);
            bdr_tr->Transform(bdr_face_ip_ref, bdr_face_qp_phys_coords);

            // Find the closest domain quadrature point (from the parent
            // element) to the current boundary face quadrature point.
            mfem_mgis::real min_dist_sq =
                std::numeric_limits<mfem_mgis::real>::max();
            int closest_domain_qp_g_idx = -1;

            for (auto dq_idx = 0; dq_idx < domain_qps_ir_mgis.GetNPoints();
                 ++dq_idx) {
              const auto &current_domain_qp_phys =
                  parent_domain_qp_phys_coords_list[dq_idx];
              mfem_mgis::real dist_sq = 0;
              // Calculate distance
              for (int d = 0; d < current_domain_qp_phys.Size(); ++d) {
                mfem_mgis::real diff =
                    current_domain_qp_phys(d) - bdr_face_qp_phys_coords(d);
                dist_sq += tfel::math::power<2>(diff);
              }

              // If this domain QP is closer than the previous closest, update
              if (dist_sq < min_dist_sq) {
                min_dist_sq = dist_sq;
                closest_domain_qp_g_idx = dq_idx;
              }
            }

            if (closest_domain_qp_g_idx != -1) {
              // The offset is (domain_elem_offset + local_qp_index) *
              // number_of_stress_components.
              const auto *const stress_at_closest_domain_qp =
                  stress_values_ptr +
                  (domain_elem_offset + closest_domain_qp_g_idx) * stress_size;

              auto sig = tfel::math::stensor<3u, mfem_mgis::real>{
                  stress_at_closest_domain_qp};

              // Compute the eigenvalues of the stress tensor -> principal
              // stress
              const auto sig_vp = sig.computeEigenValues<
                  tfel::math::stensor_common::FSESJACOBIEIGENSOLVER>();

              // Find the largest eigenvalue
              const auto mvp =
                  *(tfel::fsalgo::max_element<3>::exe(sig_vp.begin()));

              auto boundary_info = opera_hpc::StressValueAndId{};
              boundary_info.value = mvp;
              boundary_info.boundary_id = bid;
              storing_results.push_back(boundary_info);
            } else {
              // This case should ideally not happen
              std::cout
                  << "      Warning: Could not find a closest domain QP\n";
            }
          }
        }
      }
    }
    return storing_results;
  }

  std::vector<StressValueAndId> accessBoundaryClosestQPData(
      const mfem_mgis::Material &m,
      std::vector<mfem_mgis::size_type> &target_boundary_attribute) {
    const auto &s = m.getPartialQuadratureSpace();
    const auto &fed = s.getFiniteElementDiscretization();

    // Check if the finite element discretization describes a parallel
    // computation.
    if (fed.describesAParallelComputation()) {
#ifdef MFEM_USE_MPI
      return accessBoundaryClosestQPDataImplementation<true>(
          m, target_boundary_attribute);
#else  /* MFEM_USE_MPI */
      mfem_mgis::reportUnsupportedParallelComputations();
      return {};
#endif /* MFEM_USE_MPI */
    }
    return accessBoundaryClosestQPDataImplementation<false>(
        m, target_boundary_attribute);
  }

}  // namespace opera_hpc