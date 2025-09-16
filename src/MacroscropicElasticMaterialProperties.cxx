/*!
 * \file   MacroscropicElasticMaterialProperties.cxx
 * \brief
 * \author Thomas Helfer
 * \date   16/09/2025
 */

#include <cmath>
#include "MFEMMGIS/Profiler.hxx"
#include "MM_OPERA_HPC/MacroscropicElasticMaterialProperties.hxx"

namespace mm_opera_hpc {

void MacroscropicElasticMaterialProperties::update(
    const mfem_mgis::real E, const mfem_mgis::real nu,
    const mfem_mgis::real G) noexcept {
  // ===========================
  // effective moduli estimations
  // ===========================
  const mfem_mgis::real C11 = E * (nu - 1) / (2 * nu * nu + nu - 1);
  const mfem_mgis::real C12 = -(nu * E) / (2 * nu * nu + nu - 1);
  const mfem_mgis::real C44 = G;
  // ---- Bulk Modulus ----
  const mfem_mgis::real Keff = 1. / 3. * (C11 + 2 * C12);
  // ---- Shear Modulus ----
  // Reuss & Voigt bounds
  const mfem_mgis::real aniso = 2 * C44 / (C11 - C12);
  const mfem_mgis::real GReuss = 5 * C44 / (3 + 2 * aniso);
  const mfem_mgis::real GVoigt = (3 * aniso + 2) / (5 * aniso) * C44;
  // Hashin–Shtrikman bounds
  const mfem_mgis::real Gv = (C11 - C12) / 2;
  const mfem_mgis::real beta1 =
      -3 * (Keff + 2 * Gv) / (5 * Gv * (3 * Keff + 4 * Gv));
  const mfem_mgis::real beta2 =
      -3 * (Keff + 2 * C44) / (5 * C44 * (3 * Keff + 4 * C44));
  const mfem_mgis::real GHSl = Gv + 3 * ((C44 - Gv) - 4 * beta1) / 5;
  const mfem_mgis::real GHSu = C44 + 3 * ((Gv - C44) - 6 * beta2) / 5;
  // Kroener
  mfem_mgis::real m0, mus;
  mfem_mgis::real D1, sD1, D2;
  mfem_mgis::real GK = 0.;
  do {
    m0 = GK;
    D1 = Keff + 2 * GK;
    sD1 = 6 * D1;
    mus = GK * (9 * Keff + 8 * GK) / sD1;
    D2 = 5 * mus + 2 * C44 + 3 * Gv;
    GK = (3 * C44 * mus + 2 * Gv * mus + 5 * Gv * C44) / D2;
  } while (std::abs(GK - m0) / Keff > 1e-10);

  // Shear modulus choice
  mfem_mgis::Profiler::Utils::Message("GReuss =", GReuss);
  mfem_mgis::Profiler::Utils::Message("GVoigt =", GVoigt);
  mfem_mgis::Profiler::Utils::Message("GHSl =", GHSl);
  mfem_mgis::Profiler::Utils::Message("GHSu =", GHSu);
  mfem_mgis::Profiler::Utils::Message("GK =", GK);

  this->Geff = GK; // GHSu
                   // Young and Poisson effective moduli
  this->nueff = (3 * Keff - 2 * Geff) / (2 * (3 * Keff + Geff));
  this->Eeff = 9 * Keff * Geff / (3 * Keff + Geff);
} // end of update

} // end of namespace mm_opera_hpc