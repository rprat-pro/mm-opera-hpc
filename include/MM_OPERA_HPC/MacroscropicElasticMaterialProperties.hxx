/*!
 * \file   MacroscropicElasticMaterialProperties.hxx
 * \brief
 * \author Thomas Helfer
 * \date   16/09/2025
 */

#ifndef LIB_MM_OPERA_HPC_MACROSCROPICELASTICMATERIALPROPERTIES_HXX
#define LIB_MM_OPERA_HPC_MACROSCROPICELASTICMATERIALPROPERTIES_HXX

#include "MM_OPERA_HPC/Config.hxx"

namespace mm_opera_hpc {

  struct MM_OPERA_HPC_EXPORT MacroscropicElasticMaterialProperties {
    mfem_mgis::real nueff = mfem_mgis::real{};
    mfem_mgis::real Eeff = mfem_mgis::real{};
    mfem_mgis::real Geff = mfem_mgis::real{};

    /*!
     * \brief compute the macroscopic elastic material properties
     *
     * \param[in] young: Young's modulus
     * \param[in] nu: Poisson's ratio
     * \param[in] G: shear modulus
     */
    void update(const mfem_mgis::real,
                const mfem_mgis::real,
                const mfem_mgis::real) noexcept;
  };  // struct MacroscropicElasticMaterialProperties

}  // end of namespace mm_opera_hpc

#endif /* LIB_MM_OPERA_HPC_MACROSCROPICELASTICMATERIALPROPERTIES_HXX */
