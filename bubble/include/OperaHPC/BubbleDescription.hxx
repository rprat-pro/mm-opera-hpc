/*!
 * \file   BubbleDescription.hxx
 * \brief
 * \author Thomas Helfer
 * \date   29/09/2024
 */

#ifndef LIB_OPERA_HPC_BUBBLEDESCRIPTION_HXX
#define LIB_OPERA_HPC_BUBBLEDESCRIPTION_HXX

#include <array>
#include <vector>
#include <string>
#include "MFEMMGIS/Config.hxx"

namespace opera_hpc {

  /*!
   * \brief a data structure describing a bubble
   */
  struct BubbleDescription {
    /*!
     * \return the list of bubble descriptions
     * \param[in] f: file name
     */
    static std::vector<BubbleDescription> read(const std::string&);
    //! \brief boundary identifier
    mfem_mgis::size_type boundary_identifier;
    //! \brief radius of the bubble
    mfem_mgis::real radius;
    //! \brief center of the bubble
    std::array<mfem_mgis::real, 3u> center;
  };

  /*!
   * \return the distance between the given bubble and the given point
   * \param[in] b: bubble
   * \param[in] p: point
   */
  mfem_mgis::real distance(const BubbleDescription&,
                           const std::array<mfem_mgis::real, 3u>&) noexcept;

}  // end of namespace opera_hpc

#endif /* LIB_OPERA_HPC_BUBBLEDESCRIPTION_HXX */
