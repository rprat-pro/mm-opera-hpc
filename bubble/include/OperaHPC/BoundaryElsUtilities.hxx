

#ifndef LIB_OPERA_HPC_BOUNDARYELSUTILITIES_HXX
#define LIB_OPERA_HPC_BOUNDARYELSUTILITIES_HXX

#include "MFEMMGIS/Config.hxx"
#include <MFEMMGIS/Material.hxx>
#include <array>
#include <vector>

namespace opera_hpc {

    struct StressValueAndId{
     mfem_mgis::real value;
     mfem_mgis::size_type boundary_id;
    //  std::array<mfem_mgis::real, 3u> position;

    };

    std::vector<StressValueAndId> accessBoundaryClosestQPData(
        const mfem_mgis::Material &m, std::vector<mfem_mgis::size_type> &target_boundary_attribute
    );

}

#endif