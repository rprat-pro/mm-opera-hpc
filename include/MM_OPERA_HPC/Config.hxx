/*!
 * \file   Config.hxx
 * \brief    
 * \author Thomas Helfer
 * \date   16/09/2025
 */

#ifndef LIB_MM_OPERA_HPC_CONFIG_HXX
#define LIB_MM_OPERA_HPC_CONFIG_HXX

#include "MFEMMGIS/Config.hxx"

#define MM_OPERA_HPC_VISIBILITY_LOCAL MGIS_VISIBILITY_LOCAL

#if defined _WIN32 || defined _WIN64 || defined __CYGWIN__
#if defined MMOPERAHPC_EXPORTS
#define MM_OPERA_HPC_EXPORT MGIS_VISIBILITY_EXPORT
#else /* defined MMOPERAHPC_EXPORTS */
#ifndef MM_OPERA_HPC_STATIC_BUILD
#define MM_OPERA_HPC_EXPORT MGIS_VISIBILITY_IMPORT
#else /* MM_OPERA_HPC_STATIC_BUILD */
#define MM_OPERA_HPC_EXPORT
#endif /* MM_OPERA_HPC_STATIC_BUILD */
#endif /* defined MMOPERAHPC_EXPORTS */
#else  /* defined _WIN32 || defined _WIN64 || defined __CYGWIN__ */
#define MM_OPERA_HPC_EXPORT MGIS_VISIBILITY_EXPORT
#endif /* */


#endif /* LIB_MM_OPERA_HPC_CONFIG_HXX */
