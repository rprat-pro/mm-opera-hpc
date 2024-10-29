/*!
* \file   umatMonoCristal_UO2.hxx
* \brief  This file declares the umat interface for the MonoCristal_UO2 behaviour law
* \author Luc Portelette / Thomas Helfer / Etienne Castelier
* \date   
*/

#ifndef LIB_CASTEM_MONOCRISTAL_UO2_HXX
#define LIB_CASTEM_MONOCRISTAL_UO2_HXX

#include"castem.h"
#ifdef umat
#undef umat
#endif /* umat */

#include"TFEL/Config/TFELConfig.hxx"

#include"MFront/Castem/Castem.hxx"

#ifdef __cplusplus
#include"MFront/Castem/CastemTraits.hxx"
#include"MFront/Castem/CastemOrthotropicBehaviour.hxx"
#include"TFEL/Material/MonoCristal_UO2.hxx"
#endif /* __cplusplus */

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif /* NOMINMAX */
#include <windows.h>
#ifdef small
#undef small
#endif /* small */
#endif /* _WIN32 */

#ifndef MFRONT_SHAREDOBJ
#define MFRONT_SHAREDOBJ TFEL_VISIBILITY_EXPORT
#endif /* MFRONT_SHAREDOBJ */

#ifndef MFRONT_EXPORT_SYMBOL
#define MFRONT_EXPORT_SYMBOL(TYPE, NAME, VALUE) \
  MFRONT_SHAREDOBJ extern TYPE NAME;            \
  MFRONT_SHAREDOBJ TYPE NAME = VALUE
#endif /* MFRONT_EXPORT_SYMBOL*/

#ifndef MFRONT_EXPORT_ARRAY_ARGUMENTS
#define MFRONT_EXPORT_ARRAY_ARGUMENTS(...) __VA_ARGS__
#endif /* MFRONT_EXPORT_ARRAY_ARGUMENTS */

#ifndef MFRONT_EXPORT_ARRAY_OF_SYMBOLS
#define MFRONT_EXPORT_ARRAY_OF_SYMBOLS(TYPE, NAME, SIZE, VALUE) \
  MFRONT_SHAREDOBJ extern TYPE NAME[SIZE];                      \
  MFRONT_SHAREDOBJ TYPE NAME[SIZE] = {VALUE}
#endif /* MFRONT_EXPORT_ARRAY_OF_SYMBOLS*/

#ifdef __cplusplus

namespace castem{

template<typename NumericType>
struct CastemTraits<tfel::material::MonoCristal_UO2<tfel::material::ModellingHypothesis::TRIDIMENSIONAL, NumericType, false> >{
using  ModellingHypothesis = tfel::material::ModellingHypothesis;
using  ModellingHypothesisToSpaceDimension = tfel::material::ModellingHypothesisToSpaceDimension<ModellingHypothesis::TRIDIMENSIONAL>;
static constexpr ModellingHypothesis::Hypothesis H = ModellingHypothesis::TRIDIMENSIONAL;
static constexpr CastemBehaviourType btype  = STANDARDFINITESTRAINBEHAVIOUR;
// space dimension
static constexpr unsigned short N           = ModellingHypothesisToSpaceDimension::value;
// tiny vector size
static constexpr unsigned short TVectorSize = N;
// symmetric tensor size
static constexpr unsigned short StensorSize = tfel::math::StensorDimeToSize<N>::value;
// tensor size
static constexpr unsigned short TensorSize  = tfel::math::TensorDimeToSize<N>::value;
// size of the driving variable array (STRAN)
static constexpr unsigned short GradientSize = TensorSize;
// size of the thermodynamic force variable array (STRESS)
static constexpr unsigned short ThermodynamicForceVariableSize = StensorSize;
static constexpr bool useTimeSubStepping = false;
static constexpr bool doSubSteppingOnInvalidResults = true;
static constexpr unsigned short maximumSubStepping = 0u;
static constexpr bool requiresStiffnessTensor = true;
static constexpr bool requiresUnAlteredStiffnessTensor = false;
static constexpr bool requiresThermalExpansionCoefficientTensor = false;
static constexpr unsigned short material_properties_nb = 0;
static constexpr unsigned short propertiesOffset = CastemOrthotropicOffset<castem::STANDARDFINITESTRAINBEHAVIOUR,H>::value;
static constexpr CastemSymmetryType stype = castem::ORTHOTROPIC;
}; // end of class CastemTraits

} // end of namespace castem

#endif /* __cplusplus */

#ifdef __cplusplus
extern "C"{
#endif /* __cplusplus */

MFRONT_SHAREDOBJ int
umatmonocristal_uo2_Tridimensional_setParameter(const char *const,const double);

MFRONT_SHAREDOBJ int
umatmonocristal_uo2_Tridimensional_setUnsignedShortParameter(const char *const,const unsigned short);

MFRONT_SHAREDOBJ void
umatmonocristal_uo2_setOutOfBoundsPolicy(const int);

MFRONT_SHAREDOBJ void
umatmonocristal_uo2(castem::CastemReal *const,
 castem::CastemReal *const,
 castem::CastemReal *const,
 castem::CastemReal *const,
 castem::CastemReal *const,
 castem::CastemReal *const,
 castem::CastemReal *const,
 castem::CastemReal *const,
 castem::CastemReal *const,
 castem::CastemReal *const,
 const castem::CastemReal *const,
 const castem::CastemReal *const,
 const castem::CastemReal *const,
 const castem::CastemReal *const,
 const castem::CastemReal *const,
 const castem::CastemReal *const,
 const castem::CastemReal *const,
 const castem::CastemReal *const,
 const char           *const,
 const castem::CastemInt  *const,
 const castem::CastemInt  *const,
 const castem::CastemInt  *const,
 const castem::CastemInt  *const,
 const castem::CastemReal *const,
 const castem::CastemInt  *const,
 const castem::CastemReal *const,
 const castem::CastemReal *const,
       castem::CastemReal *const,
 const castem::CastemReal *const,
 const castem::CastemReal *const,
 const castem::CastemReal *const,
 const castem::CastemInt  *const,
 const castem::CastemInt  *const,
 const castem::CastemInt  *const,
 const castem::CastemInt  *const,
 const castem::CastemInt  *const,
       castem::CastemInt  *const,
const int);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* LIB_CASTEM_MONOCRISTAL_UO2_HXX */
