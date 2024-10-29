/*!
* \file   MonoCristal_UO2-generic.hxx
* \brief  This file declares the umat interface for the MonoCristal_UO2 behaviour law
* \author Luc Portelette / Thomas Helfer / Etienne Castelier
* \date   
*/

#ifndef LIB_GENERIC_MONOCRISTAL_UO2_HXX
#define LIB_GENERIC_MONOCRISTAL_UO2_HXX

#include"TFEL/Config/TFELConfig.hxx"
#include"MFront/GenericBehaviour/BehaviourData.h"

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
extern "C"{
#endif /* __cplusplus */

/*!
 * \brief rotate the gradients from the global frame to the material frame
 * \param[out] dest: gradients in the material frame
 * \param[in] src: gradients in the global frame
 * \param[in] rv: rotation matrix
 */
MFRONT_SHAREDOBJ void MonoCristal_UO2_Tridimensional_rotateGradients(mfront_gb_real* const, const mfront_gb_real* const, const mfront_gb_real* const);

/*!
 * \brief rotate an array of gradients from the global frame to the material frame
 * \param[out] dest: pointer to the gradients array in the material frame
 * \param[in] src: pointer to the gradients array in the global frame
 * \param[in] rv: rotation matrix
 * \param[in] s: number of entry of the array to be treated
 */
MFRONT_SHAREDOBJ void MonoCristal_UO2_Tridimensional_rotateArrayOfGradients(mfront_gb_real* const, const mfront_gb_real* const, const mfront_gb_real* const, const mfront_gb_size_type);

/*!
 * \brief rotate the Cauchy stress from the material  frame to the global frame
 * \param[out] dest: Cauchy stress in the global frame
 * \param[in] src: Cauchy stress in the material frame
 * \param[in] rv: rotation matrix
 */
MFRONT_SHAREDOBJ void MonoCristal_UO2_Tridimensional_rotateThermodynamicForces_CauchyStress(mfront_gb_real* const, const mfront_gb_real* const, const mfront_gb_real* const);

/*!
 * \brief rotate an array of Cauchy stresses from the material frame to the global frame
 * \param[out] dest: array of Cauchy stresses in the global frame
 * \param[in] src: array of Cauchy stresses in the material frame
 * \param[in] rv: rotation matrix
 * \param[in] s: number of entry of the array to be treated
 */
MFRONT_SHAREDOBJ void MonoCristal_UO2_Tridimensional_rotateArrayOfThermodynamicForces_CauchyStress(mfront_gb_real* const, const mfront_gb_real* const, const mfront_gb_real* const, const mfront_gb_size_type);

/*!
 * \brief rotate the second Piola-Kirchhoff stress from the material frame to the global frame
 * \param[out] dest: second Piola-Kirchhoff stress in the global frame
 * \param[in] src: second Piola-Kirchhoff stress in the material frame
 * \param[in] rv: rotation matrix
 */
MFRONT_SHAREDOBJ void MonoCristal_UO2_Tridimensional_rotateThermodynamicForces_PK2Stress(mfront_gb_real* const, const mfront_gb_real* const, const mfront_gb_real* const);

/*!
 * \brief rotate an array of second Piola-Kirchhoff stresses from the material frame to the global frame
 * \param[out] dest: array of second Piola-Kirchhoff stresses in the global frame
 * \param[out] src: array of second Piola-Kirchhoff stresses in the material frame
 * \param[in] rv: rotation matrix
 * \param[in] s: number of entry of the array to be treated
 */
MFRONT_SHAREDOBJ void MonoCristal_UO2_Tridimensional_rotateArrayOfThermodynamicForces_PK2Stress(mfront_gb_real* const, const mfront_gb_real* const, const mfront_gb_real* const, const mfront_gb_size_type);

/*!
 * \brief rotate the first Piola-Kirchhoff stress from the material frame to the global frame
 * \param[out] dest: first Piola-Kirchhoff stress in the global frame
 * \param[in] src: first Piola-Kirchhoff stress in the material frame
 * \param[in] rv: rotation matrix
 */
MFRONT_SHAREDOBJ void MonoCristal_UO2_Tridimensional_rotateThermodynamicForces_PK1Stress(mfront_gb_real* const, const mfront_gb_real* const, const mfront_gb_real* const);

/*!
 * \brief rotate an array of the first Piola-Kirchhoff stresses from the material frame to the global frame
 * \param[out] dest: array of first Piola-Kirchhoff stresses in the global frame
 * \param[out] src: array of first Piola-Kirchhoff stresses in the material frame
 * \param[in] rv: rotation matrix
 * \param[in] s: number of entry of the array to be treated
 */
MFRONT_SHAREDOBJ void MonoCristal_UO2_Tridimensional_rotateArrayOfThermodynamicForces_PK1Stress(mfront_gb_real* const, const mfront_gb_real* const, const mfront_gb_real* const, const mfront_gb_size_type);

/*!
 * \brief rotate the derivative of the Cauchy stress with respect to the deformation gradient from the material frame to the global frame
 * \param[out] dest: derivative of the Cauchy stress with respect to the deformation gradient in the global frame
 * \param[in] src: derivative of the Cauchy stress with respect to the deformation gradient in the material frame
 * \param[in] rv: rotation matrix
 */
MFRONT_SHAREDOBJ void MonoCristal_UO2_Tridimensional_rotateTangentOperatorBlocks_dsig_dF(mfront_gb_real* const, const mfront_gb_real* const, const mfront_gb_real* const);

/*!
 * \brief rotate an array of derivatives of the Cauchy stress with respect to the deformation gradient from the material frame to the global frame
 * \param[out] dest: array of derivatives of the Cauchy stress with respect to the deformation gradient in the global  frame
 * \param[in] src: array of derivatives of the Cauchy stress with respect to the deformation gradient in the material frame
 * \param[in] rv: rotation matrix
 * \param[in] s: number of entry of the array to be treated
 */
MFRONT_SHAREDOBJ void MonoCristal_UO2_Tridimensional_rotateArrayOfTangentOperatorBlocks_dsig_dF(mfront_gb_real* const, const mfront_gb_real* const, const mfront_gb_real* const, const mfront_gb_size_type);

/*!
 * \brief rotate the derivative of the first Piola-Kirchhoff stress with respect to the deformation gradient from the material frame to the global frame
 * \param[out] dest: derivative of the first Piola-Kirchhoff stress with respect to the deformation gradient in the global frame
 * \param[in] src: derivative of the first Piola-Kirchhoff stress with respect to the deformation gradient in the material frame
 * \param[in] rv: rotation matrix
 */
MFRONT_SHAREDOBJ void MonoCristal_UO2_Tridimensional_rotateTangentOperatorBlocks_dPK1_dF(mfront_gb_real* const, const mfront_gb_real* const, const mfront_gb_real* const);

/*!
 * \brief rotate an array of derivatives of the first Piola-Kirchhoff stress with respect to the deformation gradient from the material frame to the global frame
 * \param[out] dest: array of derivatives of the first Piola-Kirchhoff stress with respect to the deformation gradient in the global frame
 * \param[in] src: derivative of the first Piola-Kirchhoff stress with respect to the deformation gradient in the material frame
 * \param[in] rv: rotation matrix
 * \param[in] s: number of entry of the array to be treated
 */
MFRONT_SHAREDOBJ void MonoCristal_UO2_Tridimensional_rotateArrayOfTangentOperatorBlocks_dPK1_dF(mfront_gb_real* const, const mfront_gb_real* const, const mfront_gb_real* const, const mfront_gb_size_type);

/*!
 * \brief rotate the derivative of the second Piola-Kirchhoff stress with respect to the Green-Lagrange strain from the material frame to the global frame
 * \param[out] dest: derivative of the second Piola-Kirchhoff stress with respect to the Green-Lagrange strain in the global frame
 * \param[in] src: derivative of the second Piola-Kirchhoff stress with respect to the Green-Lagrange strain in the material frame
 * \param[in] rv: rotation matrix
 */
MFRONT_SHAREDOBJ void MonoCristal_UO2_Tridimensional_rotateTangentOperatorBlocks_dPK2_dEGL(mfront_gb_real* const, const mfront_gb_real* const, const mfront_gb_real* const);

/*!
 * \brief rotate an array of the derivatives of the second Piola-Kirchhoff stress with respect to the Green-Lagrange strain from the material frame to the global frame
 * \param[out] dest: array of derivatives of the second Piola-Kirchhoff stress with respect to the Green-Lagrange strain in the global frame
 * \param[out] src: array of derivatives of the second Piola-Kirchhoff stress with respect to the Green-Lagrange strain in the material frame
 * \param[in] rv: rotation matrix
 * \param[in] s: number of entry of the array to be treated
 */
MFRONT_SHAREDOBJ void MonoCristal_UO2_Tridimensional_rotateArrayOfTangentOperatorBlocks_dPK2_dEGL(mfront_gb_real* const, const mfront_gb_real* const, const mfront_gb_real* const, const mfront_gb_size_type);

/*!
 * \brief rotate the derivative of the Kirchhoff stress with respect to the spatial increment of the deformation gradient from the material frame to the global frame
 * \param[out] dest: derivative of the Kirchhoff stress with respect to the spatial increment of the deformation gradient in the global frame
 * \param[in] src: derivative of the Kirchhoff stress with respect to the spatial increment of the deformation gradient in the material frame
 * \param[in] rv: rotation matrix
 */
MFRONT_SHAREDOBJ void MonoCristal_UO2_Tridimensional_rotateTangentOperatorBlocks_dtau_ddF(mfront_gb_real* const, const mfront_gb_real* const, const mfront_gb_real* const);

/*!
 * \brief rotate an array of derivatives of the Kirchhoff stress with respect to the spatial increment of the deformation gradient from the material frame to the global frame
 * \param[out] dest: array of derivatives of the Kirchhoff stress with respect to the spatial increment of the deformation gradient in the global  frame
 * \param[in] src: array of derivatives of the Kirchhoff stress with respect to the spatial increment of the deformation gradient in the material frame
 * \param[in] rv: rotation matrix
 * \param[in] s: number of entry of the array to be treated
 */
MFRONT_SHAREDOBJ void MonoCristal_UO2_Tridimensional_rotateArrayOfTangentOperatorBlocks_dtau_ddF(mfront_gb_real* const, const mfront_gb_real* const, const mfront_gb_real* const, const mfront_gb_size_type);

MFRONT_SHAREDOBJ void
MonoCristal_UO2_setOutOfBoundsPolicy(const int);

MFRONT_SHAREDOBJ int
MonoCristal_UO2_Tridimensional_setParameter(const char *const,const double);

MFRONT_SHAREDOBJ int
MonoCristal_UO2_Tridimensional_setUnsignedShortParameter(const char *const,const unsigned short);

/*!
 * \param[in,out] d: material data
 */
MFRONT_SHAREDOBJ int MonoCristal_UO2_Tridimensional(mfront_gb_BehaviourData* const);

#ifdef __cplusplus
}
#endif /* __cplusplus */

#endif /* LIB_GENERIC_MONOCRISTAL_UO2_HXX */
