/*!
* \file   umatMonoCristal_UO2.cxx
* \brief  This file implements the umat interface for the MonoCristal_UO2 behaviour law
* \author Luc Portelette / Thomas Helfer / Etienne Castelier
* \date   
*/

#include<iostream>
#include<stdexcept>
#include"MFront/Castem/CastemOutOfBoundsPolicy.hxx"
#include"MFront/Castem/CastemInterface.hxx"

#include"MFront/Castem/CastemStressFreeExpansionHandler.hxx"

#include"TFEL/Material/MonoCristal_UO2.hxx"
#include"MFront/Castem/umatMonoCristal_UO2.hxx"

static tfel::material::OutOfBoundsPolicy&
umatmonocristal_uo2_getOutOfBoundsPolicy(){
static auto policy = []{
  const auto p =   castem::CastemOutOfBoundsPolicy::getCastemOutOfBoundsPolicy().  getOutOfBoundsPolicy();
  if(p.has_value()){
    return *p;
  }
  return tfel::material::None;
}();
return policy;
}

extern "C"{

MFRONT_EXPORT_SYMBOL(const char*, umatmonocristal_uo2_author, "Luc Portelette / Thomas Helfer / Etienne Castelier");

MFRONT_EXPORT_SYMBOL(const char*, umatmonocristal_uo2_date, "");

MFRONT_EXPORT_SYMBOL(const char*, umatmonocristal_uo2_description, "Loi de plasticité cristalline UO2\nglissement devier/combiné (GD/C ou GD-C)");

MFRONT_EXPORT_SYMBOL(const char*, umatmonocristal_uo2_build_id, "");

MFRONT_EXPORT_SYMBOL(const char*, umatmonocristal_uo2_mfront_ept, "umatmonocristal_uo2");

MFRONT_EXPORT_SYMBOL(const char*, umatmonocristal_uo2_tfel_version, "4.2.0-dev");

MFRONT_EXPORT_SYMBOL(const char*, umatmonocristal_uo2_unit_system, "");

MFRONT_EXPORT_SYMBOL(unsigned short, umatmonocristal_uo2_mfront_mkt, 1u);

MFRONT_EXPORT_SYMBOL(const char*, umatmonocristal_uo2_mfront_interface, "Castem");

MFRONT_EXPORT_SYMBOL(const char*, umatmonocristal_uo2_src, "Mono_UO2_CosH_Jaco.mfront");

MFRONT_EXPORT_SYMBOL(unsigned short, umatmonocristal_uo2_nModellingHypotheses, 1u);

MFRONT_EXPORT_ARRAY_OF_SYMBOLS(const char *, umatmonocristal_uo2_ModellingHypotheses, 1, MFRONT_EXPORT_ARRAY_ARGUMENTS("Tridimensional"));

MFRONT_EXPORT_SYMBOL(unsigned short, umatmonocristal_uo2_nTangentOperatorBlocks, 0u);

MFRONT_EXPORT_SYMBOL(const char * const *, umatmonocristal_uo2_TangentOperatorBlocks, nullptr);

MFRONT_EXPORT_SYMBOL(unsigned short, umatmonocristal_uo2_SymmetryType, 1u);

MFRONT_EXPORT_SYMBOL(unsigned short, umatmonocristal_uo2_ElasticSymmetryType, 1u);

MFRONT_EXPORT_SYMBOL(unsigned short, umatmonocristal_uo2_UsesGenericPlaneStressAlgorithm, 0u);

MFRONT_EXPORT_SYMBOL(unsigned short, umatmonocristal_uo2_nElasticMaterialPropertiesEntryPoints, 0u);

MFRONT_EXPORT_SYMBOL(const char * const *, umatmonocristal_uo2_ElasticMaterialPropertiesEntryPoints, nullptr);

MFRONT_EXPORT_SYMBOL(unsigned short, umatmonocristal_uo2_nLinearThermalExpansionCoefficientsEntryPoints, 0u);

MFRONT_EXPORT_SYMBOL(const char * const *, umatmonocristal_uo2_LinearThermalExpansionCoefficientsEntryPoints, nullptr);

MFRONT_EXPORT_SYMBOL(unsigned short, umatmonocristal_uo2_TemperatureRemovedFromExternalStateVariables, 1u);

MFRONT_EXPORT_SYMBOL(unsigned short, umatmonocristal_uo2_Tridimensional_UsableInPurelyImplicitResolution, 0u);

MFRONT_EXPORT_SYMBOL(unsigned short, umatmonocristal_uo2_Tridimensional_nMaterialProperties, 0u);

MFRONT_EXPORT_SYMBOL(const char * const *, umatmonocristal_uo2_Tridimensional_MaterialProperties, nullptr);

MFRONT_EXPORT_SYMBOL(unsigned short, umatmonocristal_uo2_Tridimensional_nInternalStateVariables, 25u);
MFRONT_EXPORT_ARRAY_OF_SYMBOLS(const char *, umatmonocristal_uo2_Tridimensional_InternalStateVariables, 25, MFRONT_EXPORT_ARRAY_ARGUMENTS("PlasticSlip[0]",
"PlasticSlip[1]","PlasticSlip[2]","PlasticSlip[3]","PlasticSlip[4]","PlasticSlip[5]",
"PlasticSlip[6]","PlasticSlip[7]","PlasticSlip[8]","PlasticSlip[9]","PlasticSlip[10]",
"PlasticSlip[11]","PlasticSlip[12]","PlasticSlip[13]","PlasticSlip[14]","PlasticSlip[15]",
"PlasticSlip[16]","PlasticSlip[17]","PlasticSlip[18]","PlasticSlip[19]","PlasticSlip[20]",
"PlasticSlip[21]","PlasticSlip[22]","PlasticSlip[23]","ElasticPartOfTheDeformationGradient"));

MFRONT_EXPORT_ARRAY_OF_SYMBOLS(int, umatmonocristal_uo2_Tridimensional_InternalStateVariablesTypes, 25, MFRONT_EXPORT_ARRAY_ARGUMENTS(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,3));

MFRONT_EXPORT_SYMBOL(unsigned short, umatmonocristal_uo2_Tridimensional_nExternalStateVariables, 0u);
MFRONT_EXPORT_SYMBOL(const char * const *, umatmonocristal_uo2_Tridimensional_ExternalStateVariables, nullptr);

MFRONT_EXPORT_SYMBOL(const int *, umatmonocristal_uo2_Tridimensional_ExternalStateVariablesTypes, nullptr);

MFRONT_EXPORT_SYMBOL(unsigned short, umatmonocristal_uo2_Tridimensional_nParameters, 6u);

MFRONT_EXPORT_ARRAY_OF_SYMBOLS(const char *, umatmonocristal_uo2_Tridimensional_Parameters, 6, MFRONT_EXPORT_ARRAY_ARGUMENTS("epsilon",
"theta","iterMax","minimal_time_step_scaling_factor","maximal_time_step_scaling_factor","numerical_jacobian_epsilon"));

MFRONT_EXPORT_ARRAY_OF_SYMBOLS(int, umatmonocristal_uo2_Tridimensional_ParametersTypes, 6, MFRONT_EXPORT_ARRAY_ARGUMENTS(0,0,2,0,0,0));

MFRONT_EXPORT_SYMBOL(double, umatmonocristal_uo2_Tridimensional_epsilon_ParameterDefaultValue, 1e-12);

MFRONT_EXPORT_SYMBOL(double, umatmonocristal_uo2_Tridimensional_theta_ParameterDefaultValue, 1);

MFRONT_EXPORT_SYMBOL(unsigned short, umatmonocristal_uo2_Tridimensional_iterMax_ParameterDefaultValue, 100);

MFRONT_EXPORT_SYMBOL(double, umatmonocristal_uo2_Tridimensional_minimal_time_step_scaling_factor_ParameterDefaultValue, 0.1);

MFRONT_EXPORT_SYMBOL(double, umatmonocristal_uo2_Tridimensional_maximal_time_step_scaling_factor_ParameterDefaultValue, 1.7976931348623e+308);

MFRONT_EXPORT_SYMBOL(double, umatmonocristal_uo2_Tridimensional_numerical_jacobian_epsilon_ParameterDefaultValue, 1e-13);

MFRONT_EXPORT_SYMBOL(unsigned short, umatmonocristal_uo2_Tridimensional_requiresStiffnessTensor, 1u);

MFRONT_EXPORT_SYMBOL(unsigned short, umatmonocristal_uo2_Tridimensional_requiresThermalExpansionCoefficientTensor, 0u);

MFRONT_EXPORT_SYMBOL(unsigned short, umatmonocristal_uo2_Tridimensional_ComputesInternalEnergy, 0u);

MFRONT_EXPORT_SYMBOL(unsigned short, umatmonocristal_uo2_Tridimensional_ComputesDissipatedEnergy, 0u);

MFRONT_EXPORT_SYMBOL(unsigned short, umatmonocristal_uo2_BehaviourType, 2u);

MFRONT_EXPORT_SYMBOL(unsigned short, umatmonocristal_uo2_BehaviourKinematic, 3u);

MFRONT_EXPORT_SYMBOL(unsigned short, umatmonocristal_uo2_Interface, 2u);

MFRONT_EXPORT_SYMBOL(unsigned short, umatmonocristal_uo2_nMainVariables, 1u);

MFRONT_EXPORT_SYMBOL(unsigned short, umatmonocristal_uo2_nGradients, 1u);

MFRONT_EXPORT_ARRAY_OF_SYMBOLS(int, umatmonocristal_uo2_GradientsTypes, 1, MFRONT_EXPORT_ARRAY_ARGUMENTS(3));

MFRONT_EXPORT_ARRAY_OF_SYMBOLS(const char *, umatmonocristal_uo2_Gradients, 1, MFRONT_EXPORT_ARRAY_ARGUMENTS("DeformationGradient"));

MFRONT_EXPORT_SYMBOL(unsigned short, umatmonocristal_uo2_nThermodynamicForces, 1u);

MFRONT_EXPORT_ARRAY_OF_SYMBOLS(int, umatmonocristal_uo2_ThermodynamicForcesTypes, 1, MFRONT_EXPORT_ARRAY_ARGUMENTS(1));

MFRONT_EXPORT_ARRAY_OF_SYMBOLS(const char *, umatmonocristal_uo2_ThermodynamicForces, 1, MFRONT_EXPORT_ARRAY_ARGUMENTS("Stress"));

MFRONT_SHAREDOBJ int
umatmonocristal_uo2_Tridimensional_setParameter(const char *const key,const double value){
using tfel::material::MonoCristal_UO2TridimensionalParametersInitializer;
auto& i = MonoCristal_UO2TridimensionalParametersInitializer::get();
try{
i.set(key,value);
} catch(std::runtime_error& e){
std::cerr << e.what() << std::endl;
return 0;
}
return 1;
}

MFRONT_SHAREDOBJ int
umatmonocristal_uo2_Tridimensional_setUnsignedShortParameter(const char *const key,const unsigned short value){
using tfel::material::MonoCristal_UO2TridimensionalParametersInitializer;
auto& i = MonoCristal_UO2TridimensionalParametersInitializer::get();
try{
i.set(key,value);
} catch(std::runtime_error& e){
std::cerr << e.what() << std::endl;
return 0;
}
return 1;
}

MFRONT_SHAREDOBJ void
umatmonocristal_uo2_setOutOfBoundsPolicy(const int p){
if(p==0){
umatmonocristal_uo2_getOutOfBoundsPolicy() = tfel::material::None;
} else if(p==1){
umatmonocristal_uo2_getOutOfBoundsPolicy() = tfel::material::Warning;
} else if(p==2){
umatmonocristal_uo2_getOutOfBoundsPolicy() = tfel::material::Strict;
} else {
std::cerr << "umatmonocristal_uo2_setOutOfBoundsPolicy: invalid argument\n";
}
}

static void 
umatmonocristal_uo2_base_TRIDIMENSIONAL(const castem::CastemInt *const NTENS, const castem::CastemReal *const DTIME,
const castem::CastemReal *const DROT,  castem::CastemReal *const DDSDDE,
const castem::CastemReal *const STRAN, const castem::CastemReal *const DSTRAN,
const castem::CastemReal *const TEMP,  const castem::CastemReal *const DTEMP,
const castem::CastemReal *const PROPS, const castem::CastemInt    *const NPROPS,
const castem::CastemReal *const PREDEF,const castem::CastemReal *const DPRED,
castem::CastemReal *const STATEV,const castem::CastemInt    *const NSTATV,
castem::CastemReal *const STRESS,castem::CastemReal *const PNEWDT,
castem::CastemInt  *const KINC,
const castem::StressFreeExpansionHandler& sfeh)
{
const auto op = umatmonocristal_uo2_getOutOfBoundsPolicy();
castem::CastemInterface<tfel::material::ModellingHypothesis::TRIDIMENSIONAL,tfel::material::MonoCristal_UO2>::exe(NTENS,DTIME,DROT,DDSDDE,STRAN,DSTRAN,
TEMP,DTEMP,PROPS,NPROPS,
PREDEF,DPRED,STATEV,NSTATV,
STRESS,PNEWDT,KINC,op,sfeh);
}

static void 
umatmonocristal_uo2_base_AXISYMMETRICAL(const castem::CastemInt *const, const castem::CastemReal *const,
const castem::CastemReal *const,  castem::CastemReal *const,
const castem::CastemReal *const, const castem::CastemReal *const,
const castem::CastemReal *const,  const castem::CastemReal *const,
const castem::CastemReal *const, const castem::CastemInt    *const,
const castem::CastemReal *const,const castem::CastemReal *const,
castem::CastemReal *const,const castem::CastemInt    *const,
castem::CastemReal *const,castem::CastemReal *const,
castem::CastemInt  *const KINC,
const castem::StressFreeExpansionHandler&)
{
std::cerr << "MonoCristal_UO2: unsupported modelling hypothesis\n";
*KINC = -2;
}

static void 
umatmonocristal_uo2_base_PLANESTRAIN(const castem::CastemInt *const, const castem::CastemReal *const,
const castem::CastemReal *const,  castem::CastemReal *const,
const castem::CastemReal *const, const castem::CastemReal *const,
const castem::CastemReal *const,  const castem::CastemReal *const,
const castem::CastemReal *const, const castem::CastemInt    *const,
const castem::CastemReal *const,const castem::CastemReal *const,
castem::CastemReal *const,const castem::CastemInt    *const,
castem::CastemReal *const,castem::CastemReal *const,
castem::CastemInt  *const KINC,
const castem::StressFreeExpansionHandler&)
{
std::cerr << "MonoCristal_UO2: unsupported modelling hypothesis\n";
*KINC = -2;
}

static void 
umatmonocristal_uo2_base_PLANESTRESS(const castem::CastemInt *const, const castem::CastemReal *const,
const castem::CastemReal *const,  castem::CastemReal *const,
const castem::CastemReal *const, const castem::CastemReal *const,
const castem::CastemReal *const,  const castem::CastemReal *const,
const castem::CastemReal *const, const castem::CastemInt    *const,
const castem::CastemReal *const,const castem::CastemReal *const,
castem::CastemReal *const,const castem::CastemInt    *const,
castem::CastemReal *const,castem::CastemReal *const,
castem::CastemInt  *const KINC,
const castem::StressFreeExpansionHandler&)
{
std::cerr << "MonoCristal_UO2: unsupported modelling hypothesis\n";
*KINC = -2;
}

static void 
umatmonocristal_uo2_base_GENERALISEDPLANESTRAIN(const castem::CastemInt *const, const castem::CastemReal *const,
const castem::CastemReal *const,  castem::CastemReal *const,
const castem::CastemReal *const, const castem::CastemReal *const,
const castem::CastemReal *const,  const castem::CastemReal *const,
const castem::CastemReal *const, const castem::CastemInt    *const,
const castem::CastemReal *const,const castem::CastemReal *const,
castem::CastemReal *const,const castem::CastemInt    *const,
castem::CastemReal *const,castem::CastemReal *const,
castem::CastemInt  *const KINC,
const castem::StressFreeExpansionHandler&)
{
std::cerr << "MonoCristal_UO2: unsupported modelling hypothesis\n";
*KINC = -2;
}

static void 
umatmonocristal_uo2_base_AXISYMMETRICALGENERALISEDPLANESTRAIN(const castem::CastemInt *const, const castem::CastemReal *const,
const castem::CastemReal *const,  castem::CastemReal *const,
const castem::CastemReal *const, const castem::CastemReal *const,
const castem::CastemReal *const,  const castem::CastemReal *const,
const castem::CastemReal *const, const castem::CastemInt    *const,
const castem::CastemReal *const,const castem::CastemReal *const,
castem::CastemReal *const,const castem::CastemInt    *const,
castem::CastemReal *const,castem::CastemReal *const,
castem::CastemInt  *const KINC,
const castem::StressFreeExpansionHandler&)
{
std::cerr << "MonoCristal_UO2: unsupported modelling hypothesis\n";
*KINC = -2;
}

static void 
umatmonocristal_uo2_base(const castem::CastemInt *const NTENS, const castem::CastemReal *const DTIME,
const castem::CastemReal *const DROT,  castem::CastemReal *const DDSDDE,
const castem::CastemReal *const STRAN, const castem::CastemReal *const DSTRAN,
const castem::CastemReal *const TEMP,  const castem::CastemReal *const DTEMP,
const castem::CastemReal *const PROPS, const castem::CastemInt    *const NPROPS,
const castem::CastemReal *const PREDEF,const castem::CastemReal *const DPRED,
castem::CastemReal *const STATEV,const castem::CastemInt    *const NSTATV,
castem::CastemReal *const STRESS,castem::CastemReal *const PNEWDT,
const castem::CastemInt *const NDI,
castem::CastemInt  *const KINC,
const castem::StressFreeExpansionHandler& sfeh)
{
if(*NDI==2){
	umatmonocristal_uo2_base_TRIDIMENSIONAL(NTENS,DTIME,DROT,DDSDDE,STRAN,DSTRAN,
 TEMP,DTEMP,PROPS,NPROPS,PREDEF,DPRED,
 STATEV,NSTATV,STRESS,PNEWDT,KINC,sfeh);
 } else if(*NDI==0){
	umatmonocristal_uo2_base_AXISYMMETRICAL(NTENS,DTIME,DROT,DDSDDE,STRAN,DSTRAN,
 TEMP,DTEMP,PROPS,NPROPS,PREDEF,DPRED,
 STATEV,NSTATV,STRESS,PNEWDT,KINC,sfeh);
 } else if(*NDI==-1){
	umatmonocristal_uo2_base_PLANESTRAIN(NTENS,DTIME,DROT,DDSDDE,STRAN,DSTRAN,
 TEMP,DTEMP,PROPS,NPROPS,PREDEF,DPRED,
 STATEV,NSTATV,STRESS,PNEWDT,KINC,sfeh);
 } else if(*NDI==-2){
	umatmonocristal_uo2_base_PLANESTRESS(NTENS,DTIME,DROT,DDSDDE,STRAN,DSTRAN,
 TEMP,DTEMP,PROPS,NPROPS,PREDEF,DPRED,
 STATEV,NSTATV,STRESS,PNEWDT,KINC,sfeh);
 } else if(*NDI==-3){
	umatmonocristal_uo2_base_GENERALISEDPLANESTRAIN(NTENS,DTIME,DROT,DDSDDE,STRAN,DSTRAN,
 TEMP,DTEMP,PROPS,NPROPS,PREDEF,DPRED,
 STATEV,NSTATV,STRESS,PNEWDT,KINC,sfeh);
 } else if(*NDI==14){
	umatmonocristal_uo2_base_AXISYMMETRICALGENERALISEDPLANESTRAIN(NTENS,DTIME,DROT,DDSDDE,STRAN,DSTRAN,
 TEMP,DTEMP,PROPS,NPROPS,PREDEF,DPRED,
 STATEV,NSTATV,STRESS,PNEWDT,KINC,sfeh);
 } else {
castem::CastemInterfaceExceptions::displayInvalidModellingHypothesisErrorMessage();
*KINC = -7;
}
}

MFRONT_SHAREDOBJ void
umatmonocristal_uo2(castem::CastemReal *const STRESS,
 castem::CastemReal *const STATEV,
 castem::CastemReal *const DDSDDE,
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
 const castem::CastemReal *const DTIME,
 const castem::CastemReal *const TEMP,
 const castem::CastemReal *const DTEMP,
 const castem::CastemReal *const PREDEF,
 const castem::CastemReal *const DPRED,
 const char           *const,
 const castem::CastemInt  *const NDI,
 const castem::CastemInt  *const,
 const castem::CastemInt  *const NTENS,
 const castem::CastemInt  *const NSTATV,
 const castem::CastemReal *const PROPS,
 const castem::CastemInt  *const NPROPS,
 const castem::CastemReal *const,
 const castem::CastemReal *const DROT,
       castem::CastemReal *const PNEWDT,
 const castem::CastemReal *const,
 const castem::CastemReal *const F0,
 const castem::CastemReal *const F1,
 const castem::CastemInt  *const,
 const castem::CastemInt  *const,
 const castem::CastemInt  *const,
 const castem::CastemInt  *const,
 const castem::CastemInt  *const,
       castem::CastemInt  *const KINC,
const int)
{
umatmonocristal_uo2_base(NTENS, DTIME, DROT, DDSDDE, F0, F1, TEMP, DTEMP, 
PROPS, NPROPS, PREDEF, DPRED, STATEV, NSTATV, 
STRESS, PNEWDT, NDI, KINC, nullptr);
}

} // end of extern "C"
