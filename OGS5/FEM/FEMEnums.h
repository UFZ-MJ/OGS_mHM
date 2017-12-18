/**
 * \file FEMEnums.h
 * 31/08/2010 KR inital implementation
 *
 */

#ifndef FEMENUMS_H
#define FEMENUMS_H

#include <string>
#include <limits>

/** Types of physical processes supported by OpenGeoSys
 * if you change this enum, make sure you apply the changes to
 * the functions convertPorcessType(), convertProcessTypeToString(),
 *  isFlowProcess() and isDeformationProcess()
 */
enum ProcessType
{
   INVALID_PROCESS = 0,                           //!< INVALID_PROCESS
   AIR_FLOW,                                      //!< AIR_FLOW
   /// M process, single/multi-phase flow
   DEFORMATION,                                   //!< DEFORMATION
   DEFORMATION_DYNAMIC,                           //!< ...
   /// C process, single/multi-phase flow
   DEFORMATION_FLOW,                              //!< DEFORMATION_FLOW
   /// H process, incompressible flow
   GROUNDWATER_FLOW,                              //!< GROUNDWATER_FLOW
   /// T process, single/multi-phase flow
   HEAT_TRANSPORT,                                //!< HEAT_TRANSPORT
   FLUID_FLOW,
   FLUID_MOMENTUM,                                // BC only
   FLUX,
   /// H process, incompressible flow
   LIQUID_FLOW,                                   //!< LIQUID_FLOW
   MASS_TRANSPORT,                                //!< MASS_TRANSPORT
   MULTI_PHASE_FLOW,                              //!< MULTI_PHASE_FLOW
   NO_PCS,                                        //!< NO_PCS
   /// H process, incompressible flow
   OVERLAND_FLOW,                                 //!< OVERLAND_FLOW
   PS_GLOBAL,                                     //!< PS_GLOBAL
   RANDOM_WALK,                                   //!< RANDOM_WALK
   /// H process, incompressible flow
   RICHARDS_FLOW,                                 //!< RICHARDS_FLOW
   /// H2 process, compressible flow
   TWO_PHASE_FLOW                                 //!< TWO_PHASE_FLOW
};

/**
 * convert the given string into the appropriate enum value
 * @param pcs_type_string string describing a process type
 * @return enum value describing process type
 */
ProcessType convertProcessType ( const std::string& pcs_type_string );

/**
 * convert the given enum value into the appropriate string
 * @param pcs_type process type described by the enum ProcessType
 * @return string describing the process type
 */
std::string convertProcessTypeToString ( ProcessType pcs_type );

/**
 * checks if the given pcs_type variable corresponds to a flow type of the enum ProcessType
 * @param pcs_type value of enum ProcessType
 * @return true if pcs_type describes a flow process, else false
 */
bool isFlowProcess (ProcessType pcs_type);

/**
 * checks if the given pcs_type variable corresponds to a deformation type of the enum ProcessType
 * @param pcs_type value of enum ProcessType
 * @return true if pcs_type describes a deformation process, else false
 */
bool isDeformationProcess (ProcessType pcs_type);

/**
 * \enum contains all values for primary variables actually handled by OGS
 */
enum PrimaryVariable
{
   INVALID_PV  = 0,                               //!< INVALID_PV
   /// Flow (phase)
   PRESSURE,                                      //!< PRESSURE
   PRESSURE2,                                     //!< PRESSURE2
   PRESSURE_RATE1,                                // OUT
   SATURATION,                                    //!< SATURATION
   SATURATION2,                                   //!< SATURATION2
   /// Heat transport
   TEMPERATURE,                                   //!< TEMPERATURE
   /// Deformation
   DISPLACEMENT_X,                                //!< DISPLACEMENT_X
   /// Deformation
   DISPLACEMENT_Y,                                //!< DISPLACEMENT_Y
   /// Deformation
   DISPLACEMENT_Z,                                //!< DISPLACEMENT_Z
   /// Mass transport
   CONCENTRATION,                                 //!< CONCENTRATION
   HEAD,                                          //!< HEAD
   VELOCITY_DM_X,                                 //!< VELOCITY_DM_X
   VELOCITY_DM_Y,                                 //!< VELOCITY_DM_Y
   VELOCITY_DM_Z,                                 //!< VELOCITY_DM_Z
   VELOCITY1_X,
   VELOCITY1_Y,
   VELOCITY1_Z,
   STRESS_XX,                                     // IC
   STRESS_XY,                                     // IC
   STRESS_XZ,                                     // IC
   STRESS_YY,                                     // IC
   STRESS_YZ,                                     // IC
   STRESS_ZZ,                                     // IC
   STRAIN_XX,                                     // Output
   STRAIN_XY,                                     // Output
   STRAIN_XZ,                                     // Output
   STRAIN_YY,                                     // Output
   STRAIN_YZ,                                     // Output
   STRAIN_ZZ,                                     // Output
   STRAIN_PLS,                                    // Output
   ACCELERATION_X1,                               //!< ACCELERATION_X1
   ACCELERATION_Y1,                               //!< ACCELERATION_Y1
   ACCELERATION_Z1,                               //!< ACCELERATION_Z1
   EXCAVATION,                                    // ST
};

/**
 * Converts the given string into the appropriate enum value.
 * @param pcs_pv_string string describing the primary variable
 * @return enum value describing the primary variable of the process
 */
                                                  //!< PrimaryVariable
PrimaryVariable convertPrimaryVariable ( const std::string& pcs_pv_string );

/**
 * Converts the given enum value into the appropriate string.
 * @param pcs_pv primary variable described by the enum ProcessType
 * @return string describing the process type
 */
std::string convertPrimaryVariableToString ( PrimaryVariable pcs_pv );

namespace FiniteElement {

enum DistributionType {
	INVALID_DIS_TYPE = 0, ANALYTICAL, // ST
	AVERAGE,
	CONSTANT, // IC, BC, ST
	CONSTANT_GEO,
	GRADIENT, // IC
	RESTART, // IC
	LINEAR, // BC, ST
	POINT, // BC
	CONSTANT_NEUMANN, // ST
	LINEAR_NEUMANN, // ST
	NORMALDEPTH, // ST
	CRITICALDEPTH, // ST
	GREEN_AMPT, // ST
	SYSTEM_DEPENDENT, // ST
	PRECIPITATION,
	DIRECT       //WW
// Sort of Neumann BC //WW
};

/**
 * Converts the given string into the appropriate enum value.
 * @param pcs_pv_string string describing the primary variable
 * @return enum value describing the primary variable of the process
 */
DistributionType convertDisType(const std::string& dis_type_string);

/**
 * Converts the given enum value into the appropriate string.
 * @param pcs_pv primary variable described by the enum ProcessType
 * @return string describing the process type
 */
std::string convertDisTypeToString(DistributionType dis_type);
} // end namespace FiniteElement

#endif                                            //FEMENUMS_H
