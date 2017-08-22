/**
 * \file FEMCondition.cpp
 * 25/11/2010 KR inital implementation
 *
 */

#include "FEMCondition.h"

#include "rf_bc_new.h"
#include "rf_ic_new.h"
#include "rf_st_new.h"

FEMCondition::FEMCondition(const std::string &geometry_name, CondType t) 
: _type(t), _geoObject(NULL), _geoName("[unspecified]"), _associated_geometry(geometry_name)
{
	this->setProcessType(INVALID_PROCESS);
	this->setProcessPrimaryVariable(INVALID_PV);
	this->setGeoType(GEOLIB::INVALID);
	this->setProcessDistributionType(FiniteElement::INVALID_DIS_TYPE);
}


BoundaryCondition::BoundaryCondition(const CBoundaryCondition &bc, const std::string &geometry_name)
: FEMCondition(geometry_name, FEMCondition::BOUNDARY_CONDITION)
{
	this->setProcessType(bc.getProcessType());
	this->setProcessPrimaryVariable(bc.getProcessPrimaryVariable());
	this->setGeoType(bc.getGeoType());
	this->setGeoName(bc.getGeoName());
	this->setProcessDistributionType(bc.getProcessDistributionType());
}

InitialCondition::InitialCondition(const CInitialCondition &ic, const std::string &geometry_name)
: FEMCondition(geometry_name, FEMCondition::INITIAL_CONDITION)
{
	this->setProcessType(ic.getProcessType());
	this->setProcessPrimaryVariable(ic.getProcessPrimaryVariable());
	this->setGeoType(ic.getGeoType());
	this->setGeoName("[unspecified]");//ic.getGeoName());
	this->setProcessDistributionType(ic.getProcessDistributionType());
}

SourceTerm::SourceTerm(const CSourceTerm &st, const std::string &geometry_name)
: FEMCondition(geometry_name, FEMCondition::SOURCE_TERM)
{
	this->setProcessType(st.getProcessType());
	this->setProcessPrimaryVariable(st.getProcessPrimaryVariable());
	this->setGeoType(st.getGeoType());
	this->setGeoName(st.getGeoName());
	this->setProcessDistributionType(st.getProcessDistributionType());
}
