/**
 * \file FEMCondition.h
 * 25/11/2010 KR inital implementation
 *
 */

#ifndef FEMCONDITION_H
#define FEMCONDITION_H

#include "GeoInfo.h"
#include "ProcessInfo.h"
#include "DistributionInfo.h"
#include "GeoObject.h"

#include <vector>

class CBoundaryCondition;
class CInitialCondition;
class CSourceTerm;

/** 
 * \brief Adapter class for handling FEM Conditions in the user Interface
 */
class FEMCondition : public GeoInfo, public ProcessInfo, public DistributionInfo
{
public:
	/// Specifier for types of FEM Conditions
	enum CondType {
		UNSPECIFIED        = 0,
		BOUNDARY_CONDITION = 1,
		INITIAL_CONDITION  = 2,
		SOURCE_TERM        = 3
	};

	FEMCondition(const std::string &geometry_name, CondType = UNSPECIFIED);
	~FEMCondition() {};

	/// Returns the type of the FEM Condition (i.e. BC, IC or ST)
	CondType getCondType() const { return _type; };

	/// Returns the value(s) for the distribution 
	const std::vector<double> getDisValue() const { return _disValue; };

	/// Returns the name of the geo-object the condition is assigned to. This object is part of the associated geometry.
	const std::string& getGeoName() const { return _geoName; };

	/// Returns the name of the associated geometry.
	const std::string& getAssociatedGeometryName() const { return _associated_geometry; };

	/// Sets a vector of values specifying the distribution.
	void setDisValue(std::vector<double> disValue) { for (size_t i=0; i<disValue.size(); i++) _disValue.push_back(disValue[i]); };

	/// Convenience method for setting a single value specifying the distribution.
	void setDisValue(double disValue) { _disValue.push_back(disValue); };

	/// Sets the name of the geo-object the condition is assigned to.
	void setGeoName(std::string geoName) { _geoName = geoName; };

protected:
	CondType _type;
	GEOLIB::GeoObject* _geoObject;
	std::string _geoName;
	std::vector<double> _disValue;
	std::string _associated_geometry;
};

/** 
 * \brief Adapter class for handling Boundary Conditions in the user Interface
 */
class BoundaryCondition : public FEMCondition
{
public:
	BoundaryCondition(const std::string &geometry_name) : FEMCondition(geometry_name, FEMCondition::BOUNDARY_CONDITION), _tim_type(0) {};
	BoundaryCondition(const CBoundaryCondition &bc, const std::string &geometry_name);
	~BoundaryCondition() {};

	size_t getTimType() const {return _tim_type; };
	void setTimType(size_t value) { _tim_type = value; };


private:
	size_t _tim_type;
};

/** 
 * \brief Adapter class for handling Initial Conditions in the user Interface
 */
class InitialCondition : public FEMCondition
{
public:
	InitialCondition(const std::string &geometry_name) : FEMCondition(geometry_name, FEMCondition::INITIAL_CONDITION) {};
	InitialCondition(const CInitialCondition &ic, const std::string &geometry_name);
	~InitialCondition() {};
};

/** 
 * \brief Adapter class for handling Source Terms in the user Interface
 */
class SourceTerm : public FEMCondition
{
public:
	SourceTerm(const std::string &geometry_name) : FEMCondition(geometry_name, FEMCondition::SOURCE_TERM), _tim_type(0) {};
	SourceTerm(const CSourceTerm &st, const std::string &geometry_name);
	~SourceTerm() {};

	size_t getTimType() const {return _tim_type; };
	void setTimType(size_t value) { _tim_type = value; };

private:
	size_t _tim_type;
};


#endif //FEMCONDITION_H
