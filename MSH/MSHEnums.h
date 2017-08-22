/**
 * \file MSHEnums.h
 * 15/11/2010 KR inital implementation
 *
 */

#ifndef MSHENUMS_H
#define MSHENUMS_H

#include <string>

/** 
 * \brief Types of mesh elements supported by OpenGeoSys.
 * Classification and associated int are identical with
 * identifier in Mesh_Group::CElem::geo_type
 */
struct MshElemType
{
	enum type {
		LINE = 1,
		QUAD = 2,
		HEXAHEDRON = 3,
		TRIANGLE = 4,
		TETRAHEDRON = 5,
		PRISM = 6,
		INVALID = -1
	};
};

/// Given a MshElemType this returns the appropriate string.
std::string MshElemType2String(const MshElemType::type t);

/// Given a string describing an element type this returns the corresponding MshElemType.
MshElemType::type String2MshElemType(const std::string &s);


#endif //MSHENUMS_H
