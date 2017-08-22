/**
 * \file MSHEnums.cpp
 * 16/11/2010 KR inital implementation
 *
 */

#include "MSHEnums.h"


std::string MshElemType2String(const MshElemType::type t)
{
	if (t == MshElemType::LINE)			return "line";
	if (t == MshElemType::QUAD)			return "quad";
	if (t == MshElemType::HEXAHEDRON)	return "hex";
	if (t == MshElemType::TRIANGLE)		return "tri";
	if (t == MshElemType::TETRAHEDRON)	return "tet";
	if (t == MshElemType::PRISM)		return "pris";
	return "none";
};

MshElemType::type String2MshElemType(const std::string &s)
{
	if (s.compare("line") == 0) return MshElemType::LINE;
	if (s.compare("quad") == 0) return MshElemType::QUAD;
	if (s.compare("hex")  == 0) return MshElemType::HEXAHEDRON;
	if (s.compare("tri")  == 0) return MshElemType::TRIANGLE;
	if (s.compare("tet")  == 0) return MshElemType::TETRAHEDRON;
	if (s.compare("pris") == 0) return MshElemType::PRISM;
	return MshElemType::INVALID;
};


