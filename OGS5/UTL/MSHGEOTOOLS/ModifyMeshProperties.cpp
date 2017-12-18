/*
 * ModifyMeshProperties.cpp
 *
 *  Created on: Dec 8, 2010
 *      Author: fischeth
 */

#include "ModifyMeshProperties.h"

// GEO
#include "Point.h"
#include "Polygon.h"

// MSH
#include "msh_node.h"
#include "msh_elem.h"
#include "msh_mesh.h"

// MATHLIB
#include "AnalyticalGeometry.h"

// STL
#include <fstream>

namespace Mesh_Group {

ModifyMeshProperties::ModifyMeshProperties(CFEMesh* msh) :
	_mesh (msh)
{}

ModifyMeshProperties::~ModifyMeshProperties()
{}

void ModifyMeshProperties::setMaterial (const GEOLIB::Polygon& polygon, size_t mat_id)
{
	// get all nodes of mesh
	const std::vector<Mesh_Group::CNode*>& msh_nodes (_mesh->getNodeVector());

	// *** rotate polygon to xy_plane
	// 1 copy all points
	std::vector<GEOLIB::Point*> polygon_points;
	for (size_t k(0); k<polygon.getNumberOfPoints(); k++) {
		polygon_points.push_back (new GEOLIB::Point (*(polygon[k])));
	}
	// 2 rotate points
	MATHLIB::Vector plane_normal_polygon(0.0, 0.0, 0.0);
	double d_polygon (0.0);
	MATHLIB::getNewellPlane (polygon_points, plane_normal_polygon, d_polygon);

//	std::cout << "plane normal: " << plane_normal_polygon << std::endl;
	MATHLIB::rotatePointsToXY(plane_normal_polygon, polygon_points);

	// 3 create new polygon
	GEOLIB::Polyline rot_polyline (polygon_points);
	for (size_t k(0); k<polygon.getNumberOfPoints(); k++) {
		rot_polyline.addPoint (k);
	}
	GEOLIB::Polygon rot_polygon (rot_polyline);

//	std::cout << "Polygon: " << std::endl;
//	for (size_t k(0); k<polygon.getNumberOfPoints(); k++) {
//		std::cout << k << ": " << *(polygon[k]) << std::endl;
//	}
//	std::cout << std::endl;
//	std::cout << "rotiertes Polygon: " << std::endl;
//	for (size_t k(0); k<rot_polygon.getNumberOfPoints(); k++) {
//		std::cout << k << ": " << *(rot_polygon[k]) << std::endl;
//	}
//	std::cout << std::endl;

	// *** rotate mesh nodes to xy-plane
	// 1 copy all mesh nodes to GEOLIB::Points
	std::vector<GEOLIB::Point*> mesh_nodes_as_points;
	for (size_t j(0); j<msh_nodes.size(); j++) {
		mesh_nodes_as_points.push_back (new GEOLIB::Point (msh_nodes[j]->X(), msh_nodes[j]->Y(), msh_nodes[j]->Z()));
	}
	// 2 rotate the Points
	MATHLIB::rotatePointsToXY(plane_normal_polygon, mesh_nodes_as_points);

	// get all elements of mesh
	const std::vector<Mesh_Group::CElem*>& msh_elem (_mesh->getElementVector());

	// *** perform search and modify mesh
	const size_t msh_elem_size (msh_elem.size());
	for (size_t j(0); j<msh_elem_size; j++) {
		// indices of nodes of the j-th element
		const vec<long>& nodes_indices (msh_elem[j]->GetNodeIndeces ());
		size_t k;
		for (k = 0; k<nodes_indices.Size(); k++) {
			if (! rot_polygon.isPntInPolygon(*(mesh_nodes_as_points[nodes_indices[k]]))) {
				break;
			}
		}

		if (k == nodes_indices.Size()) {
			msh_elem[j]->setPatchIndex (mat_id);
		}
	}

	for (size_t k(0); k<rot_polygon.getNumberOfPoints(); k++) {
		delete rot_polygon[k];
	}
	for (size_t j(0); j<mesh_nodes_as_points.size(); j++) {
		delete mesh_nodes_as_points[j];
	}
}

}
