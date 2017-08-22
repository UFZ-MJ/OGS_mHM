/*
 * ExtractSurface.cpp
 *
 *  Created on: Jan 26, 2011
 *      Author: TF
 */

#include "ExtractSurface.h"

// Base
#include "quicksort.h"

// GEOLIB
#include "PointWithID.h"

// MATHLIB
#include "AnalyticalGeometry.h"

// MSH
#include "msh_node.h"
#include "msh_elem.h"

namespace Mesh_Group {

ExtractSurface::ExtractSurface(CFEMesh const * msh) :
	_mesh (msh)
{}

GEOLIB::Surface* ExtractSurface::extractSurface(GEOLIB::Polygon const & polygon,
		std::vector<GEOLIB::Point*>& pnts) const
{
	if (! pnts.empty()) {
		for (size_t k(0); k<pnts.size(); k++) {
			delete pnts[k];
		}
		pnts.clear();
	}

	// get all nodes of mesh
	const std::vector<Mesh_Group::CNode*>& mesh_nodes (_mesh->getNodeVector());

	// check if nodes (!projected to x-y-plane) are inside the polygon
	std::vector<size_t> id_map;
	std::vector<GEOLIB::PointWithID*> points_inside_polygon;
	const size_t number_of_mesh_nodes (mesh_nodes.size());
	for (size_t j(0); j<number_of_mesh_nodes; j++) {
		if (polygon.isPntInPolygon (mesh_nodes[j]->X(), mesh_nodes[j]->Y(), 0.0)) {
			points_inside_polygon.push_back (new GEOLIB::PointWithID (mesh_nodes[j]->X(), mesh_nodes[j]->Y(), mesh_nodes[j]->Z(), j));
		}
		// initialize id_map
		id_map.push_back (std::numeric_limits<size_t>::max());
	}

	// lexicographical sort of point
	std::vector<size_t> perm;
	for (size_t k(0); k<points_inside_polygon.size(); k++) {
		perm.push_back(k);
	}
	Quicksort<GEOLIB::PointWithID*> (points_inside_polygon, 0, points_inside_polygon.size(), perm);

	// get surface points
	double eps (sqrt(std::numeric_limits<double>::min()));
	std::vector<GEOLIB::PointWithID*> surface_pnts;
	const size_t number_points_inside_polygon (points_inside_polygon.size());
	for (size_t k(1); k<number_points_inside_polygon; k++) {
		const GEOLIB::PointWithID& p0 (*(points_inside_polygon[k-1]));
		const GEOLIB::PointWithID& p1 (*(points_inside_polygon[k]));
		if (fabs (p0[0]-p1[0]) > eps || fabs (p0[1]-p1[1]) > eps) {
			surface_pnts.push_back (points_inside_polygon[k-1]);
		}
	}
	// last point
	surface_pnts.push_back (points_inside_polygon[number_points_inside_polygon-1]);


	// update id mapping and copy surface points in vector
	for (size_t k(0); k<surface_pnts.size(); k++){
		id_map[surface_pnts[k]->getID()] = k;
		pnts.push_back ( new GEOLIB::Point ( (*(surface_pnts[k]))[0], (*(surface_pnts[k]))[1], (*(surface_pnts[k]))[2] ) );
	}

	std::cout << "found " << pnts.size () << " mesh nodes in surface" << std::endl;

	// create surface with new point vector
	GEOLIB::Surface* sfc (new GEOLIB::Surface (pnts));

	// get all elements of mesh
	const std::vector<Mesh_Group::CElem*>& msh_elem (_mesh->getElementVector());
	const size_t msh_elem_size (msh_elem.size());
	for (size_t j(0); j<msh_elem_size; j++) {
		switch (msh_elem[j]->GetElementType()) {
		case MshElemType::TRIANGLE:
		{
			// indices of nodes of the j-th element
			const vec<long>& nodes_indices (msh_elem[j]->GetNodeIndeces ());
			size_t k;
			for (k = 0; k<nodes_indices.Size(); k++) {
				if (id_map[nodes_indices[k]] == std::numeric_limits<size_t>::max()) {
					break;
				}
			}
			if (k == nodes_indices.Size()) { // all nodes of element are points in the surface
				sfc->addTriangle (id_map[nodes_indices[0]], id_map[nodes_indices[1]], id_map[nodes_indices[2]]);
			}
			break;
		}
		case MshElemType::QUAD:
		{
			// indices of nodes of the j-th element
			const vec<long>& nodes_indices (msh_elem[j]->GetNodeIndeces ());
			size_t k;
			for (k = 0; k<nodes_indices.Size(); k++) {
				if (id_map[nodes_indices[k]] == std::numeric_limits<size_t>::max()) {
					break;
				}
			}
			if (k == nodes_indices.Size()) { // all nodes of element are points in the surface
				sfc->addTriangle (id_map[nodes_indices[0]], id_map[nodes_indices[1]], id_map[nodes_indices[2]]);
				sfc->addTriangle (id_map[nodes_indices[0]], id_map[nodes_indices[2]], id_map[nodes_indices[3]]);
			}
			break;
		}
		case MshElemType::TETRAHEDRON:
		{
			// max number of nodes for a triangle of the tetrahedron is 6 (high quality element)
			int *face_nodes (new int[6]);
			for (size_t face_number(0); face_number<4; face_number++) {
				size_t n_nodes (msh_elem[j]->GetElementFaceNodes(face_number, face_nodes));

				// indices of nodes of the j-th element
				const vec<long>& nodes_indices (msh_elem[j]->GetNodeIndeces ());
				size_t k;
				for (k = 0; k<n_nodes; k++) {
					if (id_map[nodes_indices[face_nodes[k]]] == std::numeric_limits<size_t>::max()) {
						break;
					}
				}
				if (k == n_nodes) { // all nodes of element are points in the surface
					sfc->addTriangle (id_map[nodes_indices[face_nodes[0]]], id_map[nodes_indices[face_nodes[1]]], id_map[nodes_indices[face_nodes[2]]]);
				}
			}
			delete [] face_nodes;
			break;
		}
		case MshElemType::PRISM:
		{
			// max number of nodes for a surface of a prism is 8(high quality quad element)
			int *face_nodes (new int[8]);
			for (size_t face_number(0); face_number<5; face_number++) {
				size_t n_nodes (msh_elem[j]->GetElementFaceNodes(face_number, face_nodes));

				// indices of nodes of the j-th element
				const vec<long>& nodes_indices (msh_elem[j]->GetNodeIndeces ());
				size_t k;
				for (k = 0; k<n_nodes; k++) {
					if (id_map[nodes_indices[face_nodes[k]]] == std::numeric_limits<size_t>::max()) {
						break;
					}
				}
				if (k == n_nodes) { // all nodes of element are points in the surface
					if (n_nodes == 3) {
						sfc->addTriangle (id_map[nodes_indices[face_nodes[0]]], id_map[nodes_indices[face_nodes[1]]], id_map[nodes_indices[face_nodes[2]]]);
					} else {
						// case n_nodes == 4 quad -> make two triangles
						sfc->addTriangle (id_map[nodes_indices[face_nodes[0]]], id_map[nodes_indices[face_nodes[1]]], id_map[nodes_indices[face_nodes[2]]]);
						sfc->addTriangle (id_map[nodes_indices[face_nodes[0]]], id_map[nodes_indices[face_nodes[2]]], id_map[nodes_indices[face_nodes[3]]]);
					}
				}
			}
			delete [] face_nodes;
			break;
		}
		case MshElemType::HEXAHEDRON:
		{
			// max number of nodes for a surface of a hexahedron is 8 (high quality quad element)
			int *face_nodes (new int[8]);
			for (size_t face_number(0); face_number<6; face_number++) {
				size_t n_nodes (msh_elem[j]->GetElementFaceNodes(face_number, face_nodes));

				// indices of nodes of the j-th element
				const vec<long>& nodes_indices (msh_elem[j]->GetNodeIndeces ());
				size_t k;
				for (k = 0; k<n_nodes; k++) {
					if (id_map[nodes_indices[face_nodes[k]]] == std::numeric_limits<size_t>::max()) {
						break;
					}
				}
				if (k == n_nodes) { // all nodes of element are points in the surface
					// make two triangles
					sfc->addTriangle (id_map[nodes_indices[face_nodes[0]]], id_map[nodes_indices[face_nodes[1]]], id_map[nodes_indices[face_nodes[2]]]);
					sfc->addTriangle (id_map[nodes_indices[face_nodes[0]]], id_map[nodes_indices[face_nodes[2]]], id_map[nodes_indices[face_nodes[3]]]);
				}
			}
			delete [] face_nodes;
			break;
		}
		default:
			break;
		}
	}

	return sfc;
}

} // end namespace Mesh_Group
