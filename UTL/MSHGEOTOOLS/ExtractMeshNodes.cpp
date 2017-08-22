/*
 * ExtractMeshNodes.cpp
 *
 *  Created on: Dec 3, 2010
 *      Author: TF
 */

#include "ExtractMeshNodes.h"

// BASELIB
#include "quicksort.h"

// GEO
#include "Point.h"

namespace Mesh_Group {

ExtractMeshNodes::ExtractMeshNodes(const CFEMesh* msh) :
	_msh (msh), _gli_pnt_offset (0)
{
}

void ExtractMeshNodes::writeMeshNodeIDs (std::ostream& os, std::ostream& gli_out, const GEOLIB::Polygon& polygon)
{
	// get all nodes of mesh
	const std::vector<Mesh_Group::CNode*>& msh_nodes (_msh->getNodeVector());

	std::vector<size_t> node_indices;

	for (size_t j(0); j<msh_nodes.size(); j++) {
		if (msh_nodes[j]->Interior()) {
			GEOLIB::Point pnt (msh_nodes[j]->X(), msh_nodes[j]->Y(), 0);
			if (polygon.isPntInPolygon(pnt)) {
				node_indices.push_back (j);
			}
		}
	}
	// write data
	for (size_t k(0); k<node_indices.size(); k++) {
		os << node_indices[k] << std::endl;
	}

	for (size_t k(0); k<node_indices.size(); k++) {
		gli_out << k + _gli_pnt_offset << " " << msh_nodes[node_indices[k]]->X() << " " <<  msh_nodes[node_indices[k]]->Y() << " " << msh_nodes[node_indices[k]]->Z() << std::endl;
	}
	_gli_pnt_offset += node_indices.size();
}

void ExtractMeshNodes::writeTopSurfaceMeshNodeIDs (std::ostream& os, std::ostream& gli_out, const GEOLIB::Polygon& polygon)
{
	// get all nodes of mesh
	const std::vector<Mesh_Group::CNode*>& msh_nodes (_msh->getNodeVector());

	std::vector<GEOLIB::PointWithID> nodes_as_points;

	for (size_t j(0); j<msh_nodes.size(); j++) {
		if (msh_nodes[j]->Interior()) {
			GEOLIB::Point pnt (msh_nodes[j]->X(), msh_nodes[j]->Y(), 0);
			if (polygon.isPntInPolygon(pnt)) {
				nodes_as_points.push_back (GEOLIB::PointWithID (msh_nodes[j]->X(), msh_nodes[j]->Y(), msh_nodes[j]->Z(), j));
			}
		}
	}

	std::vector<size_t> perm;
	for (size_t k(0); k<nodes_as_points.size(); k++) {
		perm.push_back(k);
	}
	Quicksort<GEOLIB::PointWithID> (nodes_as_points, 0, nodes_as_points.size(), perm);

	double eps (sqrt(std::numeric_limits<double>::min()));

	// write data
	for (size_t k(1); k<nodes_as_points.size(); k++) {
		const GEOLIB::PointWithID& p0 (nodes_as_points[k-1]);
		const GEOLIB::PointWithID& p1 (nodes_as_points[k]);
		if (fabs (p0[0]-p1[0]) > eps || fabs (p0[1]-p1[1]) > eps) {
			os << p0.getID() << std::endl;
		}
	}
	// write last point
	os << nodes_as_points[nodes_as_points.size()-1].getID() << std::endl;

	size_t n_nodes (0);
	gli_out.precision (14);
	for (size_t k(1); k<nodes_as_points.size(); k++) {
		const GEOLIB::PointWithID& p0 (nodes_as_points[k-1]);
		const GEOLIB::PointWithID& p1 (nodes_as_points[k]);
		if (fabs (p0[0]-p1[0]) > eps || fabs (p0[1]-p1[1]) > eps) {
			gli_out << n_nodes + _gli_pnt_offset << " " << std::scientific << p0 << " $NAME " << p0.getID() << std::endl;
			n_nodes++;
		}
	}
	// write last point
	gli_out << n_nodes + _gli_pnt_offset << " " << std::scientific << nodes_as_points[nodes_as_points.size()-1] << " $NAME " << nodes_as_points[nodes_as_points.size()-1].getID() << std::endl;
	n_nodes++;
	_gli_pnt_offset += n_nodes;
}


} // end namespace Mesh_Group
