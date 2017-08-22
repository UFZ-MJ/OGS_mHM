/*
 * MeshNodesAlongPolyline.h
 *
 *  Created on: Aug 9, 2010
 *      Author: TF
 */

#ifndef MESHNODESALONGPOLYLINE_H_
#define MESHNODESALONGPOLYLINE_H_

// GEOLIB
#include "Polyline.h"

namespace Mesh_Group {

// forward declaration
class CFEMesh;

class MeshNodesAlongPolyline {
public:
	MeshNodesAlongPolyline(const GEOLIB::Polyline* ply, const CFEMesh* mesh);
	const std::vector<size_t>& getNodeIDs () const;
	const GEOLIB::Polyline* getPolyline () const;
	size_t getNumberOfLinearNodes () const;
    const std::vector<double>& getDistOfProjNodeFromPlyStart() const;

private:
	const GEOLIB::Polyline* _ply;
	const CFEMesh* _mesh;
	size_t _linear_nodes;
	std::vector<size_t> _msh_node_ids;
	std::vector<double> _dist_of_proj_node_from_ply_start;
};

}

#endif /* MESHNODESALONGPOLYLINE_H_ */
