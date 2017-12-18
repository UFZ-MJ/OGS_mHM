/*
 * ExtractMeshNodes.h
 *
 *  Created on: Dec 3, 2010
 *      Author: TF
 */

#ifndef EXTRACTMESHNODES_H_
#define EXTRACTMESHNODES_H_

#include <iostream>

// MSH
#include "msh_mesh.h"
#include "msh_lib.h"

// GEO
#include "GEOObjects.h"
#include "Polygon.h"
#include "PointWithID.h"

namespace Mesh_Group {

/**
 * This class implements an algorithm to extract mesh node ids from a given (extruded) mesh.
 */
class ExtractMeshNodes {
public:
	/**
	 * constructor - take the mesh
	 * @param msh an instance of class CFEMesh
	 */
	ExtractMeshNodes(const CFEMesh* msh);
	/**
	 * This method first projects a mesh node into the x-y-plane (z=0).
	 * Then it checks if this mesh node is within the given polygon
	 * (polygon have to be located in the x-y-plane (z=0).
	 * The IDs of all (projected) mesh nodes within the polygon will
	 * be written to the stream os. For visual control the mesh nodes
	 * will be written (as gli points) to the stream gli_out.
	 * @param os output stream for IDs
	 * @param gli_out output stream for points
	 * @param polygon the polygon that have to be located in the x-y-plane (z=0)
	 */
	void writeMeshNodeIDs (std::ostream& os, std::ostream& gli_out, const GEOLIB::Polygon& polygon);
	/**
	 * This method first projects a mesh node into the x-y-plane (z=0).
	 * Then it checks if this mesh node is within the given polygon
	 * (polygon have to be located in the x-y-plane (z=0). All the mesh
	 * nodes are located within a "cylinder". The method sorts the mesh nodes
	 * lexicographical (first by x, then by y and in the end by z). The id of the
	 * mesh node with largest z coordinate and identical x and y coordinates
	 * will be written to the stream os. For visual control the associated mesh
	 * node will be written to the stream gli_out (as gli point).
	 * @param os output stream for IDs
	 * @param gli_out output stream for points
	 * @param polygon the polygon that have to be located in the x-y-plane (z=0)
	 */
	void writeTopSurfaceMeshNodeIDs (std::ostream& os, std::ostream& gli_out, const GEOLIB::Polygon& polygon);

private:
	const CFEMesh* _msh;
	/**
	 * offset for gli point index
	 */
	size_t _gli_pnt_offset;
};

}

#endif /* EXTRACTMESHNODES_H_ */
