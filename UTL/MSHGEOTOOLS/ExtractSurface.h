/*
 * ExtractSurface.h
 *
 *  Created on: Jan 26, 2011
 *      Author: TF
 */

#ifndef EXTRACTSURFACE_H_
#define EXTRACTSURFACE_H_

// STL
#include <vector>

// GEOLIB
#include "Point.h"
#include "Polygon.h"
#include "Surface.h"

// MSH
#include "msh_mesh.h"

namespace Mesh_Group {

class ExtractSurface {
public:
	ExtractSurface(CFEMesh const * msh);

	GEOLIB::Surface* extractSurface (GEOLIB::Polygon const & ply, std::vector<GEOLIB::Point*>& pnts) const;

private:
	CFEMesh const * _mesh;
};

}

#endif /* EXTRACTSURFACE_H_ */
