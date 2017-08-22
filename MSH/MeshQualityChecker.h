/*
 * MeshQualityChecker.h
 *
 *  Created on: Dec 8, 2010
 *      Author: TF
 */

#ifndef MESHQUALITYCHECKER_H_ 
#define MESHQUALITYCHECKER_H_

#include <vector>

// MSH
#include "msh_mesh.h"

namespace Mesh_Group {

class MeshQualityChecker {
public:
	MeshQualityChecker(CFEMesh const * const mesh);

	void check ();
	const std::vector<double>& getMeshQuality () const { return _mesh_quality_messure; }
	void getHistogramm (std::vector<size_t>& histogramm) const;

private:
	double checkTriangle (GEOLIB::Point const * const a, GEOLIB::Point const * const b, GEOLIB::Point const * const c) const;
	double checkQuad (GEOLIB::Point const * const a, GEOLIB::Point const * const b, GEOLIB::Point const * const c, GEOLIB::Point const * const d) const;
	double checkTetrahedron (GEOLIB::Point const * const a, GEOLIB::Point const * const b, GEOLIB::Point const * const c, GEOLIB::Point const * const d) const;
	double checkPrism (std::vector<GEOLIB::Point *> const & pnts) const;
	double checkHexahedron (std::vector<GEOLIB::Point *> const & pnts) const;
	CFEMesh const * const _mesh;
	std::vector<double> _mesh_quality_messure;
	std::vector<size_t> _static_histogramm;
};

}

#endif /* MESHQUALITYCHECKER_H_ */
