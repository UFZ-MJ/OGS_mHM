/*
 * MeshQualityChecker.cpp
 *
 *  Created on: Dec 8, 2010
 *      Author: TF
 */

#include "MeshQualityChecker.h"
#include "msh_elem.h"
#include <cmath>

namespace Mesh_Group {

MeshQualityChecker::MeshQualityChecker(CFEMesh const * const mesh) :
	_mesh (mesh), _static_histogramm (100, 0)
{
	if (_mesh) {
		_mesh_quality_messure.resize ((_mesh->getElementVector()).size(), 0);
	}
}

void MeshQualityChecker::check ()
{
	// get all elements of mesh
	const std::vector<Mesh_Group::CElem*>& msh_elem (_mesh->getElementVector());

	for (size_t k(0); k<msh_elem.size(); k++) {
		switch (msh_elem[k]->GetElementType()) {
		case MshElemType::LINE:
			// if you know a reasonable criterion let me know (TF)
			break;
		case MshElemType::TRIANGLE: {
			GEOLIB::Point* a (new GEOLIB::Point ((msh_elem[k]->GetNode(0))->getData ()));
			GEOLIB::Point* b (new GEOLIB::Point ((msh_elem[k]->GetNode(1))->getData ()));
			GEOLIB::Point* c (new GEOLIB::Point ((msh_elem[k]->GetNode(2))->getData ()));
			_mesh_quality_messure[k] = checkTriangle (a,b,c);
			delete a;
			delete b;
			delete c;
			break;
		}
		case MshElemType::QUAD: {
			GEOLIB::Point* a (new GEOLIB::Point ((msh_elem[k]->GetNode(0))->getData ()));
			GEOLIB::Point* b (new GEOLIB::Point ((msh_elem[k]->GetNode(1))->getData ()));
			GEOLIB::Point* c (new GEOLIB::Point ((msh_elem[k]->GetNode(2))->getData ()));
			GEOLIB::Point* d (new GEOLIB::Point ((msh_elem[k]->GetNode(3))->getData ()));
			_mesh_quality_messure[k] = checkQuad (a,b,c,d);
			delete a;
			delete b;
			delete c;
			delete d;
			break;
		}
		case MshElemType::TETRAHEDRON: {
			GEOLIB::Point* a (new GEOLIB::Point ((msh_elem[k]->GetNode(0))->getData ()));
			GEOLIB::Point* b (new GEOLIB::Point ((msh_elem[k]->GetNode(1))->getData ()));
			GEOLIB::Point* c (new GEOLIB::Point ((msh_elem[k]->GetNode(2))->getData ()));
			GEOLIB::Point* d (new GEOLIB::Point ((msh_elem[k]->GetNode(3))->getData ()));
			_mesh_quality_messure[k] = checkTetrahedron (a,b,c,d);
			delete a;
			delete b;
			delete c;
			delete d;
			break;
		}
		case MshElemType::PRISM: {
			std::vector<GEOLIB::Point *> pnts;
			for (size_t j(0); j<6; j++) {
				pnts.push_back (new GEOLIB::Point ((msh_elem[k]->GetNode(j))->getData ()));
			}
			_mesh_quality_messure[k] = checkPrism (pnts);
			for (size_t j(0); j<6; j++) {
				delete pnts[j];
			}
			break;
		}
		case MshElemType::HEXAHEDRON: {
			std::vector<GEOLIB::Point *> pnts;
			for (size_t j(0); j<8; j++) {
				pnts.push_back (new GEOLIB::Point ((msh_elem[k]->GetNode(j))->getData ()));
			}
			_mesh_quality_messure[k] = checkHexahedron (pnts);
			for (size_t j(0); j<8; j++) {
				delete pnts[j];
			}
			break;
		}
		default:
			std::cout << "MeshQualityChecker::check () check for element type " << MshElemType2String (msh_elem[k]->GetElementType()) << " not implemented" << std::endl;
		}

	}
}

void MeshQualityChecker::getHistogramm (std::vector<size_t>& histogramm) const
{
	// get all elements of mesh
	const std::vector<Mesh_Group::CElem*>& msh_elem (_mesh->getElementVector());

	const size_t msh_elem_size (msh_elem.size());
	const size_t histogramm_size (histogramm.size());
	for (size_t k(0); k<msh_elem_size; k++) {
		if (msh_elem[k]->GetElementType() != MshElemType::LINE) {
			histogramm[static_cast<size_t>(_mesh_quality_messure[k] * histogramm_size)]++;
		}
	}
}

double MeshQualityChecker::checkTriangle (GEOLIB::Point const * const a,
		GEOLIB::Point const * const b, GEOLIB::Point const * const c) const
{
	double len0 (sqrt(MATHLIB::sqrDist (b,a)));
	double len1 (sqrt(MATHLIB::sqrDist (b,c)));
	double len2 (sqrt(MATHLIB::sqrDist (a,c)));

	if (len0 < len1 && len0 < len2) {
		if (len1 < len2) {
			return len0/len2;
		} else {
			return len0/len1;
		}
	} else {
		if (len1 < len2) {
			if (len0 < len2) {
				return len1/len2;
			} else {
				return len1/len0;
			}
		} else {
			if (len0 < len1) {
				return len2/len1;
			} else {
				return len2/len0;
			}
		}
	}
}

double MeshQualityChecker::checkQuad (GEOLIB::Point const * const a, GEOLIB::Point const * const b,
		GEOLIB::Point const * const c, GEOLIB::Point const * const d) const
{
	double sqr_lengths[4] = {MATHLIB::sqrDist (b,a),
			MATHLIB::sqrDist (c,b),
			MATHLIB::sqrDist (d,c),
			MATHLIB::sqrDist (a,d)};

	// sort lengths - since this is a very small array we use bubble sort
	for (size_t i(0); i<4; i++) {
		for (size_t j(i+1); j<4; j++) {
			if (sqr_lengths[i] >= sqr_lengths[j]) std::swap (sqr_lengths[i], sqr_lengths[j]);
		}
	}

	return sqrt(sqr_lengths[0]) / sqrt(sqr_lengths[3]);
}

double MeshQualityChecker::checkTetrahedron (GEOLIB::Point const * const a, GEOLIB::Point const * const b,
		GEOLIB::Point const * const c, GEOLIB::Point const * const d) const
{
	double sqr_lengths[6] = {MATHLIB::sqrDist (b,a), MATHLIB::sqrDist (c,b),
			MATHLIB::sqrDist (c,a), MATHLIB::sqrDist (a,d),
			MATHLIB::sqrDist (b,d), MATHLIB::sqrDist (c,d)};

	// sort lengths - since this is a very small array we use bubble sort
	for (size_t i(0); i<6; i++) {
		for (size_t j(i+1); j<6; j++) {
			if (sqr_lengths[i] >= sqr_lengths[j]) std::swap (sqr_lengths[i], sqr_lengths[j]);
		}
	}

	return sqrt(sqr_lengths[0]) / sqrt(sqr_lengths[5]);
}

double MeshQualityChecker::checkPrism (std::vector<GEOLIB::Point *> const & pnts) const
{
	double sqr_lengths[9] = {MATHLIB::sqrDist (pnts[0],pnts[1]),
			MATHLIB::sqrDist (pnts[1],pnts[2]),
			MATHLIB::sqrDist (pnts[2],pnts[0]),
			MATHLIB::sqrDist (pnts[3],pnts[4]),
			MATHLIB::sqrDist (pnts[4],pnts[5]),
			MATHLIB::sqrDist (pnts[5],pnts[3]),
			MATHLIB::sqrDist (pnts[0],pnts[3]),
			MATHLIB::sqrDist (pnts[1],pnts[4]),
			MATHLIB::sqrDist (pnts[2],pnts[5])};

	// sort lengths - since this is a very small array we use bubble sort
	for (size_t i(0); i<9; i++) {
		for (size_t j(i+1); j<9; j++) {
			if (sqr_lengths[i] >= sqr_lengths[j]) std::swap (sqr_lengths[i], sqr_lengths[j]);
		}
	}

	return sqrt(sqr_lengths[0]) / sqrt(sqr_lengths[8]);
}

double MeshQualityChecker::checkHexahedron (std::vector<GEOLIB::Point *> const & pnts) const
{
	double sqr_lengths[12] = {MATHLIB::sqrDist (pnts[0],pnts[1]),
			MATHLIB::sqrDist (pnts[1],pnts[2]),
			MATHLIB::sqrDist (pnts[2],pnts[3]),
			MATHLIB::sqrDist (pnts[3],pnts[0]),
			MATHLIB::sqrDist (pnts[4],pnts[5]),
			MATHLIB::sqrDist (pnts[5],pnts[6]),
			MATHLIB::sqrDist (pnts[6],pnts[7]),
			MATHLIB::sqrDist (pnts[7],pnts[4]),
			MATHLIB::sqrDist (pnts[0],pnts[4]),
			MATHLIB::sqrDist (pnts[1],pnts[5]),
			MATHLIB::sqrDist (pnts[2],pnts[6]),
			MATHLIB::sqrDist (pnts[3],pnts[7])};

	// sort lengths - since this is a very small array we use bubble sort
	for (size_t i(0); i<12; i++) {
		for (size_t j(i+1); j<12; j++) {
			if (sqr_lengths[i] >= sqr_lengths[j]) std::swap (sqr_lengths[i], sqr_lengths[j]);
		}
	}

	return sqrt(sqr_lengths[0]) / sqrt(sqr_lengths[11]);
}

}
