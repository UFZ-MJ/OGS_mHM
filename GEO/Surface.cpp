/*
 * Surface.cpp
 *
 *  Created on: Apr 22, 2010
 *      Author: TF
 */

#include <list>

// GEOLIB
#include "Surface.h"
#include "AxisAlignedBoundingBox.h"
#include "Polygon.h"

// MATHLIB
#include "AnalyticalGeometry.h"
#include "EarClippingTriangulation.h"

namespace GEOLIB {

Surface::Surface (const std::vector<Point*> &pnt_vec) :
	GeoObject(), m_sfc_pnts(pnt_vec), bv (new AABB())
{}

Surface::~Surface ()
{
	for (size_t k(0); k<m_sfc_triangles.size(); k++)
		delete m_sfc_triangles[k];
	delete bv;
}

void Surface::addTriangle (size_t pnt_a, size_t pnt_b, size_t pnt_c)
{
	assert (pnt_a < m_sfc_pnts.size() && pnt_b < m_sfc_pnts.size() && pnt_c < m_sfc_pnts.size());
	m_sfc_triangles.push_back (new Triangle(m_sfc_pnts, pnt_a, pnt_b, pnt_c));
	bv->update (*m_sfc_pnts[pnt_a]);
	bv->update (*m_sfc_pnts[pnt_b]);
	bv->update (*m_sfc_pnts[pnt_c]);
}

Surface* Surface::createSurface(const Polyline &ply)
{
	if (!ply.isClosed()) {
		std::cout << "Error in Surface::createSurface() - Polyline is not closed..." << std::cout;
		return NULL;
	}

	if (ply.getNumberOfPoints() > 2) {
		// create empty surface
		Surface *sfc(new Surface(ply.getPointsVec()));

		Polygon* polygon (new Polygon (ply));
		std::list<Triangle> triangles;
		MATHLIB::EarClippingTriangulation (polygon, triangles);
		std::cout << "done - " << triangles.size () << " triangles " << std::endl;

		// add Triangles to Surface
		std::list<Triangle>::const_iterator it (triangles.begin());
		while (it != triangles.end()) {
			sfc->addTriangle ((*it)[0], (*it)[1], (*it)[2]);
			it++;
		}
		return sfc;
	} else {
		std::cout << "Error in Surface::createSurface() - Polyline consists of less than three points and therefore cannot be triangulated..." << std::cout;
		return NULL;
	}

}

size_t Surface::getNTriangles () const
{
	return m_sfc_triangles.size();
}

const Triangle* Surface::operator[] (size_t i) const
{
	assert (i < m_sfc_triangles.size());
	return m_sfc_triangles[i];
}

bool Surface::isPntInBV (const double *pnt) const
{
	return bv->containsPoint (pnt);
}

bool Surface::isPntInSfc (const double *pnt) const
{
	bool nfound (true);
	for (size_t k(0); k<m_sfc_triangles.size() && nfound; k++) {
		if (m_sfc_triangles[k]->containsPoint (pnt)) nfound = false;
	}
	return !nfound;
}

} // end namespace
