/*
 * \file Triangle.h
 *
 *  Created on: Mar 23, 2010
 *      Author: TF
 */

#ifndef TRIANGLE_H_
#define TRIANGLE_H_

#include <vector>
#include "Point.h"
#include "Vector3.h"

namespace GEOLIB {

/** \brief Class Triangle consists of a reference to a point vector and
 * a vector that stores the indices in the point vector.
 * A surface is composed by triangles. The class Surface stores the position
 * of pointers to the points of triangles in the m_sfc_pnt_ids vector.
 * */
class Triangle
{
public:
	/**
	 * construction of object, initialization of reference to point vector
	 */
	Triangle (std::vector<Point *> const &pnt_vec) : m_pnts(pnt_vec) {}
	/**
	 * construction of object, initialization of reference to point vector,
	 * saves the three indices describing a triangle
	 */
	Triangle (std::vector<Point *> const &pnt_vec, size_t pnt_a, size_t pnt_b, size_t pnt_c)
	: m_pnts(pnt_vec)
	{
		m_pnt_ids[0] = pnt_a;
		m_pnt_ids[1] = pnt_b;
		m_pnt_ids[2] = pnt_c;
	}

	/**
	 * saves three indices describing a triangle
	 * */
	void setTriangle (size_t pnt_a, size_t pnt_b, size_t pnt_c)
	{
		assert (pnt_a < m_pnts.size() && pnt_b < m_pnts.size() && pnt_c < m_pnts.size());
		m_pnt_ids[0] = pnt_a;
		m_pnt_ids[1] = pnt_b;
		m_pnt_ids[2] = pnt_c;
	}

	/** \brief const access operator to access the index
	 * of the i-th triangle point
	*/
	const size_t& operator[] (size_t i) const {
		assert (i < 3);
		return m_pnt_ids[i];
	}

	/** \brief access operator to access the index
	 * of the i-th triangle point
	 * */
	size_t& operator[] (size_t i) {
		assert (i < 3);
		return m_pnt_ids[i];
	}

	/**
	 * \brief const access operator to access the i-th triangle Point
	 */
	const Point* getPoint (size_t i) const {
		assert (i < 3);
		return m_pnts[m_pnt_ids[i]];
	}

	/**
	 * contains Point? assuming Triangle is in counterclockwise order
	 * from book Real-Time Collision detection p. 204
	 */
	bool containsPoint (const double *pnt) const
	{
		MATHLIB::Vector a(*(m_pnts[m_pnt_ids[0]])), b(*(m_pnts[m_pnt_ids[1]])), c(*(m_pnts[m_pnt_ids[2]]));
		MATHLIB::Vector p (pnt[0], pnt[1], pnt[2]);
		a -= p;
		b -= p;
		c -= p;

		// compute normal vectors for triangles pab and pbc
		MATHLIB::Vector u (a.Cross (b));
		MATHLIB::Vector v (b.Cross (c));
		// make sure they are both pointing in the same direction
		if (u.Dot(v) < 0.0) return false;
		// compute normal vector for triangle pca
		MATHLIB::Vector w (c.Cross (a));
		// make sure it points in the same direction as the first two
		if (u.Dot(w) < 0.0) return false;

		return true;
	}

	bool containsPoint (const Point &pnt) const
	{
		return containsPoint (pnt.getData());
	}

protected:
	/** a vector of pointers to points */
	const std::vector<Point*> &m_pnts;
	/** position of pointers to the geometric points */
	size_t m_pnt_ids[3];
};

}

#endif /* TRIANGLE_H_ */
