/*
 * TriangularSolve.cpp
 *
 *  Created on: May 5, 2010
 *      Author: TF
 */

#include "Matrix.h"

namespace MATHLIB {

void forwardSolve (const Matrix <double> &L, double* b)
{
	size_t m (L.getNRows());
	double t;

	for (size_t r=0; r<m; r++) {
		t = 0.0;
		for (size_t c=0; c<r; c++) {
			t += L(r,c)*b[c];
		}
		b[r] = b[r]-t;
	}
}

void backwardSolve (const Matrix <double> &U, double* y)
{
	double t;
	size_t m (U.getNRows()), n(U.getNCols());
	for (int r=m-1; r>=0; r--) {
		t = 0.0;
		for (size_t c=r+1; c<n; c++) {
			t += U(r,c)*y[c];
		}
		y[r] = (y[r]-t) / U(r,r);
	}
}

} // end namespace MATHLIB
