/*
 * GaussAlgorithm.cpp
 *
 *  Created on: May 6, 2010
 *      Author: TF
 */

#include <cmath>
#include "GaussAlgorithm.h"
#include "swap.h"

namespace MATHLIB {

GaussAlgorithm::GaussAlgorithm (Matrix <double> &A) :
	_mat (A), _n(_mat.getNRows()), _perm (new size_t [_n])
{
	size_t k, i, j, nr (_mat.getNRows()), nc(_mat.getNCols());
	double l;

	for (k=0; k<nc; k++) {
		// search pivot
		double t = fabs(_mat(k, k));
		_perm[k] = k;
		for (i=k+1; i<nr; i++) {
			if (fabs(_mat(i,k)) > t) {
				t = _mat(i,k);
				_perm[k] = i;
			}
		}

		// exchange rows
		if (_perm[k] != k) {
			for (j=0; j<nc; j++) BASELIB::swap (_mat(_perm[k],j), _mat(k,j));
		}

		// eliminate
		for (i=k+1; i<nr; i++) {
			l=_mat(i,k)/_mat(k,k);
			for (j=k; j<nc; j++) {
				_mat(i,j) -= _mat(k,j) * l;
			}
			_mat(i,k) = l;
		}
	}
}

GaussAlgorithm::~GaussAlgorithm()
{
	delete [] _perm;
}

void GaussAlgorithm::execute (double *b) const
{
	permuteRHS (b);
	forwardSolve (_mat, b); // L z = b, b will be overwritten by z
	backwardSolve (_mat, b); // U x = z, b (z) will be overwritten by x
}

void GaussAlgorithm::permuteRHS (double* b) const
{
	for (size_t i=0; i<_n; i++) {
		if (_perm[i] != i) BASELIB::swap(b[i], b[_perm[i]]);
	}
}

} // end namespace MATHLIB
