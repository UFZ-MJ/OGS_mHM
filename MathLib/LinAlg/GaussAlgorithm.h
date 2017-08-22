/*
 * GaussAlgorithm.h
 *
 *  Created on: May 6, 2010
 *      Author: TF
 */

#ifndef GAUSSALGORITHM_H_
#define GAUSSALGORITHM_H_

#include <cstddef>
#include "Matrix.h"
#include "DenseDirectLinearSolver.h"
#include "TriangularSolve.h"

namespace MATHLIB {

/**
 * This is a class for the direct solution of (dense) systems of
 * linear equations, \$fA x = b\f$. During the construction of
 * the object the matrix A is factorized in matrices L and U using
 * Gauss-Elimination with partial pivoting (rows are exchanged). In doing so
 * the entries of A change! The solution for a specific
 * right hand side is computed by the method execute().
 */
class GaussAlgorithm : public MATHLIB::DenseDirectLinearSolver {
public:
	/**
	 * A direct solver for the (dense) linear system \$fA x = b\f$.
	 * @param A at the beginning the matrix A, at the end of the construction
	 * of the object the matrix contains the factor L (without the diagonal)
	 * in the strictly lower part and the factor U in the upper part.
	 * The diagonal entries of L are all 1.0 and are not explicitly stored.
	 * Attention: the given matrix will be destroyed!
	 * @return a object of type GaussAlgorithm
	 */
	GaussAlgorithm(Matrix<double> &A);
	/**
	 * destructor, deletes the permutation
	 */
	~GaussAlgorithm();

	/**
	 * Method solves the linear system \$fA x = b\f$ (based on the LU factorization)
	 * using forward solve and backward solve
	 * @param b at the beginning the right hand side, at the end the solution
	 */
	void execute (double *b) const;

private:
	/**
	 * permute the right hand side vector according to the
	 * row permutations of the LU factorization
	 * @param b the entries of the vector b are permuted
	 */
	void permuteRHS (double* b) const;

	/**
	 * a reference to the matrix
	 */
	Matrix<double>& _mat;
	/**
	 * the size of the matrix
	 */
	size_t _n;
	/**
	 * the permutation of the rows
	 */
	size_t* _perm;
	/**
	 * the pivot index
	 */
	int* _pivot_idx;
};

} // end namespace MATHLIB

#endif /* GAUSSALGORITHM_H_ */
