/*
 * TriangularSolve.h
 *
 *  Created on: May 6, 2010
 *      Author: TF
 */

#ifndef TRIANGULARSOLVE_H_
#define TRIANGULARSOLVE_H_

namespace MATHLIB {

/**
 * solves the \$fn \times n\$f triangular linear system \$fL \cdot y = b\f$,
 * assumes \f$L_{ii} = 1.0\f$, \f$i=1,...,n\f$, \f$b\f$ is destroyed
 * @param L the lower triangular matrix
 * @param b at beginning the right hand side vector, at the end the solution vector
 */
void forwardSolve (const Matrix <double> &L, double* b);

/**
 * solves the \$fn \times n\$f triangular linear system \$fU \cdot x=y\$f,
 * \$fU\$f, where \$fU\$f is a upper triangular matrix.
 * @param U upper triangular matrix
 * @param y at beginning the right hand side, at the end the solution
 */
void backwardSolve (const Matrix <double> &U, double* y);

} // end namespace MATHLIB

#endif /* TRIANGULARSOLVE_H_ */
