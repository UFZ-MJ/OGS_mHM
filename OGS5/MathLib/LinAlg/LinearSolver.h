/*
 * LinearSolver.h
 *
 *  Created on: Jan 7, 2011
 *      Author: TF
 */

#ifndef LINEARSOLVER_H_
#define LINEARSOLVER_H_

namespace MATHLIB {

/**
 * Base class for all linear solver classes.
 */
class LinearSolver {
public:
	LinearSolver() {};
	virtual ~LinearSolver() {};
};

}

#endif /* LINEARSOLVER_H_ */
