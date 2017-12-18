/*
 * DirectLinearSolver.h
 *
 *  Created on: Jan 7, 2011
 *      Author: TF
 */

#ifndef DIRECTLINEARSOLVER_H_
#define DIRECTLINEARSOLVER_H_

#include <LinearSolver.h>

namespace MATHLIB {

class DirectLinearSolver: public MATHLIB::LinearSolver {
public:
	DirectLinearSolver() {};
	virtual ~DirectLinearSolver() {};
};

}

#endif /* DIRECTLINEARSOLVER_H_ */
