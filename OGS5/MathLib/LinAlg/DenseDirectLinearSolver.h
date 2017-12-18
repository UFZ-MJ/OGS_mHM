/*
 * DenseDirectLinearSolver.h
 *
 *  Created on: Jan 7, 2011
 *      Author: TF
 */

#ifndef DENSEDIRECTLINEARSOLVER_H_
#define DENSEDIRECTLINEARSOLVER_H_

#include <DirectLinearSolver.h>

namespace MATHLIB {

class DenseDirectLinearSolver: public MATHLIB::DirectLinearSolver {
public:
	DenseDirectLinearSolver() {};
	virtual ~DenseDirectLinearSolver() {};
};

}

#endif /* DENSEDIRECTLINEARSOLVER_H_ */
