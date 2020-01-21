/**
 * \copyright
 * Copyright (c) 2020, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 */

/*
 * DirectLinearSolver.h
 *
 *  Created on: Jan 7, 2011
 *      Author: TF
 */

#ifndef DIRECTLINEARSOLVER_H_
#define DIRECTLINEARSOLVER_H_

#include "LinearSolver.h"

namespace MathLib
{
class DirectLinearSolver : public MathLib::LinearSolver
{
public:
	DirectLinearSolver() {}
	virtual ~DirectLinearSolver() {}
};
}

#endif /* DIRECTLINEARSOLVER_H_ */
