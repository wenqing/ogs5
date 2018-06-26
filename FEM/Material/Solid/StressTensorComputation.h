/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \author Wenqing Wang
 *  \file   StressTensorComputation.h
 *  Created on June 25, 2018, 1:54 PM
 */

#pragma once

namespace StressTensorComputation
{
/// \return The first stress invariant.
double getDeviatoricStress(double* stress);
double getNorm(const double* s, const int dim);
double getSecondInvariant(const double* s1, const double* s2, const int dim);
double getThirdInvariant(const double* s1, const double* s2, const double* s3, const int dim);
}
