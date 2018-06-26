/**
 * \copyright
 * Copyright (c) 2012-2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \brief  Base on the implementation in 2003
 *  \author Wenqing Wang
 *  \file   StressTensorComputation.cpp
 *  Created on June 25, 2018, 1:54 PM
 */

#include "StressTensorComputation.h"

namespace StressTensorComputation
{
double getDeviatoricStress(double* stress)
{
	const double I1 = stress[0] + stress[1] + stress[2];
	for (int i = 0; i < 3; i++)
		stress[i] -= I1 / 3.0;

	return I1;
}

double getNorm(const double* s, const int dim)
{
	const double mean_s = (s[0] + s[1] + s[2]) / 3.0;
	const double s1 = s[0] - mean_s;
	const double s2 = s[1] - mean_s;
	const double s3 = s[2] - mean_s;
	double val = s1 * s1 + s2 * s2 + s3 * s3 + 2.0 * s[3] * s[3];
	if (dim == 3)
		val += +2.0 * s[4] * s[4] + 2.0 * s[5] * s[5];
	return val;
}

double getSecondInvariant(double const* const s, const int dim)
{
	switch (dim)
	{
		case 2:
			return s[0] * s[0] + s[1] * s[1] + s[2] * s[2] + 2.0 * s[3] * s[3];
		case 3:
			return s[0] * s[0] + s[1] * s[1] + s[2] * s[2] + 2.0 * s[3] * s[3] + 2.0 * s[4] * s[4]
			       + 2.0 * s[5] * s[5];
	}
	return 0.0; // To avoid warnings
}

double getThirdInvariant(double const* const s, const int dim)
{
	switch (dim)
	{
		case 2:
			return (  s[0] * (s[0] * s[0] + s[3] * s[3]) // s1_11*(s2_11*s3_11+s2_12*s3_21)
			        + s[3] * (s[0] * s[3] + s[3] * s[1]) // s1_12*(s2_11*s3_12+s2_12*s3_22)
			        + s[3] * (s[3] * s[0] + s[1] * s[3]) // s1_21*(s2_21*s3_11+s2_22*s3_21)
			        + s[1] * (s[3] * s[3] + s[1] * s[1]) // s1_22*(s2_21*s3_12+s2_22*s3_22)
			        + s[2] * s[2] * s[2])
			       / 3.0; // s33*s33*s33
		case 3:
			return (
			           // s1_11*(s2_11*s3_11+s2_12*s3_21+s2_13*s3_31)
			             s[0] * (s[0] * s[0] + s[3] * s[3] + s[4] * s[4])
			           // s1_12*(s2_11*s3_12+s2_12*s3_22+s2_13*s3_32)
			           + s[3] * (s[0] * s[3] + s[3] * s[1] + s[4] * s[5])
			           // s1_13*(s2_11*s3_13+s2_12*s3_23+s2_13*s3_33)
			           + s[4] * (s[0] * s[4] + s[3] * s[5] + s[4] * s[2])
			           // s1_21*(s2_21*s3_11+s2_22*s3_21+s2_23*s3_31)
			           + s[3] * (s[3] * s[0] + s[1] * s[3] + s[5] * s[4])
			           // s1_22*(s2_21*s3_12+s2_22*s3_22+s2_23*s3_32)
			           + s[1] * (s[3] * s[3] + s[1] * s[1] + s[5] * s[5])
			           // s1_23*(s2_21*s3_13+s2_22*s3_23+s2_23*s3_33)
			           + s[5] * (s[3] * s[4] + s[1] * s[5] + s[5] * s[2])
			           // s1_31*(s2_31*s3_11+s2_32*s3_21+s2_33*s3_31)
			           + s[4] * (s[4] * s[0] + s[5] * s[3] + s[2] * s[4]) // WX:bug fixed s3_31 is s3[4]
			           // s1_32*(s2_31*s3_12+s2_32*s3_22+s2_33*s3_32)
			           + s[5] * (s[4] * s[3] + s[5] * s[1] + s[2] * s[5])
			           // s1_33*(s2_31*s3_13+s2_32*s3_23+s2_33*s3_33)
			           + s[2] * (s[4] * s[4] + s[5] * s[5] + s[2] * s[2]))
			       / 3.0;
	}
	return 0.0; // To avoid warnings
}
}
