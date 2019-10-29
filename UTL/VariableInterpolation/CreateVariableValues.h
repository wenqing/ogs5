/*
 * \file
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * Created on October 28, 2019, 4:43 PM
 */

#ifndef CREATE_VARIABLE_VALUES_H
#define CREATE_VARIABLE_VALUES_H

#include <string>

namespace FiniteElement
{
class ShapeFunctionPool;
}

namespace UTL
{
class VariableValues;
VariableValues* createVariableValues(
    const std::string& file_path,
    const std::string& file_name,
    FiniteElement::ShapeFunctionPool* linear_shapefunction_pool);
}  // namespace UTL

#endif
