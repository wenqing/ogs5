/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * Created on October 28, 2019, 1:56 PM
 */

#include "VariableValues.h"

#include "fem_ele.h"
#include "msh_mesh.h"

namespace UTL
{
VariableValues::VariableValues(
    MeshLib::CFEMesh const* mesh,
    FiniteElement::CElement const* quadrature,
    std::vector<SpecifiedPoint> const& specified_points,
    std::vector<DataPVD> const pvd_data)
    : _mesh(mesh),
      _quadrature(quadrature),
      _specified_points(specified_points),
      _pvd_data(pvd_data)
{
}

VariableValues::~VariableValues()
{
    delete _mesh;
    delete _quadrature;
}
}  // namespace UTL
