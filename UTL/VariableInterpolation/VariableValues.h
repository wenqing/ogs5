/*
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * File:   VariableValues.h
 *
 * Created on October 28, 2019, 1:56 PM
 */

#ifndef VARIABLE_VALUES_H
#define VARIABLE_VALUES_H

#include <string>
#include <vector>

namespace FiniteElement
{
class CElement;
}

namespace MeshLib
{
class CFEMesh;
class CElem;
}

namespace UTL
{
struct SpecifiedPoint
{
    SpecifiedPoint(std::string const& name_, const double x_[3])
        : element_coverred_point(NULL), name(name_)
    {
        for (int i = 0; i < 3; i++)
        {
            x[i] = x_[i];
        }
    }

    MeshLib::CElem* element_coverred_point;
    const std::string name;
    double x[3];
};

struct DataPVD
{
    DataPVD(const double time_, std::string const& vtu_file_name_)
        : time(time_), vtu_file_name(vtu_file_name_)
    {
    }

    double time;
    std::string vtu_file_name;
};

class VariableValues
{
public:
    VariableValues(MeshLib::CFEMesh const* mesh,
                   FiniteElement::CElement* quadrature,
                   std::vector<SpecifiedPoint> const& specified_points,
                   std::vector<DataPVD> const pvd_data);
    ~VariableValues();

    void interpolate();

private:
    MeshLib::CFEMesh const* _mesh;
    FiniteElement::CElement* _quadrature;

    std::vector<SpecifiedPoint> const _specified_points;
    std::vector<DataPVD> const _pvd_data;
};

void subtractStringInQuatation(std::string& a_string);

}  // namespace UTL
#endif
