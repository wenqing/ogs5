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

#include <fstream>
#include <iostream>
#include <sstream>

#include "display.h"
#include "fem_ele.h"
#include "msh_mesh.h"
#include "StringTools.h"

namespace UTL
{
VariableValues::VariableValues(
    MeshLib::CFEMesh const* mesh,
    FiniteElement::CElement* quadrature,
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

void VariableValues::interpolate()
{
    double shapefunction[27];
    for (std::vector<UTL::DataPVD>::const_iterator it = _pvd_data.begin();
         it != _pvd_data.end();
         ++it)
    {
        std::ifstream ins_vtu((*it).vtu_file_name.data(), std::ios::in);

        if (!ins_vtu.good())
        {
            Display::ScreenMessage("Can not open file %s \n",
                                   (*it).vtu_file_name.data());
            exit(1);
        }

        Display::ScreenMessage("Processing data in %s ...\n",
                               (*it).vtu_file_name.data());

        std::vector<std::string> variable_names;
        std::vector<std::vector<double> > interpolated_value_set;
        while (!ins_vtu.eof())
        {
            std::string line_buffer;
            BaseLib::safeGetline(ins_vtu, line_buffer);
            if (line_buffer.find("<PointData") != std::string::npos)
            {
                for (;;)
                {
                    BaseLib::safeGetline(ins_vtu, line_buffer);
                    if (line_buffer.find("</PointData>") != std::string::npos)
                    {
                        break;
                    }

                    if (line_buffer.find("<DataArray") != std::string::npos)
                    {
                        std::stringstream ss;
                        ss.str(line_buffer);
                        std::string string_buff;
                        std::string name;
                        // The third one holds Name="..."
                        ss >> string_buff >> string_buff >> name >> string_buff;
                        ss.clear();
                        subtractStringInQuatation(name);
                        int ncomponents = 1;
                        if (string_buff.find("NumberOfComponents") !=
                            std::string::npos)
                        {
                            subtractStringInQuatation(string_buff);
                            ncomponents = std::stoi(string_buff);
                        }

                        if (std::distance(_pvd_data.begin(), it) == 0)
                        {
                            if (ncomponents == 1)
                            {
                                variable_names.push_back(name);
                            }
                            else
                            {
                                for (int i = 1; i <= ncomponents; i++)
                                {
                                    std::ostringstream s;
                                    s << i;
                                    const std::string i_as_string(s.str());

                                    variable_names.push_back(name + "_" +
                                                             i_as_string);
                                }
                            }
                        }

                        std::vector<std::vector<double> > variable_values;
                        variable_values.reserve(_mesh->GetNodesNumber(false));
                        for (std::size_t j = 0;
                             j < _mesh->GetNodesNumber(false);
                             j++)
                        {
                            std::vector<double> component_values;
                            for (int k = 0; k < ncomponents; k++)
                            {
                                double value;
                                ins_vtu >> value;
                                component_values.push_back(value);
                            }
                            variable_values.emplace_back(component_values);
                        }

                        std::vector<double> interpolated_values;
                        interpolated_values.reserve(_specified_points.size());
                        for (std::vector<UTL::SpecifiedPoint>::const_iterator
                                 it_spt = _specified_points.begin();
                             it_spt != _specified_points.end();
                             ++it_spt)
                        {
                            const UTL::SpecifiedPoint& point_info = *it_spt;

                            MeshLib::CElem* const element =
                                point_info.element_coverred_point;
                            _quadrature->ConfigElementWithoutQuature(element);
                            _quadrature->ConfigShapefunction(
                                element->GetElementType());
                            _quadrature->setUnitCoordinates(point_info.x);
                            _quadrature->ComputeShapefct(1, shapefunction);

                            for (int k = 0; k < ncomponents; k++)
                            {
                                double value = 0.0;
                                for (std::size_t i = 0;
                                     i < element->GetNodesNumber(false);
                                     i++)
                                {
                                    value +=
                                        variable_values[element->GetNodeIndex(
                                            i)][k] *
                                        shapefunction[i];
                                }
                                interpolated_values.push_back(value);
                            }
                        }

                        interpolated_value_set.emplace_back(
                            interpolated_values);
                    }
                }
            }
        }

        if (std::distance(_pvd_data.begin(), it) == 0)
        {
            for (std::vector<UTL::SpecifiedPoint>::const_iterator it_spt =
                     _specified_points.begin();
                 it_spt != _specified_points.end();
                 ++it_spt)
            {
                const UTL::SpecifiedPoint& point_info = *it_spt;
                std::ofstream ofs(point_info.name.data(), std::ios::trunc);

                ofs << "Time";

                const std::string delim = " ";
                for (std::size_t i = 0; i < variable_names.size(); i++)
                {
                    ofs << delim << variable_names[i];
                }
                ofs << std::endl;

                ofs.setf(std::ios::scientific, std::ios::floatfield);
                ofs.precision(12);
                ofs << (*it).time;

                for (std::size_t i = 0; i < interpolated_value_set.size(); i++)
                {
                    ofs << delim
                        << interpolated_value_set[i][std::distance(
                               _specified_points.begin(), it_spt)];
                }
                ofs << std::endl;

                ofs.close();
            }
        }
        else
        {
            for (std::vector<UTL::SpecifiedPoint>::const_iterator it_spt =
                     _specified_points.begin();
                 it_spt != _specified_points.end();
                 ++it_spt)
            {
                const UTL::SpecifiedPoint& point_info = *it_spt;
                std::ofstream ofs(point_info.name.data(), std::ios::app);

                ofs.setf(std::ios::scientific, std::ios::floatfield);
                ofs.precision(12);
                ofs << (*it).time;

                const std::string delim = " ";
                for (std::size_t i = 0; i < interpolated_value_set.size(); i++)
                {
                    ofs << delim
                        << interpolated_value_set[i][std::distance(
                               _specified_points.begin(), it_spt)];
                }
                ofs << std::endl;
                ofs.close();
            }
        }

        ins_vtu.clear();
        ins_vtu.close();
    }
}

void subtractStringInQuatation(std::string& a_string)
{
    std::size_t pos = a_string.find_first_of("\"");
    a_string.replace(a_string.begin(), a_string.begin() + pos + 1, "");
    pos = a_string.find("\"");
    a_string.replace(a_string.begin() + pos, a_string.end(), "");
}

}  // namespace UTL
