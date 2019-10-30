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
 * Created on October 28, 2019, 4:43 PM
 */

#include "CreateVariableValues.h"

#include <fstream>
#include <iostream>
#include <sstream>

#include "ShapeFunctionPool.h"
#include "fem_ele.h"
#include "mathlib.h"
#include "msh_mesh.h"

#include "VariableValues.h"

namespace UTL
{
void subtractStringInQuatation(std::string& a_string)
{
    std::size_t pos = a_string.find_first_of("\"");
    a_string.replace(a_string.begin(), a_string.begin() + pos + 1, "");
    pos = a_string.find("\"");
    a_string.replace(a_string.begin() + pos, a_string.end(), "");
}

double computePyramidVolume(double const* const p1, double const* const p2,
                            double const* const p3, double const* const p4,
                            double const* const p5)
{
    return ComputeDetTex(p2, p4, p1, p5) + ComputeDetTex(p2, p3, p4, p5);
}

bool isPointInElement(MeshLib::CElem const& element, const double x[3])
{
    const double tol = 1.0e-9;
    switch (element.GetElementType())
    {
        case MshElemType::HEXAHEDRON:
        {
            double const* x_node[8];
            for (int i = 0; i < 8; i++)
            {
                x_node[i] = element.GetNode(i)->getData();
            }

            const double v0 = element.GetVolume();
            const double v1 = computePyramidVolume(x_node[4], x_node[5],
                                                   x_node[1], x_node[0], x) +
                              computePyramidVolume(x_node[5], x_node[6],
                                                   x_node[2], x_node[1], x) +
                              computePyramidVolume(x_node[2], x_node[6],
                                                   x_node[7], x_node[3], x) +
                              computePyramidVolume(x_node[0], x_node[3],
                                                   x_node[7], x_node[4], x) +
                              computePyramidVolume(x_node[0], x_node[1],
                                                   x_node[2], x_node[3], x) +
                              computePyramidVolume(x_node[7], x_node[6],
                                                   x_node[5], x_node[4], x);

            return (std::fabs(v0 - v1) < tol);
        }
        break;
        case MshElemType::TETRAHEDRON:
        {
            double const* x_node[4];
            for (int i = 0; i < 4; i++)
            {
                x_node[i] = element.GetNode(i)->getData();
            }

            const double v0 = element.GetVolume();
            const double v1 =
                ComputeDetTex(x_node[0], x_node[1], x_node[2], x) +
                ComputeDetTex(x_node[3], x_node[2], x_node[1], x) +
                ComputeDetTex(x_node[0], x_node[2], x_node[3], x) +
                ComputeDetTex(x_node[0], x_node[3], x_node[1], x);

            return (std::fabs(v0 - v1) < tol);
        }
        break;
        default:
        {
            std::cout << "Only HEXAHEDRON and TETRAHEDRON are supported to "
                         "identify a point in it.\n";
            exit(1);
        }
    }
    return false;
}

VariableValues* createVariableValues(
    const std::string& file_path,
    const std::string& file_name,
    FiniteElement::ShapeFunctionPool* linear_shapefunction_pool)
{
    std::ifstream ins(file_name, std::ios::in);
    if (!ins.good())
    {
        std::cout << "Can not open file " << file_name << "\n";
        exit(1);
    }

    std::string mesh_file_name;
    ins >> mesh_file_name >> std::ws;
    std::string pvd_file_name;
    ins >> pvd_file_name >> std::ws;

    std::string string_buff;
    double tol = 1.0-10;
    ins >> string_buff >> tol;

    int num_points;
    ins >> num_points >> std::ws;

    std::vector<SpecifiedPoint> specified_points;
    specified_points.reserve(num_points);

    for (int i = 0; i < num_points; i++)
    {
        int id;
        ins >> id;
        double x[3];
        ins >> x[0] >> x[1] >> x[2];
        std::string point_name;
        ins >> point_name >> point_name >> std::ws;
        specified_points.emplace_back(SpecifiedPoint(point_name, x));
    }

    std::ifstream is_mesh(file_path + "/" + mesh_file_name, std::ios::in);

    if (!is_mesh.good())
    {
        std::cout << "Cannot open mesh file " << mesh_file_name << std::endl;
        exit(1);
    }

    std::string s_buff;
    std::getline(is_mesh, s_buff);

    MeshLib::CFEMesh* mesh = NULL;
    if (s_buff.find("#FEM_MSH") != std::string::npos)
    {
        mesh = new MeshLib::CFEMesh(NULL, &mesh_file_name);
        mesh->Read(&is_mesh);
        mesh->ConstructGrid();
    }
    else
    {
        std::cout << "Cannot find keyword #FEM_MSH in" << mesh_file_name
                  << std::endl;
        exit(1);
    }

    FiniteElement::CElement* quadrature =
        new FiniteElement::CElement(mesh->GetCoordinateFlag());
    quadrature->setOrder(1);

    // Find elements and compute the local coordinates of the specified points
    const std::vector<MeshLib::CElem*>& elements = mesh->getElementVector();

    for (std::size_t i = 0; i < elements.size(); i++)
    {
        MeshLib::CElem* element = elements[i];
        bool done_ConfigElement = false;
        for (std::size_t j = 0; j < specified_points.size(); j++)
        {
            // Done
            if (specified_points[j].element_coverred_point)
                continue;

            double* x = specified_points[j].x;
            if (isPointInElement(*element, x))
            {
                if (!done_ConfigElement)
                {
                    quadrature->ConfigElementWithoutQuature(element);
                    quadrature->ConfigShapefunction(element->GetElementType());

                    done_ConfigElement = true;
                }
                specified_points[j].element_coverred_point = element;

                quadrature->UnitCoordinates(specified_points[j].x, tol);
            }
        }
    }

    // Read PVD file.
    std::ifstream is_pvd(file_path + "/" + pvd_file_name, std::ios::in);

    if (!is_pvd.good())
    {
        std::cout << "Cannot open PVD file " << pvd_file_name << std::endl;
        exit(1);
    }

    std::vector<DataPVD> pvd_data;
    while (!is_pvd.eof())
    {
        std::string line_buffer;
        std::getline(is_pvd, line_buffer);
        if (line_buffer.find("<Collection>") != std::string::npos)
        {
            for (;;)
            {
                std::getline(is_pvd, line_buffer);
                if (line_buffer.find("</Collection>") != std::string::npos)
                {
                    break;
                }

                std::stringstream ss;
                ss.str(line_buffer);
                double time;
                std::string vtu_file_name;
                while (!ss.eof())
                {
                    std::string sub_string;
                    ss >> sub_string;
                    if (sub_string.find("timestep=\"") != std::string::npos)
                    {
                        subtractStringInQuatation(sub_string);
                        time = std::atof(sub_string.c_str());
                    }
                    if (sub_string.find("file=\"") != std::string::npos)
                    {
                        subtractStringInQuatation(sub_string);
                        vtu_file_name = sub_string;
                    }
                }
                ss.clear();

                pvd_data.push_back(DataPVD(time, vtu_file_name));
            }
        }
    }

    return new VariableValues(mesh, quadrature, specified_points, pvd_data);
}
}  // namespace UTL
