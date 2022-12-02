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

    // For some members of std containers, e.g. std::vector::erase();
    SpecifiedPoint& operator=(const SpecifiedPoint& specified_point);

    MeshLib::CElem* element_coverred_point;
    std::string name;
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

struct Excavation
{
    /*
    Excavation(const int direction,
    const double start_position,
    const double depth,
    const double start_time,
    const double end_tim,
    std::vector<std::size_t> const excavated_dom_ids) :
    direction(excv.direction),
    start_position(excv.start_position),
    depth(excv.depth),
    start_time(excv.start_time),
    end_tim(excv.end_tim),
    excavated_dom_ids(excv.excavated_dom_ids)
    {}


    Excavation(const Excavation& excv) :
    direction(excv.direction),
    start_position(excv.start_position),
    depth(excv.depth),
    start_time(excv.start_time),
    end_tim(excv.end_tim),
    excavated_dom_ids(excv.excavated_dom_ids)
    {}

    */

    int direction = 0;  // 0, 1, 2 for x,y,z
    double start_position = 0;
    double depth = 0;
    double start_time = 0;
    double end_tim = 0;
    std::vector<std::size_t> excavated_dom_ids;
};

struct DeactivatedDoms
{
    std::size_t dom_id;
    double start_time;
    double end_time;
};

class VariableValues
{
public:
    VariableValues(MeshLib::CFEMesh const* mesh,
                   FiniteElement::CElement* quadrature,
                   std::vector<SpecifiedPoint> const& specified_points,
                   std::vector<DataPVD> const pvd_data,
                   Excavation const& excavation,
                   std::vector<DeactivatedDoms> const& deactivated_doms);
    ~VariableValues();

    void interpolate(const std::string& output_path);

private:
    MeshLib::CFEMesh const* _mesh;
    FiniteElement::CElement* _quadrature;

    std::vector<SpecifiedPoint> const _specified_points;
    std::vector<DataPVD> const _pvd_data;

    Excavation const _excavation;
    std::vector<DeactivatedDoms> const _deactivated_doms;
};

void subtractStringInQuatation(std::string& a_string);

}  // namespace UTL
#endif
