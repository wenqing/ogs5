/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file Excavation.h
 *
 * Created on June 12, 2019, 3:52 PM
 *
 */

#pragma once

#include <vector>

namespace MeshLib
{
class CElem;

class Excavation
{
public:
    Excavation(const int zone_id, const double start_position[3],
               const double end_position[3], const double start_time,
               const double end_time)
        : _zone_id(zone_id),
          _start_position{start_position[0], start_position[1],
                          start_position[2]},
          _end_position{end_position[0], end_position[1], end_position[2]},
          _start_time(start_time),
          _end_time(end_time)
    {
    }

    void deactivateElements(const double t,
                            std::vector<MeshLib::CElem*>& element_vector) const;

private:
    const int _zone_id;
    const double _start_position[3];
    const double _end_position[3];
    const double _start_time;
    const double _end_time;

    bool isInExcavatedZone(const double t,
                           const double current_position[]) const;
};
}  // namespace MeshLib
