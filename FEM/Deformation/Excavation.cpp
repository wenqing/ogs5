/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file Excavation.cpp
 *
 * Created on June 12, 2019, 3:52 PM
 *
 */

#include "Excavation.h"

#include <cmath>
#include <limits>

#include "MathTools.h"
#include "msh_mesh.h"

namespace MeshLib
{
void Excavation::deactivateElements(
    const double t, std::vector<MeshLib::CElem*>& element_vector) const
{
    if (t < _start_time || t > _end_time)
        return;
    for (std::size_t i = 0; i < element_vector.size(); i++)
    {
        MeshLib::CElem* element = element_vector[i];
        if (element->GetPatchIndex() != _zone_id)
            continue;

        double const* center = element->GetGravityCenter();
        const bool status = isInExcavatedZone(t, center);
        element->SetMark(status);
        element->SetExcavState(status ? 1 : -1);
    }
}

bool Excavation::isInExcavatedZone(const double t,
                                   const double current_position[]) const
{
    // For 2D
    if ((_start_time == _end_time) &&
        std::fabs(_end_time - t) < std::numeric_limits<double>::epsilon())
        return true;

    // 3D

    if (t > _end_time)  // Finished excavation.
        return true;

    const double k = (t - _start_time) / (_end_time - _start_time);
    static double current_end[3];
    for (int i = 0; i < 3; i++)
    {
        current_end[i] =
            k * (_end_position[i] - _start_position[i]) + _start_position[i];
    }

    if (!MathLib::isAcuteAngle(current_position, _start_position, current_end))
        return false;

    return MathLib::isAcuteAngle(
        current_position, current_end, _start_position);
}

}  // namespace MeshLib