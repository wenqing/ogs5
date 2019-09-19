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
void Excavation::deactivateElementsForExcavation(const double t,
                                                 MeshLib::CFEMesh& mesh) const
{
    if (t < _start_time || t > _end_time)
        return;
    if (mesh._elements_deactivation_status.empty())
        mesh._elements_deactivation_status.resize(mesh.ele_vector.size());

    const double k = (t - _start_time) / (_end_time - _start_time);
    static double current_end[3];
    for (int i = 0; i < 3; i++)
    {
        current_end[i] =
            k * (_end_position[i] - _start_position[i]) + _start_position[i];
    }

    for (std::size_t i = 0; i < mesh.ele_vector.size(); i++)
    {
        MeshLib::CElem* element = mesh.ele_vector[i];
        if (element->GetPatchIndex() != _zone_id)
        {
            continue;
        }

        double const* center = element->GetGravityCenter();
        mesh._elements_deactivation_status[i] =
            isInExcavatedZone(t, current_end, center);
        element->_elements_deactivation_status =
            &mesh._elements_deactivation_status;
    }

    mesh.markDeactivatedNodes();
}

bool Excavation::isInExcavatedZone(const double t, const double current_end[],
                                   const double current_position[]) const
{
    // For 2D
    if ((_start_time == _end_time) &&
        std::fabs(_end_time - t) < std::numeric_limits<double>::epsilon())
        return true;

    // 3D
    if (t > _end_time)  // Finished excavation.
        return true;


    if (!MathLib::isAcuteAngle(current_position, _start_position, current_end))
        return false;

    return MathLib::isAcuteAngle(current_position, current_end,
                                 _start_position);
}

}  // namespace MeshLib