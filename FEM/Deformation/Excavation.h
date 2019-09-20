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
class CFEMesh;

class Excavation
{
public:
    Excavation(const int zone_id, const double start_position[3],
               const double end_position[3], const double start_time,
               const double end_time, const double ambient_temperature,
               const double ambient_pore_pressure)
        : _zone_id(zone_id),
          _start_time(start_time),
          _end_time(end_time),
          _ambient_temperature(ambient_temperature),
          _ambient_pore_pressure(ambient_pore_pressure)
    {
        for (int i = 0; i < 3; i++)
        {
            _start_position[i] = start_position[i];
            _end_position[i] = end_position[i];
        }
    }

    void deactivateElementsForExcavation(const double t, MeshLib::CFEMesh& mesh,
                                         const int rank) const;

    double getAmbientTemperature() const { return _ambient_temperature; }
    double getAmbientPorePressure() const { return _ambient_pore_pressure; }

private:
    const int _zone_id;
    double _start_position[3];
    double _end_position[3];
    const double _start_time;
    const double _end_time;
    const double _ambient_temperature;
    const double _ambient_pore_pressure;

    bool isInExcavatedZone(const double t, const double current_end[],
                           const double current_position[]) const;
};
}  // namespace MeshLib
