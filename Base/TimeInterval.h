/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file TimeInterval.h
 *
 * Created on July 17, 2019, 11:07 AM
 *
 */

#ifndef TIME_INTERVAL_H
#define TIME_INTERVAL_H

#include <cassert>
#include <vector>

namespace BaseLib
{
class TimeInterval
{
public:
    TimeInterval(const double start_time_, const double end_time_)
        : _start_time(start_time_), _end_time(end_time_)
    {
        assert(_end_time > _start_time);
    }

    bool isInTimeInterval(const double time) const
    {
        if (time < _start_time)
            return false;
        if (time > _end_time)
            return false;

        return true;
    }

private:
    const double _start_time;
    const double _end_time;
};

bool isInTimeInterval(const double time,
                      std::vector<TimeInterval*> const& time_intervals);
}  // end of namespace BaseLib

#endif
