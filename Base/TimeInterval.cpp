/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file TimeInterval.cpp
 *
 * Created on Oct 10, 2019, 14:07 AM
 *
 */

#include "TimeInterval.h"

namespace BaseLib
{
bool isInTimeInterval(const double time,
                      std::vector<TimeInterval*> const& time_intervals)
{
    // No period defined. That means the time is always in period.
    if (time_intervals.empty())
        return true;

    for (std::size_t i = 0; i < time_intervals.size(); i++)
    {
        if (time_intervals[i]->isInTimeInterval(time))
        {
            return true;
        }
    }
    return false;
}
}  // end of namespace BaseLib
