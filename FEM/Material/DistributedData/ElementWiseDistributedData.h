/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file ElementWiseDistributedData.h
 *
 * Created on March 21, 2019, 10:23 AM
 *
 */

#pragma once

#include <algorithm>
#include <string>
#include <vector>

namespace MaterialLib
{
class ElementWiseDistributedData
{
public:
    explicit ElementWiseDistributedData(const std::string& file_name)
    {
        _anisotropic_factor[0] = 1.0;
        _anisotropic_factor[1] = 1.0;
        _anisotropic_factor[2] = 1.0;

        readData(file_name);
    }
    ElementWiseDistributedData(const std::string& file_name,
                               const double anisotropic_factor[3])
    {
        _anisotropic_factor[0] = anisotropic_factor[0];
        _anisotropic_factor[1] = anisotropic_factor[1];
        _anisotropic_factor[2] = anisotropic_factor[2];
        readData(file_name);
    }

    double getParameterAtElement(const std::size_t element_id) const
    {
        return _data[element_id];
    }

    double getAnisotropicFactor(const int component) const
    {
        return _anisotropic_factor[component];
    }

private:
    double _anisotropic_factor[3];
    std::vector<double> _data;

    void readData(const std::string& file_name);
};
}  // end of namespace
