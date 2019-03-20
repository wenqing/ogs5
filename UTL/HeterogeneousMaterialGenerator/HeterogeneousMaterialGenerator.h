/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * \file HeterogeneousMaterialGenerator.h
 *
 * by Wenqing Wang
 * Created on March 18, 2019, 3:18 PM
 *
 */

#pragma once

#include <array>
#include <cmath>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "MSH/msh_mesh.h"

namespace MeshLib
{
/**
 *   It generates element wise heterogeneous material parameters for
 *    task 4, DECOVALEX 2019 by using the Gaussian distribution
 */
class HeterogeneousMaterialGenerator : public CFEMesh
{
public:
    HeterogeneousMaterialGenerator(
        std::string* file_name, const std::string variable_name,
        std::vector<std::unique_ptr<std::normal_distribution<double> > >
            gaussian,
        std::vector<std::array<double, 2> >&& ranges)
        : CFEMesh(NULL, file_name),
          _variable_name(variable_name),
          _guassian(std::move(gaussian)),
          _ranges(std::move(ranges))

    {
    }

    ~HeterogeneousMaterialGenerator(){};

    void generate(const std::string& base_file_name);

private:
    const std::string _variable_name;
    std::vector<std::unique_ptr<std::normal_distribution<double> > > _guassian;
    std::vector<std::array<double, 2> > const _ranges;
};

std::unique_ptr<HeterogeneousMaterialGenerator>
createHeterogeneousMaterialGenerator(const std::string file_name);
}  // end of namespace MeshLib
