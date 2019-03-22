/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file ElementWiseDistributedData.cpp
 *
 * Created on March 21, 2019, 10:23 AM
 *
 */

#include "ElementWiseDistributedData.h"

#include <cstdlib>
#include <fstream>

#include "display.h"

namespace MaterialLib
{
void readData(const std::string& file_name,
              const std::size_t number_of_elements,
              MaterialParameter::Name const parameter_name,
              std::map<MaterialParameter::Name, std::vector<double> >&
                  heterogeneous_material_data)
{
    std::ifstream ins(file_name.c_str());
    if (!ins.good())
    {
        Display::ScreenMessage("Cannot open file %s ", file_name.c_str());
        exit(1);
    }

    if (heterogeneous_material_data.find(parameter_name) ==
        heterogeneous_material_data.end())
    {
        heterogeneous_material_data[parameter_name] =
            std::vector<double>(number_of_elements);
    }

    std::vector<double>& data = heterogeneous_material_data[parameter_name];

    int size;
    ins >> size;

    for (int i = 0; i < size; i++)
    {
        int id;
        ins >> id;
        double value;
        ins >> value;
        data[id] = value;
    }
}

}  // end of namespace
