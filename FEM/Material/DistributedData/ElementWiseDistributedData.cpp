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

#include "display.h"

#include <fstream>

namespace MaterialLib
{
void readData(const std::string& file_name, std::vector<double>& data)
{
    std::ifstream ins(file_name.c_str());
    if (!ins.good())
    {
        Display::ScreenMessage("Cannot open file %s ", file_name);
        exit(1);
    }
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
