/**
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *  \file HeterogeneousMaterialGenerator.cpp
 *
 * Created on March 18, 2019, 3:18 PM
 *
 */

#include "HeterogeneousMaterialGenerator.h"

#include <iostream>
#include <fstream>
#include <limits>

#include "display.h"

namespace MeshLib
{
std::unique_ptr<HeterogeneousMaterialGenerator>
createHeterogeneousMaterialGenerator(const std::string file_name)
{
    std::ifstream ins(file_name.c_str());
    if (!ins.good())
    {
        Display::ScreenMessage("Cannot open file %s ", file_name);
        exit(1);
    }

    //
    std::vector<std::unique_ptr<std::normal_distribution<double> > > guassian;
    std::vector<std::array<double, 2> > ranges;

    std::string variable_name;
    std::string s_buff;
    while (!ins.eof())
    {
        ins >> s_buff;
        if (s_buff.find("</data>") != std::string::npos)
            break;

        if (s_buff.find("<variable_name>") != std::string::npos)
        {
            ins >> variable_name;
        }

        if (s_buff.find("<gaussian>") != std::string::npos)
        {
            double mean;
            ins >> mean;
            double deviation;
            ins >> deviation;
            guassian.push_back(
                std::make_unique<std::normal_distribution<double> >(mean,
                                                                    deviation));
            double low;
            ins >> low;
            double high;
            ins >> high;

            assert(low < high);

            ranges.push_back(std::array<double, 2>{{low, high}});
        }
    }

    std::string fname = file_name;
    return std::make_unique<HeterogeneousMaterialGenerator>(
        &fname, variable_name, std::move(guassian), std::move(ranges));
}

void HeterogeneousMaterialGenerator::generate(const std::string& base_file_name)
{
    std::cout << "Generating random " << _variable_name << std::endl;

    std::random_device rd{};
    std::mt19937 gen{rd()};

    const std::string f_name_mat = base_file_name + "_" + ".txt";
    std::ofstream os_mat(f_name_mat.c_str(), std::ios::trunc);
    os_mat.setf(std::ios::scientific, std::ios::floatfield);
    os_mat.precision(12);

    const std::string f_name_vtu_cell = base_file_name + "_" + ".vtu.cell";
    std::ofstream os_vtu(f_name_vtu_cell.c_str(), std::ios::trunc);
    os_vtu.setf(std::ios::scientific, std::ios::floatfield);
    os_vtu.precision(12);

    const std::string delim = " ";
    std::vector<double> values(ele_vector.size());

    double min_val = std::numeric_limits<double>::max();
    double max_val = 0.0;

    os_mat << ele_vector.size() << std::endl;
    for (std::size_t i = 0; i < ele_vector.size(); i++)
    {
        auto const& element = *ele_vector[i];
        const auto mat_id = element.GetPatchIndex();
        auto& gaussian = *_guassian[mat_id];
        auto const& range = _ranges[mat_id];
        auto t = gaussian(gen);
        while (t < range[0] || t > range[1])
        {
            t = gaussian(gen);
        }

        min_val = std::min(t, min_val);
        max_val = std::max(t, max_val);
        values[i] = t;
        os_mat << i << delim << t << "\n";
    }

    os_vtu << "<DataArray type=\"Float64\" Name=\"" << _variable_name
           << "\" format=\"ascii\" RangeMin=\"" << min_val << "\" RangeMax=\""
           << max_val << "\">\n";

    os_vtu << values[0];
    for (std::size_t i = 1; i < ele_vector.size(); i++)
    {
        os_vtu << delim << values[i];
    }
    os_vtu << "\n</DataArray>\n";
}

}  // end of namespace MeshLib
