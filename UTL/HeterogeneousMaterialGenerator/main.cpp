/**
 * \copyright
 * Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * File:   main.cpp
 *
 * Created on March 18, 2019, 3:09 PM
 */

#include <cmath>
#include <fstream>
#include <string>
#include <vector>

#include "HeterogeneousMaterialGenerator.h"

#include "Base/FileTools.h"

void DisplayOption()
{
    std::string opt =
        "Options:\n"
        "    --version:  display version number.\n"
        "    --help:     display help info.\n"
        "    --o: set output directory.\n";
    std::cout << opt << std::endl;
}
void DisplayMessage()
{
    std::string info =
        "It generates element wise heterogeneous material parameters for\n "
        "task 4, DECOVALEX 2019 by using the Gaussian distribution\n "
        "Command: hete_mat [mesh file name] [material file name]\n";
    std::cout << info << std::endl;
}

void DisplayMessageConsole()
{
    DisplayMessage();

    std::cout << "Run in console.\nInput file name or with option:"
              << std::endl;
}
void DisplayVersion()
{
    std::string ver = "Version: 1.0.";
    std::cout << ver << std::endl;
}

int main(int argc, char* argv[])
{
    std::vector<std::string> arg_strings;
    if (argc > 1)
    {
        for (int i = 1; i < argc; i++)
        {
            arg_strings.push_back(std::string(argv[i]));
        }
    }
    else
    {
        DisplayMessageConsole();

        std::string s_buff;
        std::stringstream ss;
        getline(std::cin, s_buff);
        ss.str(s_buff);
        while (!ss.eof())
        {
            ss >> s_buff;
            arg_strings.push_back(s_buff);
        }
    }

    std::string file_name_msh;
    std::string file_name_mat;
    std::string o_path;
    for (std::size_t i = 0; i < arg_strings.size(); i++)
    {
        const std::string anArg = arg_strings[i];
        if (anArg == "--help" || anArg == "-h")
        {
            DisplayMessage();
            DisplayOption();
            exit(0);
        }
        if (anArg == "--version")
        {
            DisplayVersion();
            exit(0);
        }
        if (anArg == "--o")
        {
            if (i + 1 >= arg_strings.size())
            {
                std::cerr << "Error: Parameter " << anArg
                          << " needs a path for output files" << std::endl;
                std::exit(EXIT_FAILURE);
            }
            std::string path = arg_strings[++i];

            if (!path.empty())
                o_path = path;
            continue;
        }
        else
        {
            if (file_name_msh.empty() && file_name_msh.empty())
            {
                file_name_msh = arg_strings[i];
                file_name_mat = arg_strings[i + 1];
            }
        }
    }

    if (file_name_msh.empty())
    {
        std::cout << "The name of mesh file is empty.";
        return EXIT_FAILURE;
    }
    if (file_name_mat.empty())
    {
        std::cout << "The name of material file is empty.";
        return EXIT_FAILURE;
    }

    //
    const std::string file_path = pathDirname(file_name_msh);
    const std::string base_file_name_msh0 =
        pathJoin(file_path, pathBasename(file_name_msh));

    const auto pos_msh = base_file_name_msh0.find_last_of('.');
    const std::string base_file_name_msh =
        pos_msh > 0 ? base_file_name_msh0.substr(0, pos_msh)
                    : base_file_name_msh0;

    const std::string base_file_name_mat =
        pathJoin(file_path, pathBasename(file_name_mat));

    if (o_path.empty())
        o_path = file_path;

    std::string fname_msh = base_file_name_msh + ".msh";
    std::ifstream is_mesh(fname_msh.c_str(), std::ios::in);

    std::string s_buff;
    std::getline(is_mesh, s_buff);

    if (s_buff.find("#FEM_MSH") == std::string::npos)
    {
        std::cout << "Cannot open mesh file " << fname_msh << std::endl;
        return EXIT_FAILURE;
    }

    auto hete_mat_generator =
        MeshLib::createHeterogeneousMaterialGenerator(base_file_name_mat);
    hete_mat_generator->Read(&is_mesh);

    hete_mat_generator->generate(pathJoin(o_path, base_file_name_mat));

    std::cout << "Terminate normally ^O^ ." << std::endl;
    return EXIT_SUCCESS;
}
