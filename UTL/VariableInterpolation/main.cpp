/**
 * \file
 *
 * \copyright
 * Copyright (c) 2012-2019, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 *
 * Created on October 28, 2019, 1:50 PM
 */

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

#include "Base/FileTools.h"

#include "CreateVariableValues.h"
#include "ShapeFunctionPool.h"
#include "VariableValues.h"

using namespace std;

void DisplayMessage()
{
    static std::string info =
        "*  Extract point data from pvd and vtu in order to generate curves "
        "over time.\n"
        "*  Syntax of input file:\n"
        "*    foo_mesh.msh   // name of pvd file\n"
        "*    foo_pvd.pvd    // name of pvd file\n"
        "*    3         // number of specified points\n"
        "*    2  0.1 0.2 0.1 $NAME foo_point1  // ID x y z $NAME file_name\n"
        "*    10 0.1 0.2 0.1 $NAME foo_point2\n"
        "*    12 0.1 0.2 0.1 $NAME foo_point3\n";
    std::cout << info << std::endl;
}

void DisplayMessageConsole()
{
    DisplayMessage();

    std::cout << "Run in console.\nInput file name:" << std::endl;
}

void DisplayOption()
{
    std::string opt =
        "Options:\n"
        "    --help:     display help info.\n"
        "    -o: set output directory.\n";
    std::cout << opt << std::endl;
}

/*
 *  Extract point data from pvd and vtu in order to generate curves over time.
 *  Syntax of input file:
 *    foo.pvd   // name of pvd file
 *    3         // number of specified points
 *    2  0.1 0.2 0.1 $NAME foo_point1  // ID x y z $NAME file_name
 *    10 0.1 0.2 0.1 $NAME foo_point2
 *    12 0.1 0.2 0.1 $NAME foo_point3
 *
 */
int main(int argc, char** argv)
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

    std::string file_name;
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
        if (anArg == "-o")
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
            file_name = arg_strings[i];
        }
    }

    if (argc > 1)
    {
        DisplayMessage();
    }

    const std::string file_path = pathDirname(file_name);
    if (o_path.empty())
        o_path = file_path;

    FiniteElement::ShapeFunctionPool* linear_shapefunction_pool = NULL;
    UTL::VariableValues* variable_values = UTL::createVariableValues(
        file_path, file_name, linear_shapefunction_pool);

    delete linear_shapefunction_pool;
    return 0;
}
