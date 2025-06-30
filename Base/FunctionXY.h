/**
 * \file
 * \copyright
 * Copyright (c) 2012-2025, OpenGeoSys Community (http://www.opengeosys.org)
 *            Distributed under a Modified BSD License.
 *              See accompanying file LICENSE.txt or
 *              http://www.opengeosys.org/project/license
 *
 * Created on 2025-06-30 13:36:27
 */

#pragma once

#include <cmath>
#include <exprtk.hpp>
#include <string>
#include <vector>

namespace BaseLib
{
class FunctionXY
{
public:
    using T = double;  // Using double precision by default
    using symbol_table_t = exprtk::symbol_table<T>;
    using expression_t = exprtk::expression<T>;
    using parser_t = exprtk::parser<T>;

    FunctionXY(const std::string& expression_str,
               const std::string& variable_x = "x",
               const std::string& variable_y = "y")
        : _variable_x(variable_x), _variable_y(variable_y)
    {
        // Setup symbol table and bind references directly
        _symbol_table = exprtk::symbol_table<T>();

        // This is the correct way to get references to variables
        _symbol_table.create_variable(_variable_x);
        _symbol_table.create_variable(_variable_y);

        // Add constants if needed (pi, epsilon, etc.)
        _symbol_table.add_constants();

        // Setup expression
        _expression = exprtk::expression<T>();
        _expression.register_symbol_table(_symbol_table);

        // Compile expression
        exprtk::parser<T> parser;
        if (!parser.compile(expression_str, _expression))
        {
            throw std::runtime_error("Compilation error: " + parser.error());
        }
    }

    // Single evaluation (function-style syntax)
    T operator()(T const x_val, T const y_val) const
    {
        auto& x = _symbol_table.get_variable(_variable_x)->ref();
        auto& y = _symbol_table.get_variable(_variable_y)->ref();
        x = x_val;
        y = y_val;
        return _expression.value();
    }

private:
    std::string const _variable_x;
    std::string const _variable_y;
    symbol_table_t _symbol_table;
    expression_t _expression;
};

}  // namespace BaseLib
