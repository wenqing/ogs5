/*! \file ShapeFunctionPool.cpp
    \brief Compute shape functions and their gradients with respect to
     the local coodinates, and store the results

     \author Wenqing Wang
     \date Feb. 2015

     \copyright
      Copyright (c) 2018, OpenGeoSys Community (http://www.opengeosys.org)
             Distributed under a Modified BSD License.
             See accompanying file LICENSE.txt or
             http://www.opengeosys.org/project/license
*/

#include "ShapeFunctionPool.h"

#include <cassert> /* assert */

#include "display.h"
#include "fem_ele.h"
#include "matrix_class.h"

#include "LinAlg/GaussAlgorithm.h"

namespace FiniteElement
{
ShapeFunctionPool::ShapeFunctionPool(
    const std::vector<MshElemType::type>& elem_types, CElement& quadrature,
    const int num_sample_gs_pnts)
{
    int num_elem_nodes[2][MshElemType::NUM_ELEM_TYPES];
    int dim_elem[MshElemType::NUM_ELEM_TYPES];

    int id = static_cast<int>(MshElemType::LINE) - 1;
    num_elem_nodes[0][id] = 2;
    num_elem_nodes[1][id] = 3;
    dim_elem[id] = 1;

    id = static_cast<int>(MshElemType::QUAD) - 1;
    num_elem_nodes[0][id] = 4;
    num_elem_nodes[1][id] = 9;
    dim_elem[id] = 2;

    id = static_cast<int>(MshElemType::QUAD8) - 1;
    num_elem_nodes[0][id] = 4;
    num_elem_nodes[1][id] = 8;
    dim_elem[id] = 2;

    id = static_cast<int>(MshElemType::TRIANGLE) - 1;
    num_elem_nodes[0][id] = 3;
    num_elem_nodes[1][id] = 6;
    dim_elem[id] = 2;

    id = static_cast<int>(MshElemType::HEXAHEDRON) - 1;
    num_elem_nodes[0][id] = 8;
    num_elem_nodes[1][id] = 20;
    dim_elem[id] = 3;

    id = static_cast<int>(MshElemType::TETRAHEDRON) - 1;
    num_elem_nodes[0][id] = 4;
    num_elem_nodes[1][id] = 10;
    dim_elem[id] = 3;

    id = static_cast<int>(MshElemType::PRISM) - 1;
    num_elem_nodes[0][id] = 6;
    num_elem_nodes[1][id] = 15;
    dim_elem[id] = 3;

    id = static_cast<int>(MshElemType::PYRAMID) - 1;
    num_elem_nodes[0][id] = 5;
    num_elem_nodes[1][id] = 13;
    dim_elem[id] = 3;

    // std::vector<int> elem_type_ids
    const std::size_t n_ele_types = elem_types.size();
    _shape_function.resize(n_ele_types);
    _shape_function_center.resize(n_ele_types);
    _grad_shape_function.resize(n_ele_types);
    _grad_shape_function_center.resize(n_ele_types);
    for (std::size_t i = 0; i < n_ele_types; i++)
    {
        const MshElemType::type e_type = elem_types[i];
        if (e_type == MshElemType::INVALID)
            continue;

        // Set number of integration points.
        quadrature.SetGaussPointNumber(num_sample_gs_pnts);
        quadrature.SetIntegrationPointNumber(e_type);

        const int type_id = static_cast<int>(e_type) - 1;
        int num_int_pnts = quadrature.GetNumGaussPoints();
        const int num_nodes =
            num_elem_nodes[quadrature.getOrder() - 1][type_id];

        std::vector<double> elem_shape_function_center(num_nodes);
        _shape_function_center[type_id] = elem_shape_function_center;

        const int size_shape_fct = num_nodes * num_int_pnts;
        std::vector<double> elem_shape_function(size_shape_fct);
        _shape_function[type_id] = elem_shape_function;

        std::vector<double> elem_grad_shape_function(dim_elem[type_id] *
                                                     size_shape_fct);
        _grad_shape_function[type_id] = elem_grad_shape_function;

        std::vector<double> elem_grad_shape_function_center(dim_elem[type_id] *
                                                            num_nodes);
        _grad_shape_function_center[type_id] = elem_grad_shape_function_center;
    }

    computeQuadratures(elem_types, num_elem_nodes, dim_elem, quadrature,
                       num_sample_gs_pnts);
    computeInverseExtrapolationMatrices(elem_types, num_elem_nodes, quadrature,
                       num_sample_gs_pnts);
}

ShapeFunctionPool::~ShapeFunctionPool()
{
    for (std::size_t i = 0; i < _inverse_extrapolation_matrices.size(); i++)
    {
        if (_inverse_extrapolation_matrices[i])
            delete _inverse_extrapolation_matrices[i];
    }
}

void ShapeFunctionPool::computeQuadratures(
    const std::vector<MshElemType::type>& elem_types,
    const int num_elem_nodes[2][MshElemType::NUM_ELEM_TYPES],
    const int dim_elem[], CElement& quadrature, const int num_sample_gs_pnts)
{
    const int order = quadrature.getOrder();
    for (std::size_t i = 0; i < elem_types.size(); i++)
    {
        const MshElemType::type e_type = elem_types[i];
        if (e_type == MshElemType::INVALID)
            continue;

        const int type_id = static_cast<int>(e_type) - 1;
        quadrature.ConfigShapefunction(e_type);

        const int nnodes = num_elem_nodes[order - 1][type_id];
        const int elem_dim = dim_elem[type_id];

        double* shape_function_center_values =
            _shape_function_center[type_id].data();
        quadrature.SetCenterGP(e_type);
        quadrature.ComputeShapefct(order, shape_function_center_values);
        double* grad_shape_function_center_values =
            _grad_shape_function_center[type_id].data();
        quadrature.computeGradShapefctLocal(order,
                                            grad_shape_function_center_values);

        double* shape_function_values = _shape_function[type_id].data();
        double* dshape_function_values = _grad_shape_function[type_id].data();
        // Set number of integration points.
        quadrature.SetGaussPointNumber(num_sample_gs_pnts);
        quadrature.SetIntegrationPointNumber(e_type);

        for (int gp = 0; gp < quadrature.GetNumGaussPoints(); gp++)
        {
            int gp_r, gp_s, gp_t;
            quadrature.SetGaussPoint(e_type, gp, gp_r, gp_s, gp_t);
            double* shape_function_values_gs =
                &shape_function_values[gp * nnodes];
            quadrature.ComputeShapefct(order, shape_function_values_gs);

            double* dshape_function_values_gs =
                &dshape_function_values[gp * nnodes * elem_dim];
            quadrature.computeGradShapefctLocal(order,
                                                dshape_function_values_gs);
        }
    }
}

void ShapeFunctionPool::computeInverseExtrapolationMatrices(
    const std::vector<MshElemType::type>& elem_types,
    const int num_elem_nodes[2][MshElemType::NUM_ELEM_TYPES],
    CElement& quadrature, const int num_sample_gs_pnts)
{
    _inverse_extrapolation_matrices.resize(elem_types.size());
    for (std::size_t i = 0; i < elem_types.size(); i++)
        _inverse_extrapolation_matrices[i] = NULL;

    for (std::size_t i = 0; i < elem_types.size(); i++)
    {
        const MshElemType::type e_type = elem_types[i];
        const int type_id = static_cast<int>(e_type) - 1;
        if (e_type == MshElemType::INVALID)
        {
            continue;
        }

        const int nnodes = num_elem_nodes[0][type_id];

        double const* const shape_function_at_ips =
            getShapeFunctionValues(e_type);

        // Set number of integration points.
        quadrature.SetGaussPointNumber(num_sample_gs_pnts);
        quadrature.SetIntegrationPointNumber(e_type);

        if (quadrature.GetNumGaussPoints() < nnodes)
        {
            Display::ScreenMessage(
                "Number of integrations points is smaller than the number of "
                "element nodes");
            abort();
        }

        Math_Group::Matrix A(nnodes, nnodes);

        // Loop over integration points to a limited number of nnodes
        for (int i = 0; i < nnodes; i++)
        {
            double const* const shapefct_at_ip =
                &shape_function_at_ips[nnodes * i];
            for (int j = 0; j < nnodes; j++)
                A(i, j) = shapefct_at_ip[j];
        }

        MathLib::GaussAlgorithm<Math_Group::Matrix> linear_solver(
            A, static_cast<std::size_t>(nnodes));

        Math_Group::Matrix* inverse_extrapolation_matrix =
            new Math_Group::Matrix(nnodes, nnodes);

        _inverse_extrapolation_matrices[type_id] = inverse_extrapolation_matrix;

        double b[20];
        for (int i = 0; i < nnodes; i++)
        {
            b[i] = 0.0;
        }
        b[0] = 1.0;

        // Compute first column, and get the LU decomposition
        linear_solver.execute(b);
        for (int i = 0; i < nnodes; i++)
            (*inverse_extrapolation_matrix)(i, 0) = b[i];

        for (int j = 1; j < nnodes; j++)
        {
            for (int i = 0; i < nnodes; i++)
            {
                b[i] = (i == j) ? 1.0 : 0.0;
            }
            linear_solver.executeWithExistedElimination(b);
            for (int i = 0; i < nnodes; i++)
                (*inverse_extrapolation_matrix)(i, j) = b[i];
        }
    }
}

const double* ShapeFunctionPool::getShapeFunctionValues(
    const MshElemType::type elem_type) const
{
    return _shape_function[static_cast<int>(elem_type) - 1].data();
}

unsigned ShapeFunctionPool::getShapeFunctionArraySize(
    const MshElemType::type elem_type) const
{
    return _shape_function[static_cast<int>(elem_type) - 1].size();
}

const double* ShapeFunctionPool::getShapeFunctionCenterValues(
    const MshElemType::type elem_type) const
{
    return _shape_function_center[static_cast<int>(elem_type) - 1].data();
}

const double* ShapeFunctionPool::getGradShapeFunctionValues(
    const MshElemType::type elem_type) const
{
    return _grad_shape_function[static_cast<int>(elem_type) - 1].data();
}

const double* ShapeFunctionPool::getGradShapeFunctionCenterValues(
    const MshElemType::type elem_type) const
{
    return _grad_shape_function_center[static_cast<int>(elem_type) - 1].data();
}

const Math_Group::Matrix* ShapeFunctionPool::getInverseExtrapolationMatrix(
    const MshElemType::type elem_type) const
{
    return _inverse_extrapolation_matrices[static_cast<int>(elem_type) - 1];
}

}  // namespace FiniteElement
