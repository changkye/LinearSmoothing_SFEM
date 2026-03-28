#pragma once

#include "nonlinear_esfem_t6.hpp"

#include <Eigen/Dense>

#include <array>
#include <vector>

namespace nonlinear_esfem_t6::shared
{

class T6MeshFactory
{
public:
    Mesh make_uniform(const Eigen::Matrix<double, 2, 2> &limits, const Eigen::Vector2i &num_els) const;
};

class PatchTestField
{
public:
    Eigen::MatrixXd exact_not_so_simple_shear(const Eigen::MatrixXd &nodes) const;
    Eigen::MatrixXd exact_bending_block(const Eigen::MatrixXd &nodes) const;
};

class BoundaryConditionBuilder
{
public:
    BoundaryCondition make_from_exact_field(const Eigen::MatrixXd &nodes,
                                            const Eigen::MatrixXd &exact_field) const;
};

class GeometryUtils
{
public:
    static std::vector<double> side_lengths(const Eigen::MatrixXd &triangle_vertices);
    static void outward_normals(const Eigen::MatrixXd &triangle_vertices,
                                const std::vector<double> &side,
                                std::vector<double> &nx,
                                std::vector<double> &ny);
    static double triangle_area(const Eigen::MatrixXd &triangle_vertices);
    static Eigen::MatrixXd element_coordinates(const Mesh &mesh, const std::array<int, 6> &element);
    static Eigen::MatrixXd nodal_displacements(const Eigen::VectorXd &uu);
};

} // namespace nonlinear_esfem_t6::shared
