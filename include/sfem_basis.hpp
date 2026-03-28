#pragma once

#include <Eigen/Dense>

#include <string>

namespace nonlinear_esfem_t6::shared
{

struct QuadratureRule
{
    Eigen::VectorXd weights;
    Eigen::MatrixXd points;
};

class QuadratureLibrary
{
public:
    QuadratureRule gauss_line(int order) const;
    QuadratureRule gauss_quad(int order) const;
    QuadratureRule triangle(int order) const;
};

class BasisEvaluator
{
public:
    void lagrange(const std::string &type,
                  const Eigen::Vector2d &coord,
                  Eigen::VectorXd &n,
                  Eigen::MatrixXd &dndxi) const;
    Eigen::Vector2d inverse_mapping(const std::string &type,
                                    const Eigen::MatrixXd &nodes,
                                    const Eigen::Vector2d &point) const;
    Eigen::VectorXd serendipity_shape(const std::string &type,
                                      const Eigen::MatrixXd &nodes,
                                      const Eigen::Vector2d &point) const;
};

} // namespace nonlinear_esfem_t6::shared
