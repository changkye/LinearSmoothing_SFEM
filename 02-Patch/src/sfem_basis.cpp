#include "sfem_basis.hpp"

#include <stdexcept>

namespace nonlinear_esfem_t6::shared
{

    QuadratureRule QuadratureLibrary::gauss_line(int order) const
    {
        QuadratureRule rule;
        if (order == 2)
        {
            rule.weights.resize(2);
            rule.points.resize(2, 1);
            rule.weights << 1.0, 1.0;
            rule.points << 0.5773502691896257, -0.5773502691896257;
            return rule;
        }
        if (order == 3)
        {
            rule.weights.resize(3);
            rule.points.resize(3, 1);
            rule.weights << 0.5555555555555556, 0.5555555555555556, 0.8888888888888888;
            rule.points << 0.7745966692414834, -0.7745966692414834, 0.0;
            return rule;
        }
        throw std::invalid_argument("Unsupported 1D quadrature order");
    }

    QuadratureRule QuadratureLibrary::gauss_quad(int order) const
    {
        QuadratureRule rule;
        if (order == 3)
        {
            const double pts[3] = {0.7745966692414834, -0.7745966692414834, 0.0};
            const double wts[3] = {0.5555555555555556, 0.5555555555555556, 0.8888888888888888};
            rule.weights.resize(9);
            rule.points.resize(9, 2);
            int k = 0;
            for (int i = 0; i < 3; ++i)
            {
                for (int j = 0; j < 3; ++j)
                {
                    rule.weights(k) = wts[i] * wts[j];
                    rule.points.row(k) << pts[i], pts[j];
                    ++k;
                }
            }
            return rule;
        }
        throw std::invalid_argument("Unsupported 2D quadrilateral quadrature order");
    }

    QuadratureRule QuadratureLibrary::triangle(int order) const
    {
        QuadratureRule rule;
        if (order == 2)
        {
            rule.weights = Eigen::Vector3d::Constant(1.0 / 6.0);
            rule.points.resize(3, 2);
            rule.points << 1.0 / 6.0, 1.0 / 6.0,
                2.0 / 3.0, 1.0 / 6.0,
                1.0 / 6.0, 2.0 / 3.0;
            return rule;
        }
        if (order == 3)
        {
            rule.weights.resize(4);
            rule.points.resize(4, 2);
            rule.points << 1.0 / 3.0, 1.0 / 3.0,
                0.2, 0.2,
                0.2, 0.6,
                0.6, 0.2;
            rule.weights << -0.28125, 0.2604166666666667, 0.2604166666666667,
                0.2604166666666667;
            return rule;
        }
        if (order == 4)
        {
            rule.weights.resize(6);
            rule.points.resize(6, 2);
            rule.points << 0.44594849091597, 0.44594849091597,
                0.44594849091597, 0.10810301816807,
                0.10810301816807, 0.44594849091597,
                0.09157621350977, 0.09157621350977,
                0.09157621350977, 0.81684757298046,
                0.81684757298046, 0.09157621350977;
            rule.weights << 0.111690794839005, 0.111690794839005, 0.111690794839005,
                0.05497587182766, 0.05497587182766, 0.05497587182766;
            return rule;
        }
        throw std::invalid_argument("Unsupported 2D triangular quadrature order");
    }

    void BasisEvaluator::lagrange(const std::string &type,
                                  const Eigen::Vector2d &coord,
                                  Eigen::VectorXd &n,
                                  Eigen::MatrixXd &dndxi) const
    {
        const double x = coord(0);
        const double y = coord(1);

        if (type == "L3")
        {
            n.resize(3);
            dndxi.resize(3, 1);
            n << -0.5 * (1.0 - x) * x, 0.5 * (1.0 + x) * x, 1.0 - x * x;
            dndxi << x - 0.5, x + 0.5, -2.0 * x;
            return;
        }

        if (type == "T6")
        {
            n.resize(6);
            dndxi.resize(6, 2);
            n << (1.0 - x - y) * (1.0 - 2.0 * x - 2.0 * y),
                x * (2.0 * x - 1.0),
                y * (2.0 * y - 1.0),
                4.0 * x * (1.0 - x - y),
                4.0 * x * y,
                4.0 * y * (1.0 - x - y);
            dndxi.col(0) << 4.0 * (x + y) - 3.0, 4.0 * x - 1.0, 0.0,
                -4.0 * (2.0 * x + y - 1.0), 4.0 * y, -4.0 * y;
            dndxi.col(1) << 4.0 * (x + y) - 3.0, 0.0, 4.0 * y - 1.0, -4.0 * x,
                4.0 * x, -4.0 * (x + 2.0 * y - 1.0);
            return;
        }

        if (type == "T3")
        {
            n.resize(3);
            dndxi.resize(3, 2);
            n << 1.0 - x - y, x, y;
            dndxi << -1.0, -1.0,
                1.0, 0.0,
                0.0, 1.0;
            return;
        }

        if (type == "Q4")
        {
            n.resize(4);
            dndxi.resize(4, 2);
            n << 0.25 * (1.0 - x) * (1.0 - y),
                0.25 * (1.0 + x) * (1.0 - y),
                0.25 * (1.0 + x) * (1.0 + y),
                0.25 * (1.0 - x) * (1.0 + y);
            dndxi << 0.25 * (y - 1.0), 0.25 * (x - 1.0),
                0.25 * (1.0 - y), -0.25 * (x + 1.0),
                0.25 * (y + 1.0), 0.25 * (x + 1.0),
                -0.25 * (y + 1.0), 0.25 * (1.0 - x);
            return;
        }

        if (type == "Q9")
        {
            n.resize(9);
            dndxi.resize(9, 2);
            n << 0.25 * x * y * (1.0 - x) * (1.0 - y),
                -0.25 * x * y * (1.0 + x) * (1.0 - y),
                0.25 * x * y * (1.0 + x) * (1.0 + y),
                -0.25 * x * y * (1.0 - x) * (1.0 + y),
                -0.5 * (1.0 - x * x) * (1.0 - y) * y,
                0.5 * (1.0 - y * y) * (1.0 + x) * x,
                0.5 * (1.0 - x * x) * (1.0 + y) * y,
                -0.5 * (1.0 - y * y) * (1.0 - x) * x,
                (1.0 - x * x) * (1.0 - y * y);
            dndxi.col(0) << 0.25 * y * (2.0 * x - 1.0) * (y - 1.0),
                0.25 * y * (2.0 * x + 1.0) * (y - 1.0),
                0.25 * y * (2.0 * x + 1.0) * (y + 1.0),
                0.25 * y * (2.0 * x - 1.0) * (y + 1.0),
                -x * (y - 1.0) * y,
                -0.5 * (2.0 * x + 1.0) * (y * y - 1.0),
                -x * (y + 1.0) * y,
                -0.5 * (2.0 * x - 1.0) * (y * y - 1.0),
                2.0 * x * (y * y - 1.0);
            dndxi.col(1) << 0.25 * x * (x - 1.0) * (2.0 * y - 1.0),
                0.25 * x * (x + 1.0) * (2.0 * y - 1.0),
                0.25 * x * (x + 1.0) * (2.0 * y + 1.0),
                0.25 * x * (x - 1.0) * (2.0 * y + 1.0),
                -0.5 * (x * x - 1.0) * (2.0 * y - 1.0),
                -x * (x + 1.0) * y,
                -0.5 * (x * x - 1.0) * (2.0 * y + 1.0),
                -x * (x - 1.0) * y,
                2.0 * (x * x - 1.0) * y;
            return;
        }

        throw std::invalid_argument("Unsupported basis type");
    }

    Eigen::Vector2d BasisEvaluator::inverse_mapping(const std::string &type,
                                                    const Eigen::MatrixXd &nodes,
                                                    const Eigen::Vector2d &point) const
    {
        Eigen::Vector2d xi_eta = Eigen::Vector2d::Zero();
        for (int iter = 0; iter < 10; ++iter)
        {
            Eigen::VectorXd n;
            Eigen::MatrixXd dndxi;
            lagrange(type, xi_eta, n, dndxi);
            const Eigen::Vector2d mapped = nodes.transpose() * n;
            const Eigen::Matrix2d jac = nodes.transpose() * dndxi;
            xi_eta -= jac.inverse() * (mapped - point);
        }
        return xi_eta;
    }

    Eigen::VectorXd BasisEvaluator::serendipity_shape(const std::string &type,
                                                      const Eigen::MatrixXd &nodes,
                                                      const Eigen::Vector2d &point) const
    {
        const Eigen::Vector2d xi_eta = inverse_mapping(type, nodes, point);
        Eigen::VectorXd n;
        Eigen::MatrixXd dndxi;
        lagrange(type, xi_eta, n, dndxi);
        return n;
    }

} // namespace nonlinear_esfem_t6::shared
