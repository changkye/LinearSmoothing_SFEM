#include "sfem_mesh.hpp"

#include <cmath>

namespace nonlinear_esfem_t6::shared
{

Mesh T6MeshFactory::make_uniform(const Eigen::Matrix<double, 2, 2> &limits, const Eigen::Vector2i &num_els) const
{
    const int xn = 2 * num_els(0) + 1;
    const int yn = 2 * num_els(1) + 1;

    Mesh mesh;
    mesh.nodes.resize(xn * yn, 2);

    int index = 0;
    for (int j = 0; j < yn; ++j)
    {
        const double y = limits(1, 0) + (limits(1, 1) - limits(1, 0)) * static_cast<double>(j) / static_cast<double>(yn - 1);
        for (int i = 0; i < xn; ++i)
        {
            const double x = limits(0, 0) + (limits(0, 1) - limits(0, 0)) * static_cast<double>(i) / static_cast<double>(xn - 1);
            mesh.nodes.row(index++) << x, y;
        }
    }

    const std::array<int, 6> idx1 = {0, 2, 2 * xn, 1, xn + 1, xn};
    const std::array<int, 6> idx2 = {2, 2 * xn + 2, 2 * xn, xn + 2, 2 * xn + 1, xn + 1};

    auto add_family = [&](const std::array<int, 6> &base) {
        int inc = 0;
        for (int row = 0; row < num_els(1); ++row)
        {
            for (int col = 0; col < num_els(0); ++col)
            {
                std::array<int, 6> element{};
                for (int a = 0; a < 6; ++a)
                {
                    element[static_cast<std::size_t>(a)] = base[static_cast<std::size_t>(a)] + inc;
                }
                mesh.elements.push_back(element);
                inc += 2;
            }
            inc = (row + 1) * (2 * xn);
        }
    };

    add_family(idx1);
    add_family(idx2);
    return mesh;
}

Eigen::MatrixXd PatchTestField::exact_not_so_simple_shear(const Eigen::MatrixXd &nodes) const
{
    Eigen::MatrixXd exact = Eigen::MatrixXd::Zero(nodes.rows(), 2);
    for (int i = 0; i < nodes.rows(); ++i)
    {
        const double y = nodes(i, 1);
        exact(i, 0) = 0.5 * y * y;
        exact(i, 1) = 0.0;
    }
    return exact;
}

Eigen::MatrixXd PatchTestField::exact_bending_block(const Eigen::MatrixXd &nodes) const
{
    constexpr double alpha = 0.9;
    Eigen::MatrixXd exact = Eigen::MatrixXd::Zero(nodes.rows(), 2);
    for (int i = 0; i < nodes.rows(); ++i)
    {
        const double x = nodes(i, 0);
        const double y = nodes(i, 1);
        exact(i, 0) = std::sqrt(2.0 * alpha * x) * std::cos(y / alpha) - x;
        exact(i, 1) = std::sqrt(2.0 * alpha * x) * std::sin(y / alpha) - y;
    }
    return exact;
}

BoundaryCondition BoundaryConditionBuilder::make_from_exact_field(const Eigen::MatrixXd &nodes,
                                                                  const Eigen::MatrixXd &exact) const
{
    const double xmin = nodes.col(0).minCoeff();
    const double xmax = nodes.col(0).maxCoeff();
    const double ymin = nodes.col(1).minCoeff();
    const double ymax = nodes.col(1).maxCoeff();

    std::vector<int> fixed_nodes;
    for (int i = 0; i < nodes.rows(); ++i)
    {
        const double x = nodes(i, 0);
        const double y = nodes(i, 1);
        if (std::abs(x - xmin) < 1.0e-12 || std::abs(x - xmax) < 1.0e-12 ||
            std::abs(y - ymin) < 1.0e-12 || std::abs(y - ymax) < 1.0e-12)
        {
            fixed_nodes.push_back(i);
        }
    }
    BoundaryCondition bc;
    bc.values.resize(static_cast<int>(fixed_nodes.size() * 2));
    bc.dofs.reserve(fixed_nodes.size() * 2);

    int cursor = 0;
    for (const int node : fixed_nodes)
    {
        bc.dofs.push_back(2 * node);
        bc.dofs.push_back(2 * node + 1);
        bc.values(cursor++) = exact(node, 0);
        bc.values(cursor++) = exact(node, 1);
    }
    return bc;
}

std::vector<double> GeometryUtils::side_lengths(const Eigen::MatrixXd &triangle_vertices)
{
    std::vector<double> side(3);
    for (int i = 0; i < 3; ++i)
    {
        const int j = (i + 1) % 3;
        side[static_cast<std::size_t>(i)] = (triangle_vertices.row(j) - triangle_vertices.row(i)).norm();
    }
    return side;
}

void GeometryUtils::outward_normals(const Eigen::MatrixXd &triangle_vertices,
                                    const std::vector<double> &side,
                                    std::vector<double> &nx,
                                    std::vector<double> &ny)
{
    nx.resize(3);
    ny.resize(3);
    for (int i = 0; i < 3; ++i)
    {
        const int j = (i + 1) % 3;
        nx[static_cast<std::size_t>(i)] = (triangle_vertices(j, 1) - triangle_vertices(i, 1)) / side[static_cast<std::size_t>(i)];
        ny[static_cast<std::size_t>(i)] = -(triangle_vertices(j, 0) - triangle_vertices(i, 0)) / side[static_cast<std::size_t>(i)];
    }
}

double GeometryUtils::triangle_area(const Eigen::MatrixXd &triangle_vertices)
{
    return 0.5 * std::abs((triangle_vertices(1, 0) - triangle_vertices(0, 0)) * (triangle_vertices(2, 1) - triangle_vertices(0, 1)) -
                          (triangle_vertices(2, 0) - triangle_vertices(0, 0)) * (triangle_vertices(1, 1) - triangle_vertices(0, 1)));
}

Eigen::MatrixXd GeometryUtils::element_coordinates(const Mesh &mesh, const std::array<int, 6> &element)
{
    Eigen::MatrixXd coords(6, 2);
    for (int a = 0; a < 6; ++a)
    {
        coords.row(a) = mesh.nodes.row(element[static_cast<std::size_t>(a)]);
    }
    return coords;
}

Eigen::MatrixXd GeometryUtils::nodal_displacements(const Eigen::VectorXd &uu)
{
    Eigen::MatrixXd u(uu.size() / 2, 2);
    for (int i = 0; i < u.rows(); ++i)
    {
        u(i, 0) = uu(2 * i);
        u(i, 1) = uu(2 * i + 1);
    }
    return u;
}

} // namespace nonlinear_esfem_t6::shared
