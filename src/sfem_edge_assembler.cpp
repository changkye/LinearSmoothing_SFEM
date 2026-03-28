#include "sfem_assemblers.hpp"
#include "sfem_basis.hpp"
#include "sfem_mesh.hpp"

#include <algorithm>

namespace nonlinear_esfem_t6::shared
{
namespace
{

struct TargetEdge
{
    int nodes[2] = {0, 0};
    int elem1 = 0;
    int elem2 = -1;
    int local_edge1 = 0;
    int local_edge2 = -1;
};

} // namespace

std::vector<SupportDomain> EdgeSmoothedDomainAssembler::build_domains(const Problem &problem) const
{
    QuadratureLibrary quadrature;
    BasisEvaluator basis;

    std::vector<std::vector<int>> support_nodes;
    support_nodes.reserve(problem.mesh.elements.size());
    for (const auto &element : problem.mesh.elements)
    {
        support_nodes.emplace_back(element.begin(), element.end());
    }

    std::vector<TargetEdge> edges;
    edges.reserve(problem.mesh.elements.size() * 2);
    for (int ie = 0; ie < static_cast<int>(problem.mesh.elements.size()); ++ie)
    {
        const auto &element = problem.mesh.elements[static_cast<std::size_t>(ie)];
        for (int j = 0; j < 3; ++j)
        {
            const int n1 = element[static_cast<std::size_t>(j)];
            const int n2 = element[static_cast<std::size_t>((j + 1) % 3)];
            bool found = false;
            for (auto &edge : edges)
            {
                if ((edge.nodes[0] == n1 && edge.nodes[1] == n2) || (edge.nodes[0] == n2 && edge.nodes[1] == n1))
                {
                    edge.elem2 = ie;
                    edge.local_edge2 = j;
                    found = true;
                    break;
                }
            }
            if (!found)
            {
                TargetEdge edge;
                edge.nodes[0] = n1;
                edge.nodes[1] = n2;
                edge.elem1 = ie;
                edge.local_edge1 = j;
                edges.push_back(edge);
            }
        }
    }

    std::vector<double> tri_area(problem.mesh.elements.size(), 0.0);
    for (int ie = 0; ie < static_cast<int>(problem.mesh.elements.size()); ++ie)
    {
        tri_area[static_cast<std::size_t>(ie)] =
            GeometryUtils::triangle_area(GeometryUtils::element_coordinates(problem.mesh, problem.mesh.elements[static_cast<std::size_t>(ie)]).topRows(3));
    }

    std::vector<double> sub_area(edges.size(), 0.0);
    for (int i = 0; i < static_cast<int>(edges.size()); ++i)
    {
        sub_area[static_cast<std::size_t>(i)] = tri_area[static_cast<std::size_t>(edges[static_cast<std::size_t>(i)].elem1)] / 3.0;
        if (edges[static_cast<std::size_t>(i)].elem2 >= 0)
        {
            sub_area[static_cast<std::size_t>(i)] += tri_area[static_cast<std::size_t>(edges[static_cast<std::size_t>(i)].elem2)] / 3.0;
        }
    }

    std::vector<SupportDomain> domains;
    domains.reserve(edges.size());
    const int bound[3][3] = {{0, 1, 3}, {1, 2, 4}, {2, 0, 5}};
    const QuadratureRule boundary_rule = quadrature.gauss_line(problem.edge_boundary_ng);

    for (int edge_index = 0; edge_index < static_cast<int>(edges.size()); ++edge_index)
    {
        const TargetEdge &target = edges[static_cast<std::size_t>(edge_index)];
        std::vector<int> neighbour = {target.elem1};
        if (target.elem2 >= 0)
        {
            neighbour.push_back(target.elem2);
        }

        Eigen::MatrixXd fx;
        Eigen::MatrixXd fy;
        std::vector<int> nodl;
        for (int ic = 0; ic < static_cast<int>(neighbour.size()); ++ic)
        {
            const auto &element = problem.mesh.elements[static_cast<std::size_t>(neighbour[static_cast<std::size_t>(ic)])];
            const Eigen::MatrixXd wkx = GeometryUtils::element_coordinates(problem.mesh, element);
            std::vector<double> nx;
            std::vector<double> ny;
            GeometryUtils::outward_normals(wkx.topRows(3), GeometryUtils::side_lengths(wkx.topRows(3)), nx, ny);

            Eigen::MatrixXd local_fx = Eigen::MatrixXd::Zero(6, 6);
            Eigen::MatrixXd local_fy = Eigen::MatrixXd::Zero(6, 6);
            for (int is = 0; is < 3; ++is)
            {
                Eigen::MatrixXd bxy = Eigen::MatrixXd::Zero(6, 6);
                Eigen::MatrixXd edge_coords(3, 2);
                for (int a = 0; a < 3; ++a)
                {
                    edge_coords.row(a) = wkx.row(bound[is][a]);
                }

                for (int ig = 0; ig < boundary_rule.weights.size(); ++ig)
                {
                    Eigen::VectorXd ng;
                    Eigen::MatrixXd dndxi;
                    basis.lagrange("L3", Eigen::Vector2d(boundary_rule.points(ig, 0), 0.0), ng, dndxi);
                    const Eigen::Vector2d gp = edge_coords.transpose() * ng;
                    const double detj = (edge_coords.transpose() * dndxi).norm();
                    const Eigen::VectorXd nt = basis.serendipity_shape("T6", wkx, gp);
                    Eigen::VectorXd p(6);
                    p << 1.0, gp(0), gp(1), gp(0) * gp(0), gp(0) * gp(1), gp(1) * gp(1);
                    bxy += p * (nt.transpose() * detj * boundary_rule.weights(ig));
                }
                local_fx += nx[static_cast<std::size_t>(is)] * bxy;
                local_fy += ny[static_cast<std::size_t>(is)] * bxy;
            }

            if (ic == 0)
            {
                nodl = support_nodes[static_cast<std::size_t>(neighbour[static_cast<std::size_t>(ic)])];
                fx = local_fx;
                fy = local_fy;
            }
            else
            {
                for (int local = 0; local < static_cast<int>(support_nodes[static_cast<std::size_t>(neighbour[static_cast<std::size_t>(ic)])].size()); ++local)
                {
                    const int node = support_nodes[static_cast<std::size_t>(neighbour[static_cast<std::size_t>(ic)])][static_cast<std::size_t>(local)];
                    auto pos = std::find(nodl.begin(), nodl.end(), node);
                    if (pos != nodl.end())
                    {
                        const int index = static_cast<int>(std::distance(nodl.begin(), pos));
                        fx.col(index) += local_fx.col(local);
                        fy.col(index) += local_fy.col(local);
                    }
                    else
                    {
                        nodl.push_back(node);
                        fx.conservativeResize(Eigen::NoChange, fx.cols() + 1);
                        fy.conservativeResize(Eigen::NoChange, fy.cols() + 1);
                        fx.col(fx.cols() - 1) = local_fx.col(local);
                        fy.col(fy.cols() - 1) = local_fy.col(local);
                    }
                }
            }
        }

        std::vector<int> reorder;
        std::string parent_element;
        QuadratureRule interior_rule;
        if (target.elem2 < 0)
        {
            reorder = {0, 1, 2, 3, 4, 5};
            parent_element = "T6";
            interior_rule = quadrature.triangle(problem.edge_tri_quad_order);
        }
        else
        {
            if (target.local_edge1 == 0)
            {
                reorder = {0, 6, 1, 2, 8, 7, 4, 5, 3};
            }
            else if (target.local_edge1 == 1)
            {
                reorder = {0, 1, 6, 2, 3, 7, 8, 5, 4};
            }
            else
            {
                reorder = {0, 1, 2, 6, 3, 4, 7, 8, 5};
            }
            parent_element = "Q9";
            interior_rule = quadrature.gauss_quad(problem.edge_quad_quad_order);
        }

        std::vector<int> reordered_nodes;
        reordered_nodes.reserve(reorder.size());
        for (const int index : reorder)
        {
            reordered_nodes.push_back(nodl[static_cast<std::size_t>(index)]);
        }

        Eigen::MatrixXd gcoord(reordered_nodes.size(), 2);
        Eigen::MatrixXd reordered_fx(6, reordered_nodes.size());
        Eigen::MatrixXd reordered_fy(6, reordered_nodes.size());
        for (int i = 0; i < static_cast<int>(reordered_nodes.size()); ++i)
        {
            gcoord.row(i) = problem.mesh.nodes.row(reordered_nodes[static_cast<std::size_t>(i)]);
            reordered_fx.col(i) = fx.col(reorder[static_cast<std::size_t>(i)]);
            reordered_fy.col(i) = fy.col(reorder[static_cast<std::size_t>(i)]);
        }

        Eigen::MatrixXd m = Eigen::MatrixXd::Zero(6, 6);
        Eigen::MatrixXd p(interior_rule.weights.size(), 6);
        Eigen::MatrixXd ni = Eigen::MatrixXd::Zero(5, gcoord.rows());
        std::vector<double> weights(static_cast<std::size_t>(interior_rule.weights.size()));

        for (int ig = 0; ig < interior_rule.weights.size(); ++ig)
        {
            Eigen::VectorXd n1;
            Eigen::MatrixXd dndxi1;
            Eigen::MatrixXd parent_coords = (target.elem2 < 0) ? gcoord.topRows(6) : gcoord.topRows(9);
            basis.lagrange(target.elem2 < 0 ? "T6" : "Q9", interior_rule.points.row(ig), n1, dndxi1);

            const double detj = (dndxi1.transpose() * parent_coords).determinant();
            const Eigen::Vector2d gp = parent_coords.transpose() * n1;
            const double w = interior_rule.weights(ig) * detj;
            weights[static_cast<std::size_t>(ig)] = w;

            Eigen::VectorXd poly(6);
            poly << 1.0, gp(0), gp(1), gp(0) * gp(0), gp(0) * gp(1), gp(1) * gp(1);
            p.row(ig) = poly.transpose();
            m += w * (poly * poly.transpose());

            const Eigen::VectorXd n = basis.serendipity_shape(parent_element, gcoord, gp);
            ni.row(0) += (n * w).transpose();
            ni.row(1) += (2.0 * n * w * gp(0)).transpose();
            ni.row(2) += (2.0 * n * w * gp(1)).transpose();
            ni.row(3) += (n * w * gp(0)).transpose();
            ni.row(4) += (n * w * gp(1)).transpose();
        }

        reordered_fx.row(1) -= ni.row(0);
        reordered_fx.row(3) -= ni.row(1);
        reordered_fx.row(4) -= ni.row(4);
        reordered_fy.row(2) -= ni.row(0);
        reordered_fy.row(4) -= ni.row(3);
        reordered_fy.row(5) -= ni.row(2);

        const double alpha = problem.edge_reg_param * std::max(m.trace() / 6.0, 1.0);
        const Eigen::MatrixXd mreg = m + alpha * Eigen::MatrixXd::Identity(6, 6);
        const auto mreg_lu = mreg.fullPivLu();

        SupportDomain domain;
        domain.nodes = reordered_nodes;
        domain.edof.resize(2 * reordered_nodes.size());
        for (int i = 0; i < static_cast<int>(reordered_nodes.size()); ++i)
        {
            domain.edof[static_cast<std::size_t>(2 * i)] = 2 * reordered_nodes[static_cast<std::size_t>(i)];
            domain.edof[static_cast<std::size_t>(2 * i + 1)] = 2 * reordered_nodes[static_cast<std::size_t>(i)] + 1;
        }
        domain.dx = p * mreg_lu.solve(reordered_fx);
        domain.dy = p * mreg_lu.solve(reordered_fy);
        domain.weights = std::move(weights);
        domain.sub_area = sub_area[static_cast<std::size_t>(edge_index)];
        domains.push_back(std::move(domain));
    }

    return domains;
}

std::vector<SupportDomain> CellSmoothedDomainAssembler::build_domains(const Problem &problem) const
{
    QuadratureLibrary quadrature;
    BasisEvaluator basis;

    const QuadratureRule interior_rule = quadrature.triangle(2);
    const QuadratureRule boundary_rule = quadrature.gauss_line(2);
    const int bound[3][3] = {{0, 1, 3}, {1, 2, 4}, {2, 0, 5}};

    std::vector<SupportDomain> domains;
    domains.reserve(problem.mesh.elements.size());

    for (const auto &element : problem.mesh.elements)
    {
        const Eigen::MatrixXd wkx = GeometryUtils::element_coordinates(problem.mesh, element);
        Eigen::RowVectorXd ni = Eigen::RowVectorXd::Zero(6);
        Eigen::MatrixXd wmat(3, interior_rule.weights.size());
        std::vector<double> weights(static_cast<std::size_t>(interior_rule.weights.size()));
        Eigen::MatrixXd points(interior_rule.weights.size(), 2);

        for (int ig = 0; ig < interior_rule.weights.size(); ++ig)
        {
            Eigen::VectorXd n1;
            Eigen::MatrixXd dndxi1;
            basis.lagrange("T3", interior_rule.points.row(ig), n1, dndxi1);
            const double detj = (dndxi1.transpose() * wkx.topRows(3)).determinant();
            points.row(ig) = (wkx.topRows(3).transpose() * n1).transpose();
            weights[static_cast<std::size_t>(ig)] = interior_rule.weights(ig) * detj;
            const Eigen::VectorXd n = basis.serendipity_shape("T6", wkx, points.row(ig));
            ni += (n * interior_rule.weights(ig) * detj).transpose();
            wmat.col(ig) << interior_rule.weights(ig) * detj,
                interior_rule.weights(ig) * detj * points(ig, 0),
                interior_rule.weights(ig) * detj * points(ig, 1);
        }

        std::vector<double> nx;
        std::vector<double> ny;
        GeometryUtils::outward_normals(wkx.topRows(3), GeometryUtils::side_lengths(wkx.topRows(3)), nx, ny);

        Eigen::MatrixXd fx = Eigen::MatrixXd::Zero(3, 6);
        Eigen::MatrixXd fy = Eigen::MatrixXd::Zero(3, 6);
        for (int is = 0; is < 3; ++is)
        {
            Eigen::MatrixXd bxy = Eigen::MatrixXd::Zero(3, 6);
            Eigen::MatrixXd edge_coords(3, 2);
            for (int a = 0; a < 3; ++a)
            {
                edge_coords.row(a) = wkx.row(bound[is][a]);
            }

            for (int ig = 0; ig < boundary_rule.weights.size(); ++ig)
            {
                Eigen::VectorXd ng;
                Eigen::MatrixXd dndxi;
                basis.lagrange("L3", Eigen::Vector2d(boundary_rule.points(ig, 0), 0.0), ng, dndxi);
                const Eigen::Vector2d gp = edge_coords.transpose() * ng;
                const double detj = (edge_coords.transpose() * dndxi).norm();
                const Eigen::VectorXd nt = basis.serendipity_shape("T6", wkx, gp);
                bxy.row(0) += (nt * detj * boundary_rule.weights(ig)).transpose();
                bxy.row(1) += (nt * detj * boundary_rule.weights(ig) * gp(0)).transpose();
                bxy.row(2) += (nt * detj * boundary_rule.weights(ig) * gp(1)).transpose();
            }

            fx += nx[static_cast<std::size_t>(is)] * bxy;
            fy += ny[static_cast<std::size_t>(is)] * bxy;
        }
        fx.row(1) -= ni;
        fy.row(2) -= ni;

        const auto wmat_lu = wmat.fullPivLu();

        SupportDomain domain;
        domain.nodes.assign(element.begin(), element.end());
        domain.edof.resize(12);
        for (int i = 0; i < 6; ++i)
        {
            domain.edof[static_cast<std::size_t>(2 * i)] = 2 * element[static_cast<std::size_t>(i)];
            domain.edof[static_cast<std::size_t>(2 * i + 1)] = 2 * element[static_cast<std::size_t>(i)] + 1;
        }
        domain.dx = wmat_lu.solve(fx);
        domain.dy = wmat_lu.solve(fy);
        domain.weights = std::move(weights);
        domain.sub_area = GeometryUtils::triangle_area(wkx.topRows(3));
        domains.push_back(std::move(domain));
    }

    return domains;
}

} // namespace nonlinear_esfem_t6::shared
