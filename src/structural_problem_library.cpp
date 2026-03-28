#include "model_parameters.hpp"
#include "structural_problem_library.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <stdexcept>

namespace fem
{
namespace
{

constexpr double kCoordTol = 1.0e-10;

Mesh make_t6_mesh(const Eigen::Matrix<double, 2, 2> &limits, const Eigen::Vector2i &num_els)
{
    const int xn = 2 * num_els(0) + 1;
    const int yn = 2 * num_els(1) + 1;

    Mesh mesh;
    mesh.nodes.resize(xn * yn, 2);

    int idx = 0;
    for (int j = 0; j < yn; ++j)
    {
        const double y = limits(1, 0) + (limits(1, 1) - limits(1, 0)) * static_cast<double>(j) / static_cast<double>(yn - 1);
        for (int i = 0; i < xn; ++i)
        {
            const double x = limits(0, 0) + (limits(0, 1) - limits(0, 0)) * static_cast<double>(i) / static_cast<double>(xn - 1);
            mesh.nodes.row(idx++) << x, y;
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

Eigen::Matrix3d linear_constitutive(double e, double nu)
{
    Eigen::Matrix3d cmat = Eigen::Matrix3d::Zero();
    const double factor = e / (1.0 - nu * nu);
    cmat << 1.0, nu, 0.0,
        nu, 1.0, 0.0,
        0.0, 0.0, 0.5 * (1.0 - nu);
    return factor * cmat;
}

std::string to_lower(std::string value)
{
    std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    });
    return value;
}

std::vector<int> boundary_nodes(const Eigen::MatrixXd &nodes)
{
    const double xmin = nodes.col(0).minCoeff();
    const double xmax = nodes.col(0).maxCoeff();
    const double ymin = nodes.col(1).minCoeff();
    const double ymax = nodes.col(1).maxCoeff();
    std::vector<int> selected;
    for (int i = 0; i < nodes.rows(); ++i)
    {
        const double x = nodes(i, 0);
        const double y = nodes(i, 1);
        if (std::abs(x - xmin) < kCoordTol || std::abs(x - xmax) < kCoordTol ||
            std::abs(y - ymin) < kCoordTol || std::abs(y - ymax) < kCoordTol)
        {
            selected.push_back(i);
        }
    }
    return selected;
}

bool is_nonlinear_scenario(Scenario scenario)
{
    return scenario == Scenario::NonlinearPatch || scenario == Scenario::BendingBlock;
}

std::vector<int> left_edge_nodes(const Eigen::MatrixXd &nodes)
{
    const double xmin = nodes.col(0).minCoeff();
    std::vector<int> selected;
    for (int i = 0; i < nodes.rows(); ++i)
    {
        if (std::abs(nodes(i, 0) - xmin) < kCoordTol)
        {
            selected.push_back(i);
        }
    }
    return selected;
}

std::vector<int> right_edge_nodes(const Eigen::MatrixXd &nodes)
{
    const double xmax = nodes.col(0).maxCoeff();
    std::vector<int> selected;
    for (int i = 0; i < nodes.rows(); ++i)
    {
        if (std::abs(nodes(i, 0) - xmax) < kCoordTol)
        {
            selected.push_back(i);
        }
    }
    return selected;
}

Eigen::VectorXd make_force_vector(Scenario scenario, const Mesh &mesh)
{
    Eigen::VectorXd force = Eigen::VectorXd::Zero(2 * mesh.nodes.rows());
    if (scenario != Scenario::Cantilever)
    {
        return force;
    }

    const std::vector<int> loaded_nodes = right_edge_nodes(mesh.nodes);
    if (loaded_nodes.empty())
    {
        return force;
    }

    const double nodal_force = -1.0 / static_cast<double>(loaded_nodes.size());
    for (const int node : loaded_nodes)
    {
        force(2 * node + 1) += nodal_force;
    }
    return force;
}

} // namespace

std::vector<Scenario> ProblemLibrary::available()
{
    return {Scenario::LinearPatch, Scenario::NonlinearPatch, Scenario::Cantilever, Scenario::BendingBlock};
}

std::string ProblemLibrary::name(Scenario scenario)
{
    switch (scenario)
    {
    case Scenario::LinearPatch:
        return "linear_patch";
    case Scenario::NonlinearPatch:
        return "nonlinear_patch";
    case Scenario::Cantilever:
        return "cantilever";
    case Scenario::BendingBlock:
        return "bending_block";
    }
    throw std::invalid_argument("Unknown scenario");
}

Scenario ProblemLibrary::parse(const std::string &value)
{
    const std::string lower = to_lower(value);
    if (lower == "linear_patch")
    {
        return Scenario::LinearPatch;
    }
    if (lower == "nonlinear_patch")
    {
        return Scenario::NonlinearPatch;
    }
    if (lower == "cantilever")
    {
        return Scenario::Cantilever;
    }
    if (lower == "bending_block")
    {
        return Scenario::BendingBlock;
    }
    throw std::invalid_argument("Unsupported scenario: " + value);
}

Eigen::MatrixXd ProblemLibrary::exact_field(Scenario scenario, const Eigen::MatrixXd &nodes)
{
    Eigen::MatrixXd u = Eigen::MatrixXd::Zero(nodes.rows(), 2);
    for (int i = 0; i < nodes.rows(); ++i)
    {
        const double x = nodes(i, 0);
        const double y = nodes(i, 1);
        if (scenario == Scenario::LinearPatch)
        {
            u(i, 0) = 0.1 + 0.1 * x + 0.2 * y;
            u(i, 1) = 0.05 + 0.15 * x + 0.1 * y;
        }
        else if (scenario == Scenario::NonlinearPatch)
        {
            u(i, 0) = 0.5 * y * y;
            u(i, 1) = 0.0;
        }
        else if (scenario == Scenario::BendingBlock)
        {
            constexpr double alpha = 0.9;
            u(i, 0) = std::sqrt(2.0 * alpha * x) * std::cos(y / alpha) - x;
            u(i, 1) = std::sqrt(2.0 * alpha * x) * std::sin(y / alpha) - y;
        }
    }
    return u;
}

BoundaryCondition ProblemLibrary::make_boundary_condition(Scenario scenario, const Eigen::MatrixXd &nodes)
{
    BoundaryCondition bc;
    if (scenario == Scenario::Cantilever)
    {
        const std::vector<int> fixed_nodes = left_edge_nodes(nodes);
        bc.dofs.reserve(fixed_nodes.size() * 2);
        bc.values = Eigen::VectorXd::Zero(static_cast<int>(fixed_nodes.size() * 2));
        int idx = 0;
        for (const int node : fixed_nodes)
        {
            bc.dofs.push_back(2 * node);
            bc.dofs.push_back(2 * node + 1);
            bc.values(idx++) = 0.0;
            bc.values(idx++) = 0.0;
        }
        return bc;
    }

    const std::vector<int> fixed_nodes = boundary_nodes(nodes);
    const Eigen::MatrixXd exact = exact_field(scenario, nodes);
    bc.dofs.reserve(fixed_nodes.size() * 2);
    bc.values.resize(static_cast<int>(fixed_nodes.size() * 2));
    int idx = 0;
    for (const int node : fixed_nodes)
    {
        bc.dofs.push_back(2 * node);
        bc.dofs.push_back(2 * node + 1);
        bc.values(idx++) = exact(node, 0);
        bc.values(idx++) = exact(node, 1);
    }
    return bc;
}

Model ProblemLibrary::build_model(Method method, Scenario scenario, const Eigen::Vector2i &num_els)
{
    Model model;
    model.method = method;
    model.scenario = scenario;
    model.num_els = num_els;
    model.name = name(scenario) + "_" + method_name(method);

    const std::string problem_type = name(scenario);
    const model_parameters::Geometry geometry = model_parameters::geometry_for(problem_type);
    const model_parameters::LinearElasticMaterial linear_material = model_parameters::linear_material_for(problem_type);
    const model_parameters::NeoHookeanMaterial nonlinear_material = model_parameters::nonlinear_material_for(problem_type);

    Eigen::Matrix<double, 2, 2> limits;
    if (scenario == Scenario::LinearPatch)
    {
        limits << 0.0, geometry.L,
            0.0, geometry.H;
        model.linear_material << linear_material.young, linear_material.poisson;
        model.has_exact_solution = true;
    }
    else if (is_nonlinear_scenario(scenario))
    {
        if (scenario == Scenario::BendingBlock)
        {
            limits << 2.0, 3.0,
                -2.0, 2.0;
        }
        else
        {
            limits << 0.0, geometry.L,
                0.0, geometry.H;
        }
        model.nonlinear_material << nonlinear_material.mu, nonlinear_material.lambda_like;
        model.has_exact_solution = true;
    }
    else
    {
        limits << 0.0, geometry.L,
            0.0, geometry.H;
        model.linear_material << linear_material.young, linear_material.poisson;
        model.has_exact_solution = false;
    }

    model.mesh = make_t6_mesh(limits, num_els);
    model.bc = make_boundary_condition(scenario, model.mesh.nodes);
    model.force = make_force_vector(scenario, model.mesh);
    if (!is_nonlinear_scenario(scenario))
    {
        model.cmat = linear_constitutive(model.linear_material(0), model.linear_material(1));
    }
    return model;
}

std::string scenario_name(Scenario scenario)
{
    return ProblemLibrary::name(scenario);
}

Scenario parse_scenario(const std::string &value)
{
    return ProblemLibrary::parse(value);
}

std::vector<Scenario> available_scenarios()
{
    return ProblemLibrary::available();
}

Eigen::MatrixXd exact_field(Scenario scenario, const Eigen::MatrixXd &nodes)
{
    return ProblemLibrary::exact_field(scenario, nodes);
}

BoundaryCondition make_boundary_condition(Scenario scenario, const Eigen::MatrixXd &nodes)
{
    return ProblemLibrary::make_boundary_condition(scenario, nodes);
}

Model build_model(Method method, Scenario scenario, const Eigen::Vector2i &num_els)
{
    return ProblemLibrary::build_model(method, scenario, num_els);
}

} // namespace fem
