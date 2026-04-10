#include "model_parameters.hpp"
#include "structural_problem_library.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <set>
#include <unordered_map>

namespace sfem
{
    namespace
    {

        constexpr double kCoordTol = 1.0e-10;
        constexpr double kCookTotalLoadY = 1.0 / 16.0;

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

            auto add_family = [&](const std::array<int, 6> &base)
            {
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

        int gmsh_nodes_per_element(const int element_type)
        {
            switch (element_type)
            {
            case 1:
                return 2; // 2-node line
            case 2:
                return 3; // 3-node triangle
            case 8:
                return 3; // 3-node second-order line
            case 9:
                return 6; // 6-node second-order triangle
            case 15:
                return 1; // point
            default:
                throw std::runtime_error("Unsupported Gmsh element type in Cook mesh: " + std::to_string(element_type));
            }
        }

        std::filesystem::path locate_data_file(const std::string &file_name)
        {
            const std::array<std::filesystem::path, 4> candidates = {
                std::filesystem::path(file_name),
                std::filesystem::path("data") / file_name,
                std::filesystem::path("..") / "data" / file_name,
                std::filesystem::path("..") / ".." / "data" / file_name};

            for (const auto &candidate : candidates)
            {
                if (std::filesystem::exists(candidate))
                {
                    return candidate;
                }
            }

            throw std::runtime_error(
                "Could not find '" + file_name + "'. Generate it first with 'cd data && ./run_mesh.sh' and choose 'Cook.geo'.");
        }

        Mesh load_gmsh_t6_mesh(const std::filesystem::path &file_path)
        {
            std::ifstream in(file_path);
            if (!in)
            {
                throw std::runtime_error("Failed to open Gmsh mesh: " + file_path.string());
            }

            std::string token;
            auto require_token = [&](const std::string &expected)
            {
                if (!(in >> token) || token != expected)
                {
                    throw std::runtime_error("Malformed Gmsh mesh '" + file_path.string() + "': expected " + expected);
                }
            };

            require_token("$MeshFormat");
            double version = 0.0;
            int file_type = -1;
            int data_size = 0;
            in >> version >> file_type >> data_size;
            if (!in)
            {
                throw std::runtime_error("Malformed Gmsh mesh header in: " + file_path.string());
            }
            if (file_type != 0)
            {
                throw std::runtime_error("Only ASCII Gmsh .msh files are supported: " + file_path.string());
            }
            require_token("$EndMeshFormat");

            std::vector<Eigen::Vector2d> coordinates;
            std::unordered_map<std::int64_t, int> node_index;
            std::vector<std::array<int, 6>> elements;

            while (in >> token)
            {
                if (token == "$Nodes")
                {
                    std::size_t num_entity_blocks = 0;
                    std::size_t num_nodes = 0;
                    std::int64_t min_node_tag = 0;
                    std::int64_t max_node_tag = 0;
                    in >> num_entity_blocks >> num_nodes >> min_node_tag >> max_node_tag;
                    if (!in)
                    {
                        throw std::runtime_error("Malformed $Nodes header in: " + file_path.string());
                    }

                    coordinates.reserve(num_nodes);
                    node_index.reserve(num_nodes);

                    for (std::size_t block = 0; block < num_entity_blocks; ++block)
                    {
                        int entity_dim = 0;
                        int entity_tag = 0;
                        int parametric = 0;
                        std::size_t num_nodes_in_block = 0;
                        in >> entity_dim >> entity_tag >> parametric >> num_nodes_in_block;
                        if (!in)
                        {
                            throw std::runtime_error("Malformed node block in: " + file_path.string());
                        }

                        std::vector<std::int64_t> block_tags(num_nodes_in_block);
                        for (std::size_t i = 0; i < num_nodes_in_block; ++i)
                        {
                            in >> block_tags[i];
                        }

                        for (std::size_t i = 0; i < num_nodes_in_block; ++i)
                        {
                            double x = 0.0;
                            double y = 0.0;
                            double z = 0.0;
                            in >> x >> y >> z;
                            for (int p = 0; p < (parametric ? entity_dim : 0); ++p)
                            {
                                double ignored = 0.0;
                                in >> ignored;
                            }
                            if (!in)
                            {
                                throw std::runtime_error("Malformed node coordinates in: " + file_path.string());
                            }

                            const int index = static_cast<int>(coordinates.size());
                            coordinates.emplace_back(x, y);
                            node_index.emplace(block_tags[i], index);
                        }
                    }

                    require_token("$EndNodes");
                    continue;
                }

                if (token == "$Elements")
                {
                    std::size_t num_entity_blocks = 0;
                    std::size_t num_elements = 0;
                    std::int64_t min_element_tag = 0;
                    std::int64_t max_element_tag = 0;
                    in >> num_entity_blocks >> num_elements >> min_element_tag >> max_element_tag;
                    if (!in)
                    {
                        throw std::runtime_error("Malformed $Elements header in: " + file_path.string());
                    }

                    elements.reserve(num_elements);

                    for (std::size_t block = 0; block < num_entity_blocks; ++block)
                    {
                        int entity_dim = 0;
                        int entity_tag = 0;
                        int element_type = 0;
                        std::size_t num_elements_in_block = 0;
                        in >> entity_dim >> entity_tag >> element_type >> num_elements_in_block;
                        if (!in)
                        {
                            throw std::runtime_error("Malformed element block in: " + file_path.string());
                        }

                        const int nodes_per_element = gmsh_nodes_per_element(element_type);
                        for (std::size_t i = 0; i < num_elements_in_block; ++i)
                        {
                            std::int64_t element_tag = 0;
                            in >> element_tag;
                            if (!in)
                            {
                                throw std::runtime_error("Malformed element entry in: " + file_path.string());
                            }

                            if (element_type == 9)
                            {
                                std::array<int, 6> element{};
                                for (int a = 0; a < 6; ++a)
                                {
                                    std::int64_t node_tag = 0;
                                    in >> node_tag;
                                    const auto it = node_index.find(node_tag);
                                    if (it == node_index.end())
                                    {
                                        throw std::runtime_error("Element references unknown node tag in: " + file_path.string());
                                    }
                                    element[static_cast<std::size_t>(a)] = it->second;
                                }
                                elements.push_back(element);
                            }
                            else
                            {
                                for (int a = 0; a < nodes_per_element; ++a)
                                {
                                    std::int64_t ignored_node = 0;
                                    in >> ignored_node;
                                }
                            }
                        }
                    }

                    require_token("$EndElements");
                    continue;
                }

                if (token.rfind("$End", 0) == 0)
                {
                    continue;
                }

                if (!token.empty() && token[0] == '$')
                {
                    const std::string end_token = "$End" + token.substr(1);
                    do
                    {
                        if (!(in >> token))
                        {
                            throw std::runtime_error("Unterminated Gmsh section in: " + file_path.string());
                        }
                    } while (token != end_token);
                }
            }

            if (coordinates.empty())
            {
                throw std::runtime_error("No nodes were found in Gmsh mesh: " + file_path.string());
            }
            if (elements.empty())
            {
                throw std::runtime_error(
                    "No T6 elements were found in Gmsh mesh: " + file_path.string() +
                    ". Make sure run_mesh.sh generates a quadratic mesh with '-order 2'.");
            }

            for (auto &element : elements)
            {
                const Eigen::Vector2d &x1 = coordinates[static_cast<std::size_t>(element[0])];
                const Eigen::Vector2d &x2 = coordinates[static_cast<std::size_t>(element[1])];
                const Eigen::Vector2d &x3 = coordinates[static_cast<std::size_t>(element[2])];
                const double signed_area2 =
                    (x2(0) - x1(0)) * (x3(1) - x1(1)) - (x3(0) - x1(0)) * (x2(1) - x1(1));
                if (signed_area2 < 0.0)
                {
                    element = {element[0], element[2], element[1], element[5], element[4], element[3]};
                }
            }

            std::set<int> used_node_indices;
            for (const auto &element : elements)
            {
                for (const int node : element)
                {
                    used_node_indices.insert(node);
                }
            }

            std::unordered_map<int, int> compact_index;
            compact_index.reserve(used_node_indices.size());
            Mesh mesh;
            mesh.nodes.resize(static_cast<int>(used_node_indices.size()), 2);

            int compact_counter = 0;
            for (const int old_index : used_node_indices)
            {
                compact_index.emplace(old_index, compact_counter);
                mesh.nodes.row(compact_counter) = coordinates[static_cast<std::size_t>(old_index)];
                ++compact_counter;
            }

            mesh.elements.reserve(elements.size());
            for (const auto &element : elements)
            {
                std::array<int, 6> compact_element{};
                for (int a = 0; a < 6; ++a)
                {
                    compact_element[static_cast<std::size_t>(a)] = compact_index.at(element[static_cast<std::size_t>(a)]);
                }
                mesh.elements.push_back(compact_element);
            }

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
            std::transform(value.begin(), value.end(), value.begin(), [](unsigned char c)
                           { return static_cast<char>(std::tolower(c)); });
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

        BoundaryCondition make_zero_boundary_condition(const std::vector<int> &fixed_nodes)
        {
            BoundaryCondition bc;
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

        Eigen::VectorXd make_force_vector(const Mesh &mesh)
        {
            return Eigen::VectorXd::Zero(2 * mesh.nodes.rows());
        }

    } // namespace

    std::vector<Scenario> ProblemLibrary::available()
    {
        return {Scenario::LinearPatch, Scenario::HighOrderPatch, Scenario::ModeTest};
    }

    std::string ProblemLibrary::name(Scenario scenario)
    {
        switch (scenario)
        {
        case Scenario::LinearPatch:
            return "linear_patch";
        case Scenario::HighOrderPatch:
            return "high_order_patch";
        case Scenario::ModeTest:
            return "mode_test";
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
        if (lower == "high_order_patch")
        {
            return Scenario::HighOrderPatch;
        }
        if (lower == "mode_test")
        {
            return Scenario::ModeTest;
        }
        throw std::invalid_argument("Unsupported scenario: " + value);
    }

    Eigen::MatrixXd ProblemLibrary::exact_field(Scenario scenario, const Eigen::MatrixXd &nodes)
    {
        if (scenario == Scenario::ModeTest)
        {
            throw std::logic_error("ModeTest has no exact displacement field");
        }
        Eigen::MatrixXd u = Eigen::MatrixXd::Zero(nodes.rows(), 2);
        for (int i = 0; i < nodes.rows(); ++i)
        {
            const double x = nodes(i, 0);
            const double y = nodes(i, 1);
            if (scenario == Scenario::HighOrderPatch)
            {
                // Quadratic displacement field satisfying equilibrium (E=1, nu=0, no body forces):
                //   d(sigma_xx)/dx + d(sigma_xy)/dy = 0.2 + (-0.2) = 0
                //   d(sigma_xy)/dx + d(sigma_yy)/dy = 0   + 0      = 0
                u(i, 0) = 0.1 + 0.1 * x + 0.2 * y + 0.1 * x * x + 0.05 * y * y;
                u(i, 1) = 0.05 + 0.15 * x + 0.1 * y - 0.5 * x * y;
            }
            else
            {
                // Linear displacement field (LinearPatch)
                u(i, 0) = 0.1 + 0.1 * x + 0.2 * y;
                u(i, 1) = 0.05 + 0.15 * x + 0.1 * y;
            }
        }
        return u;
    }

    BoundaryCondition ProblemLibrary::make_boundary_condition(Scenario scenario, const Eigen::MatrixXd &nodes)
    {
        BoundaryCondition bc;
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

        const model_parameters::LinearElasticMaterial linear_material =
            model_parameters::linear_material_for(name(scenario));

        if (scenario == Scenario::ModeTest)
        {
            model.mesh = load_gmsh_t6_mesh(locate_data_file("mode_test.msh"));
            model.num_els = Eigen::Vector2i(static_cast<int>(model.mesh.elements.size()), 0);
            model.linear_material << linear_material.young, linear_material.poisson;
            model.cmat = linear_constitutive(linear_material.young, linear_material.poisson);
            model.has_exact_solution = false;
            // No BCs — free body to expose all zero-energy modes
            model.bc = BoundaryCondition{};
            model.force = make_force_vector(model.mesh);
            return model;
        }

        if (scenario == Scenario::HighOrderPatch)
        {
            model.mesh = load_gmsh_t6_mesh(locate_data_file("high_order_patch.msh"));
            model.num_els = Eigen::Vector2i(static_cast<int>(model.mesh.elements.size()), 0);
            model.linear_material << linear_material.young, linear_material.poisson;
            model.cmat = linear_constitutive(linear_material.young, linear_material.poisson);
            model.has_exact_solution = true;
            model.bc = make_boundary_condition(scenario, model.mesh.nodes);
            model.force = make_force_vector(model.mesh);
            return model;
        }

        model.mesh = load_gmsh_t6_mesh(locate_data_file("linear_patch.msh"));
        model.num_els = Eigen::Vector2i(static_cast<int>(model.mesh.elements.size()), 0);
        model.linear_material << linear_material.young, linear_material.poisson;
        model.cmat = linear_constitutive(linear_material.young, linear_material.poisson);
        model.has_exact_solution = true;
        model.bc = make_boundary_condition(scenario, model.mesh.nodes);
        model.force = make_force_vector(model.mesh);
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
