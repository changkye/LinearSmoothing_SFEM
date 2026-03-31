#pragma once

#include <Eigen/Dense>

#include <array>
#include <stdexcept>
#include <string_view>

namespace model_parameters
{
    struct Geometry
    {
        double L;
        double H;
    };
    struct LinearElasticMaterial
    {
        double young;
        double poisson;
    };
    struct NeoHookeanMaterial
    {
        double mu;
        double lambda_like;
    };

    // 1) Problem Type
    inline constexpr char default_problem_type[] = "linear_patch";

    // 2) Geometry
    inline Geometry geometry_for(std::string_view problem_type);

    inline Geometry default_geometry() { return geometry_for(default_problem_type); }

    inline Geometry geometry_for(std::string_view problem_type)
    {
        if (problem_type == "linear_patch")
        {
            return Geometry{1.0, 1.0};
        }
        if (problem_type == "nonlinear_patch")
        {
            return Geometry{2.0, 2.0};
        }
        if (problem_type == "cantilever")
        {
            return Geometry{4.0, 1.0};
        }
        if (problem_type == "bending_block")
        {
            return Geometry{1.0, 4.0};
        }
        throw std::invalid_argument("Unsupported problem type in geometry_for");
    }

    // 3) Number Of Elements
    // Users can edit this array directly.
    inline constexpr std::array<int, 2> default_num_els{{20, 20}};

    inline Eigen::Vector2i default_num_els_vector() { return Eigen::Vector2i(default_num_els[0], default_num_els[1]); }

    // 4) Material Properties

    // Default material values shared across problems unless a case-specific
    // override is added below.
    inline constexpr LinearElasticMaterial default_linear_material{1.0, 0.3};
    inline constexpr NeoHookeanMaterial default_nonlinear_material{0.6, 100};

    inline LinearElasticMaterial linear_material_for(std::string_view problem_type)
    {
        LinearElasticMaterial material = default_linear_material;
        if (problem_type == "linear_patch")
        {
            material.poisson = 0.0;
        }
        // Cook* problems: mu=1.0, kappa=1000 => E=9*kappa*mu/(3*kappa+mu), nu=(3*kappa-2*mu)/(2*(3*kappa+mu))
        if (problem_type.substr(0, 4) == "cook")
        {
            material.young   = 9000.0 / 3001.0;
            material.poisson = 1499.0 / 3001.0;
        }
        return material;
    }

    inline NeoHookeanMaterial nonlinear_material_for(std::string_view problem_type)
    {
        if (problem_type.substr(0, 4) == "cook")
        {
            return NeoHookeanMaterial{1.0, 1000.0};
        }
        if (problem_type == "cantilever")
        {
            const LinearElasticMaterial linear = linear_material_for(problem_type);
            const double mu = linear.young / (2.0 * (1.0 + linear.poisson));
            const double lambda_like =
                linear.young * linear.poisson / ((1.0 + linear.poisson) * (1.0 - 2.0 * linear.poisson));
            return NeoHookeanMaterial{mu, lambda_like};
        }
        return default_nonlinear_material;
    }

    // 5) ES-FEM Stabilization / Quadrature Parameters
    inline constexpr int esfem_edge_boundary_ng = 3;
    inline constexpr int esfem_edge_tri_quad_order = 4;
    inline constexpr int esfem_edge_quad_quad_order = 3;
    inline constexpr double esfem_edge_reg_param = 1.0e-10;

    // 6) Newton-Raphson Iteration Settings
    inline constexpr int default_nstep = 100;
    inline constexpr int default_maxiter = 10;
    inline constexpr double default_tolerance = 1.0e-9;
    inline constexpr double default_residual_tolerance = 1.0e-6;
    inline constexpr int default_line_search_max_backtracks = 6;
    inline constexpr double default_line_search_reduction = 0.5;
    inline constexpr double default_line_search_min_alpha = 1.0 / 64.0;
    inline constexpr bool default_adaptive_load_stepping = true;
    inline constexpr int default_max_step_cuts = 10;

    // 7) Others
    inline constexpr char default_structural_linear_solver[] = "sparselu";
    inline constexpr int default_structural_iterative_maxiter = 2000;
    inline constexpr double default_structural_iterative_tolerance = 1.0e-10;
    // why esfem has its own default linear solver and others?
    inline constexpr char default_esfem_linear_solver[] = "sparselu";
    inline constexpr int default_esfem_iterative_maxiter = 2000;
    inline constexpr double default_esfem_iterative_tolerance = 1.0e-10;
    inline constexpr int default_large_problem_element_threshold = 200;
    inline constexpr char default_large_problem_linear_solver[] = "bicgstab";
    inline constexpr int default_large_problem_iterative_maxiter = 8000;
    inline constexpr double default_large_problem_iterative_tolerance = 1.0e-9;
    inline constexpr bool default_compute_exact_solution = true;
    inline constexpr bool default_write_vtu = true;
    inline constexpr bool default_write_postprocess = true;
    inline constexpr bool default_write_benchmark = false;

} // namespace model_parameters
