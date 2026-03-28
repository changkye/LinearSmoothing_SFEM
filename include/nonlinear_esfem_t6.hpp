#pragma once

#include <Eigen/Dense>

#include <array>
#include <filesystem>
#include <string>
#include <vector>

namespace nonlinear_esfem_t6
{

enum class LinearSolverType
{
    SparseLU,
    BiCGSTAB,
    UmfPack
};

struct Mesh
{
    Eigen::MatrixXd nodes;
    std::vector<std::array<int, 6>> elements;
};

struct BoundaryCondition
{
    std::vector<int> dofs;
    Eigen::VectorXd values;
};

struct Problem
{
    Mesh mesh;
    BoundaryCondition bc;
    Eigen::VectorXd force;
    std::string problem_type;
    Eigen::Vector2d material = Eigen::Vector2d::Zero();
    Eigen::Vector2i num_els = Eigen::Vector2i::Ones();
    int edge_boundary_ng = 3;
    int edge_tri_quad_order = 4;
    int edge_quad_quad_order = 3;
    double edge_reg_param = 1.0e-10;
};

struct SolverOptions
{
    int nstep = 100;
    int maxiter = 80;
    double tolerance = 1.0e-9;
    LinearSolverType linear_solver = LinearSolverType::SparseLU;
    int iterative_maxiter = 2000;
    double iterative_tolerance = 1.0e-10;
    double residual_tolerance = 1.0e-6;
    int line_search_max_backtracks = 6;
    double line_search_reduction = 0.5;
    double line_search_min_alpha = 1.0 / 64.0;
    bool adaptive_load_stepping = true;
    int max_step_cuts = 10;
    bool allow_step_growth = true;
    bool aggressive_stagnation_control = false;
    bool compute_exact_solution = true;
    bool debug_csfem_bending = false;
    std::filesystem::path debug_output_dir;
};

struct Result
{
    Eigen::VectorXd uu;
    Eigen::MatrixXd nodal_u;
    Eigen::MatrixXd exact_u;
    std::vector<double> strain_energy_history;
    double relative_error = 0.0;
    double strain_energy = 0.0;
};

class NonlinearEsfemT6Problem
{
public:
    static NonlinearEsfemT6Problem make_not_so_simple_shear(const Eigen::Vector2i &num_els,
                                                            const Eigen::Vector2d &material);

    const Problem &data() const;
    Problem &data();

private:
    explicit NonlinearEsfemT6Problem(Problem problem);
    Problem problem_;
};

class NonlinearEsfemT6Solver
{
public:
    Result solve(const NonlinearEsfemT6Problem &problem,
                 const SolverOptions &options) const;
};

class NonlinearCsfemT6Solver
{
public:
    Result solve(const NonlinearEsfemT6Problem &problem,
                 const SolverOptions &options) const;
};

void write_vtu(const std::filesystem::path &file_path,
               const Mesh &mesh,
               const Eigen::MatrixXd &displacement,
               const Eigen::MatrixXd &exact_displacement);

} // namespace nonlinear_esfem_t6
