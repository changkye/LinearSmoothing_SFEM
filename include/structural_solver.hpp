#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <array>
#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace sfem
{

    enum class Method
    {
        FEM,
        CSFEM,
        ESFEM
    };

    enum class Scenario
    {
        LinearPatch,
        NonlinearPatch,
        Cantilever,
        BendingBlock,
        Cook
    };

    enum class LinearSolverType
    {
        SparseLU,
        BiCGSTAB,
        UmfPack
    };

    struct BoundaryCondition
    {
        std::vector<int> dofs;
        Eigen::VectorXd values;
    };

    struct Mesh
    {
        Eigen::MatrixXd nodes;
        std::vector<std::array<int, 6>> elements;
    };

    struct Model
    {
        Method method = Method::FEM;
        Scenario scenario = Scenario::LinearPatch;
        Mesh mesh;
        BoundaryCondition bc;
        Eigen::VectorXd force;
        Eigen::Matrix3d cmat = Eigen::Matrix3d::Zero();
        Eigen::Vector2d linear_material = Eigen::Vector2d::Zero();
        Eigen::Vector2d nonlinear_material = Eigen::Vector2d::Zero();
        Eigen::Vector2i num_els = Eigen::Vector2i::Ones();
        bool has_exact_solution = true;
        std::string name;
    };

    struct NewtonOptions
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

    class StructuralProblem
    {
    public:
        StructuralProblem() = default;
        explicit StructuralProblem(Model model);

        const Model &data() const;
        Model &data();

        Method method() const;
        Scenario scenario() const;
        const std::string &name() const;

        static StructuralProblem make_patch_test(Method method,
                                                 Scenario scenario,
                                                 const Eigen::Vector2i &num_els);

    private:
        Model model_;
    };

    class StructuralSolver
    {
    public:
        virtual ~StructuralSolver() = default;
        virtual Method method() const = 0;
        virtual Result solve(const StructuralProblem &problem,
                             const NewtonOptions &options) const = 0;
    };

    class FemSolver final : public StructuralSolver
    {
    public:
        Method method() const override;
        Result solve(const StructuralProblem &problem,
                     const NewtonOptions &options) const override;
    };

    class CellSmoothedFemSolver final : public StructuralSolver
    {
    public:
        Method method() const override;
        Result solve(const StructuralProblem &problem,
                     const NewtonOptions &options) const override;
    };

    class EdgeSmoothedFemSolver final : public StructuralSolver
    {
    public:
        Method method() const override;
        Result solve(const StructuralProblem &problem,
                     const NewtonOptions &options) const override;
    };

    std::unique_ptr<StructuralSolver> make_solver(Method method);

    std::string method_name(Method method);
    std::string scenario_name(Scenario scenario);
    Method parse_method(const std::string &value);
    Scenario parse_scenario(const std::string &value);
    std::vector<Scenario> available_scenarios();
    Eigen::MatrixXd exact_field(Scenario scenario, const Eigen::MatrixXd &nodes);
    BoundaryCondition make_boundary_condition(Scenario scenario, const Eigen::MatrixXd &nodes);

    Model build_model(Method method, Scenario scenario, const Eigen::Vector2i &num_els);
    Result solve_linear_patch(const Model &model, const NewtonOptions &options);
    Result solve_nonlinear_patch(const Model &model, const NewtonOptions &options);

    void write_vtu(const std::filesystem::path &file_path,
                   const Mesh &mesh,
                   const Eigen::MatrixXd &displacement,
                   const Eigen::MatrixXd &exact_displacement);

} // namespace fem
