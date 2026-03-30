#include "benchmark_utils.hpp"
#include "model_parameters.hpp"
#include "structural_solver.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

using namespace sfem;

namespace
{

    Method fixed_method()
    {
#if defined(FIXED_METHOD_CSFEM)
        return Method::CSFEM;
#elif defined(FIXED_METHOD_ESFEM)
        return Method::ESFEM;
#else
        throw std::invalid_argument("No fixed method is configured for this build");
#endif
    }

    struct OutputOptions
    {
        bool write_vtu = model_parameters::default_write_vtu;
        bool write_postprocess = model_parameters::default_write_postprocess;
        bool write_benchmark = model_parameters::default_write_benchmark;
    };

    LinearSolverType parse_solver_type(const std::string &name)
    {
        if (name == "sparselu")
        {
            return LinearSolverType::SparseLU;
        }
        if (name == "bicgstab")
        {
            return LinearSolverType::BiCGSTAB;
        }
        if (name == "umfpack")
        {
            return LinearSolverType::UmfPack;
        }
        throw std::invalid_argument("Unsupported solver type: " + name);
    }

    Eigen::Vector2i parse_num_els_value(const std::string &value)
    {
        const auto xpos = value.find_first_of("xX,");
        if (xpos == std::string::npos)
        {
            const int n = std::stoi(value);
            if (n <= 0)
            {
                throw std::invalid_argument("num-els must be positive");
            }
            return Eigen::Vector2i::Constant(n);
        }

        const std::string xpart = value.substr(0, xpos);
        const std::string ypart = value.substr(xpos + 1);
        if (xpart.empty() || ypart.empty())
        {
            throw std::invalid_argument("num-els must be formatted as n or nxm");
        }

        const int nx = std::stoi(xpart);
        const int ny = std::stoi(ypart);
        if (nx <= 0 || ny <= 0)
        {
            throw std::invalid_argument("num-els entries must be positive");
        }
        return Eigen::Vector2i(nx, ny);
    }

    void print_help()
    {
#if defined(FIXED_METHOD_CSFEM)
        constexpr const char *exe_name = "./csfem";
#elif defined(FIXED_METHOD_ESFEM)
        constexpr const char *exe_name = "./esfem";
#else
        constexpr const char *exe_name = "./sfem";
#endif
        std::cout
            << "Usage:\n"
            << "  " << exe_name << " [linear_patch|nonlinear_patch|cantilever|bending_block|cook] [options]\n\n"
            << "Options:\n"
            << "  --num-els <n|nxm>\n"
            << "  --num-els-x <n>\n"
            << "  --num-els-y <m>\n"
            << "  --output-dir <path>\n"
            << "  --nstep <value>\n"
            << "  --maxiter <value>\n"
            << "  --tol <value>\n"
            << "  --solver <sparselu|bicgstab|umfpack>\n"
            << "  --iter-maxiter <value>\n"
            << "  --iter-tol <value>\n"
            << "  --no-exact\n"
            << "  --no-vtu\n"
            << "  --summary-only\n"
            << "  --benchmark\n"
            << "  --benchmark-log <path>\n"
            << "  --debug-csfem-bending\n"
            << "  --help\n\n"
            << "Notes:\n"
            << "  cook reads data/Cook.msh directly; generate it first with data/run_mesh.sh.\n";
    }

    std::string solver_name(LinearSolverType solver)
    {
        switch (solver)
        {
        case LinearSolverType::SparseLU:
            return "sparselu";
        case LinearSolverType::BiCGSTAB:
            return "bicgstab";
        case LinearSolverType::UmfPack:
            return "umfpack";
        }
        throw std::invalid_argument("Unknown solver type");
    }

    std::string method_prefix(Method method)
    {
        switch (method)
        {
        case Method::FEM:
            return "FEM";
        case Method::CSFEM:
            return "CSFEM";
        case Method::ESFEM:
            return "ESFEM";
        }
        throw std::invalid_argument("Unknown method for output prefix");
    }

    std::string result_base_name(Method method,
                                 Scenario scenario,
                                 const Eigen::Vector2i &num_els)
    {
        if (scenario == Scenario::Cook)
        {
            return method_prefix(method) + "_cook";
        }
        return method_prefix(method) + "_" + scenario_name(scenario) + "_" +
               std::to_string(num_els(0)) + "x" + std::to_string(num_els(1));
    }

    void write_res_json(const std::filesystem::path &output_dir,
                        const std::string &base_name,
                        const Mesh &mesh,
                        const Result &result)
    {
        const std::filesystem::path json_file = output_dir / (base_name + "_res.json");
        std::ofstream out(json_file);
        if (!out)
        {
            throw std::runtime_error("Failed to open json result: " + json_file.string());
        }

        out << std::setprecision(16);
        out << "{\n";
        out << "  \"name\": \"" << base_name << "\",\n";
        out << "  \"nodes\": [\n";
        for (int i = 0; i < mesh.nodes.rows(); ++i)
        {
            out << "    {\"id\": " << i + 1
                << ", \"x\": " << mesh.nodes(i, 0)
                << ", \"y\": " << mesh.nodes(i, 1) << "}";
            out << (i + 1 == mesh.nodes.rows() ? "\n" : ",\n");
        }
        out << "  ],\n";

        out << "  \"elements\": [\n";
        for (std::size_t i = 0; i < mesh.elements.size(); ++i)
        {
            out << "    {\"id\": " << i + 1 << ", \"connectivity\": [";
            for (int a = 0; a < 6; ++a)
            {
                out << mesh.elements[i][static_cast<std::size_t>(a)] + 1;
                out << (a == 5 ? "" : ", ");
            }
            out << "]}";
            out << (i + 1 == mesh.elements.size() ? "\n" : ",\n");
        }
        out << "  ],\n";

        out << "  \"displacement\": [\n";
        for (int i = 0; i < result.nodal_u.rows(); ++i)
        {
            out << "    {\"id\": " << i + 1
                << ", \"ux\": " << result.nodal_u(i, 0)
                << ", \"uy\": " << result.nodal_u(i, 1) << "}";
            out << (i + 1 == result.nodal_u.rows() ? "\n" : ",\n");
        }
        out << "  ]\n";
        out << "}\n";
    }

    void apply_problem_specific_newton_tuning(Method method,
                                              Scenario scenario,
                                              const Eigen::Vector2i &num_els,
                                              NewtonOptions &newton)
    {
        if (method != Method::CSFEM || scenario != Scenario::BendingBlock)
        {
            return;
        }

        // Match the original MATLAB CS-FEM bending-block driver as closely as possible:
        // fixed-load increments, plain Newton updates, and convergence based on du/u.
        newton.nstep = std::max(newton.nstep, 100);
        newton.maxiter = 80;
        newton.tolerance = std::min(newton.tolerance, 1.0e-9);
        newton.residual_tolerance = std::numeric_limits<double>::infinity();
        newton.line_search_max_backtracks = 0;
        newton.adaptive_load_stepping = false;
        newton.allow_step_growth = false;
        newton.aggressive_stagnation_control = false;
    }

    void apply_large_problem_tuning(Method method,
                                    Scenario scenario,
                                    const Eigen::Vector2i &num_els,
                                    NewtonOptions &newton)
    {
        if (num_els.prod() < model_parameters::default_large_problem_element_threshold)
        {
            return;
        }

        if (scenario != Scenario::NonlinearPatch && scenario != Scenario::BendingBlock)
        {
            return;
        }

        if (method == Method::ESFEM || method == Method::CSFEM)
        {
            newton.linear_solver = parse_solver_type(
                model_parameters::default_large_problem_linear_solver);
            newton.iterative_maxiter =
                std::max(newton.iterative_maxiter,
                         model_parameters::default_large_problem_iterative_maxiter);
            newton.iterative_tolerance =
                std::max(newton.iterative_tolerance,
                         model_parameters::default_large_problem_iterative_tolerance);
        }
    }

    void write_matlab_bundle(const std::filesystem::path &output_dir,
                             const std::string &base_name,
                             const Mesh &mesh,
                             const Result &result,
                             bool has_exact_solution)
    {
        const std::filesystem::path bundle_file = output_dir / (base_name + "_result.m");
        std::ofstream out(bundle_file);
        if (!out)
        {
            throw std::runtime_error("Failed to open MATLAB bundle: " + bundle_file.string());
        }

        out << std::setprecision(16);
        out << "% Auto-generated result bundle\n";
        out << "result.name = '" << base_name << "';\n";

        out << "result.nodes = [\n";
        for (int i = 0; i < mesh.nodes.rows(); ++i)
        {
            out << "    " << i + 1 << ", " << mesh.nodes(i, 0) << ", " << mesh.nodes(i, 1) << ";\n";
        }
        out << "];\n";

        out << "result.elements = [\n";
        for (std::size_t i = 0; i < mesh.elements.size(); ++i)
        {
            out << "    " << i + 1;
            for (int a = 0; a < 6; ++a)
            {
                out << ", " << mesh.elements[i][static_cast<std::size_t>(a)] + 1;
            }
            out << ";\n";
        }
        out << "];\n";

        if (has_exact_solution && result.exact_u.rows() == result.nodal_u.rows())
        {
            out << "result.nodal_result_columns = {'node','x','y','ux','uy','ux_exact','uy_exact','ux_error','uy_error'};\n";
            out << "result.nodal_result = [\n";
            for (int i = 0; i < result.nodal_u.rows(); ++i)
            {
                out << "    " << i + 1
                    << ", " << mesh.nodes(i, 0)
                    << ", " << mesh.nodes(i, 1)
                    << ", " << result.nodal_u(i, 0)
                    << ", " << result.nodal_u(i, 1)
                    << ", " << result.exact_u(i, 0)
                    << ", " << result.exact_u(i, 1)
                    << ", " << (result.nodal_u(i, 0) - result.exact_u(i, 0))
                    << ", " << (result.nodal_u(i, 1) - result.exact_u(i, 1))
                    << ";\n";
            }
        }
        else
        {
            out << "result.nodal_result_columns = {'node','x','y','ux','uy'};\n";
            out << "result.nodal_result = [\n";
            for (int i = 0; i < result.nodal_u.rows(); ++i)
            {
                out << "    " << i + 1
                    << ", " << mesh.nodes(i, 0)
                    << ", " << mesh.nodes(i, 1)
                    << ", " << result.nodal_u(i, 0)
                    << ", " << result.nodal_u(i, 1)
                    << ";\n";
            }
        }
        out << "];\n";
    }

    void write_comparison_files(const std::filesystem::path &output_dir,
                                const std::string &base_name,
                                Method method,
                                Scenario scenario,
                                const Eigen::Vector2i &num_els,
                                const Result &result,
                                const Mesh &mesh,
                                bool has_exact_solution,
                                const OutputOptions &output_options)
    {
        (void)method;
        (void)scenario;
        (void)num_els;
        (void)has_exact_solution;
        (void)output_options;
        std::filesystem::create_directories(output_dir);

        const std::filesystem::path summary_file = output_dir / (base_name + "_output.txt");
        std::ofstream summary(summary_file);
        if (!summary)
        {
            throw std::runtime_error("Failed to open output summary: " + summary_file.string());
        }
        summary << std::setprecision(16);
        summary << "relative_error " << result.relative_error << '\n';
        summary << "strain_energy " << result.strain_energy << '\n';
        summary << "num_nodes " << result.nodal_u.rows() << '\n';
        summary << "num_steps " << result.strain_energy_history.size() << '\n';
        summary << "has_exact_solution " << (has_exact_solution ? 1 : 0) << '\n';
        for (std::size_t i = 0; i < result.strain_energy_history.size(); ++i)
        {
            summary << "strain_energy_history[" << i + 1 << "] " << result.strain_energy_history[i] << '\n';
        }

        write_res_json(output_dir, base_name, mesh, result);
    }

} // namespace

int main(int argc, char **argv)
{
    try
    {
        std::vector<std::string> args(argv + 1, argv + argc);
        for (const auto &arg : args)
        {
            if (arg == "--help" || arg == "-h")
            {
                print_help();
                return 0;
            }
        }

        Scenario scenario = parse_scenario(model_parameters::default_problem_type);
        bool scenario_set = false;
        const std::vector<Method> methods = {fixed_method()};
        Eigen::Vector2i num_els = model_parameters::default_num_els_vector();
        NewtonOptions newton;
        newton.nstep = model_parameters::default_nstep;
        newton.maxiter = model_parameters::default_maxiter;
        newton.tolerance = model_parameters::default_tolerance;
        newton.residual_tolerance = model_parameters::default_residual_tolerance;
        newton.line_search_max_backtracks =
            model_parameters::default_line_search_max_backtracks;
        newton.line_search_reduction = model_parameters::default_line_search_reduction;
        newton.line_search_min_alpha = model_parameters::default_line_search_min_alpha;
        newton.adaptive_load_stepping = model_parameters::default_adaptive_load_stepping;
        newton.max_step_cuts = model_parameters::default_max_step_cuts;
        newton.linear_solver =
            parse_solver_type(model_parameters::default_structural_linear_solver);
        newton.iterative_maxiter = model_parameters::default_structural_iterative_maxiter;
        newton.iterative_tolerance =
            model_parameters::default_structural_iterative_tolerance;
        newton.compute_exact_solution = model_parameters::default_compute_exact_solution;
        std::filesystem::path output_dir = std::filesystem::path("res");
        OutputOptions output_options;
        std::optional<std::filesystem::path> benchmark_log;
        bool solver_overridden = false;
        bool iterative_maxiter_overridden = false;
        bool iterative_tolerance_overridden = false;

        for (std::size_t i = 0; i < args.size(); ++i)
        {
            const std::string &arg = args[i];
            if (!arg.empty() && arg[0] == '-')
            {
                auto require_value = [&](const char *name) -> const std::string &
                {
                    if (i + 1 >= args.size())
                    {
                        throw std::invalid_argument(
                            std::string("Missing value for ") + name);
                    }
                    ++i;
                    return args[i];
                };

                if (arg == "--num-els")
                {
                    num_els = parse_num_els_value(require_value("--num-els"));
                }
                else if (arg == "--num-els-x")
                {
                    num_els(0) = std::stoi(require_value("--num-els-x"));
                    if (num_els(0) <= 0)
                    {
                        throw std::invalid_argument("num-els-x must be positive");
                    }
                }
                else if (arg == "--num-els-y")
                {
                    num_els(1) = std::stoi(require_value("--num-els-y"));
                    if (num_els(1) <= 0)
                    {
                        throw std::invalid_argument("num-els-y must be positive");
                    }
                }
                else if (arg == "--output-dir")
                {
                    output_dir = require_value("--output-dir");
                }
                else if (arg == "--nstep")
                {
                    newton.nstep = std::stoi(require_value("--nstep"));
                }
                else if (arg == "--maxiter")
                {
                    newton.maxiter = std::stoi(require_value("--maxiter"));
                }
                else if (arg == "--tol")
                {
                    newton.tolerance = std::stod(require_value("--tol"));
                }
                else if (arg == "--solver")
                {
                    newton.linear_solver = parse_solver_type(require_value("--solver"));
                    solver_overridden = true;
                }
                else if (arg == "--iter-maxiter")
                {
                    newton.iterative_maxiter = std::stoi(require_value("--iter-maxiter"));
                    iterative_maxiter_overridden = true;
                }
                else if (arg == "--iter-tol")
                {
                    newton.iterative_tolerance = std::stod(require_value("--iter-tol"));
                    iterative_tolerance_overridden = true;
                }
                else if (arg == "--no-exact")
                {
                    newton.compute_exact_solution = false;
                }
                else if (arg == "--no-vtu")
                {
                    output_options.write_vtu = false;
                }
                else if (arg == "--summary-only")
                {
                    output_options.write_postprocess = false;
                    output_options.write_vtu = false;
                }
                else if (arg == "--benchmark")
                {
                    output_options.write_benchmark = true;
                }
                else if (arg == "--benchmark-log")
                {
                    benchmark_log = require_value("--benchmark-log");
                    output_options.write_benchmark = true;
                }
                else if (arg == "--debug-csfem-bending")
                {
                    newton.debug_csfem_bending = true;
                }
                else
                {
                    throw std::invalid_argument("Unknown option: " + arg);
                }
            }
            else if (!scenario_set)
            {
                scenario = parse_scenario(arg);
                scenario_set = true;
            }
            else
            {
                throw std::invalid_argument("Too many positional arguments");
            }
        }

        apply_problem_specific_newton_tuning(fixed_method(), scenario, num_els, newton);
#ifdef USE_EIGEN_UMFPACK
        if (!solver_overridden &&
            scenario == Scenario::Cook &&
            (fixed_method() == Method::CSFEM || fixed_method() == Method::ESFEM))
        {
            newton.linear_solver = LinearSolverType::UmfPack;
        }
#endif
        NewtonOptions large_problem_tuned = newton;
        apply_large_problem_tuning(fixed_method(), scenario, num_els,
                                   large_problem_tuned);
        if (!solver_overridden)
        {
            newton.linear_solver = large_problem_tuned.linear_solver;
        }
        if (!iterative_maxiter_overridden)
        {
            newton.iterative_maxiter = large_problem_tuned.iterative_maxiter;
        }
        if (!iterative_tolerance_overridden)
        {
            newton.iterative_tolerance = large_problem_tuned.iterative_tolerance;
        }

        for (const Method method : methods)
        {
            const auto total_start = std::chrono::steady_clock::now();
            const StructuralProblem problem = StructuralProblem::make_patch_test(
                method, scenario, num_els);
            const auto solver = make_solver(method);
            const std::filesystem::path run_dir =
                output_dir / scenario_name(scenario) / method_name(method);
            std::filesystem::create_directories(run_dir);
            if (newton.debug_csfem_bending && newton.debug_output_dir.empty())
            {
                newton.debug_output_dir = run_dir;
            }

            const auto solve_start = std::chrono::steady_clock::now();
            const Result result = solver->solve(problem, newton);
            const auto solve_end = std::chrono::steady_clock::now();

            const std::string base_name = result_base_name(method, scenario, problem.data().num_els);
            const std::filesystem::path vtk_file = run_dir / (base_name + ".vtu");
            if (output_options.write_vtu)
            {
                write_vtu(vtk_file, problem.data().mesh, result.nodal_u, result.exact_u);
            }
            if (method == Method::CSFEM || method == Method::ESFEM)
            {
                write_comparison_files(run_dir,
                                       base_name,
                                       method,
                                       scenario,
                                       problem.data().num_els,
                                       result,
                                       problem.data().mesh,
                                       problem.data().has_exact_solution &&
                                           newton.compute_exact_solution,
                                       output_options);
            }

            const auto total_end = std::chrono::steady_clock::now();
            if (output_options.write_benchmark)
            {
                benchmark::Record record;
                record.timestamp_utc = benchmark::utc_timestamp_now();
                record.problem_type = scenario_name(scenario);
                record.method = method_name(method);
                record.num_els_x = problem.data().num_els(0);
                record.num_els_y = problem.data().num_els(1);
                record.nstep = newton.nstep;
                record.maxiter = newton.maxiter;
                record.solver = solver_name(newton.linear_solver);
                record.compute_exact_solution = newton.compute_exact_solution;
                record.write_vtu = output_options.write_vtu;
                record.write_postprocess = output_options.write_postprocess;
                record.relative_error = result.relative_error;
                record.strain_energy = result.strain_energy;
                record.solve_wall_time_sec =
                    benchmark::seconds_between(solve_start, solve_end);
                record.total_wall_time_sec =
                    benchmark::seconds_between(total_start, total_end);
                record.peak_rss_mb = benchmark::peak_rss_mb();
                const std::filesystem::path log_path =
                    benchmark_log.value_or(run_dir / (base_name + "_benchmark.csv"));
                benchmark::append_csv(log_path, record);
            }

            std::cout << "--------------------------------------------------\n";
            std::cout << std::fixed << std::setprecision(6)
                      << "Final Timing [s] solve="
                      << benchmark::seconds_between(solve_start, solve_end)
                      << " total="
                      << benchmark::seconds_between(total_start, total_end)
                      << '\n';
            std::cout << scenario_name(scenario) << " / " << method_name(method)
                      << ": rel_error=" << result.relative_error
                      << ", strain_energy=" << result.strain_energy
                      << ", output=" << (output_options.write_vtu ? vtk_file.string() : run_dir.string()) << '\n';
        }

        return 0;
    }
    catch (const std::exception &ex)
    {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }
}
