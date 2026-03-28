#include "benchmark_utils.hpp"
#include "model_parameters.hpp"
#include "nonlinear_esfem_t6.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

struct OutputOptions
{
    bool write_vtu = model_parameters::default_write_vtu;
    bool write_postprocess = model_parameters::default_write_postprocess;
    bool write_benchmark = model_parameters::default_write_benchmark;
};

nonlinear_esfem_t6::LinearSolverType parse_solver_type(const std::string &name)
{
    if (name == "sparselu")
    {
        return nonlinear_esfem_t6::LinearSolverType::SparseLU;
    }
    if (name == "bicgstab")
    {
        return nonlinear_esfem_t6::LinearSolverType::BiCGSTAB;
    }
    if (name == "umfpack")
    {
        return nonlinear_esfem_t6::LinearSolverType::UmfPack;
    }
    throw std::invalid_argument("Unsupported solver type: " + name);
}

void print_help()
{
    std::cout
        << "Usage:\n"
        << "  ./nonlinear_csfem_t6 [nonlinear_patch] [options]\n\n"
        << "Options:\n"
        << "  --num-els <n>\n"
        << "  --mu <value>\n"
        << "  --lambda_like <value>\n"
        << "  --nstep <value>\n"
        << "  --maxiter <value>\n"
        << "  --tol <value>\n"
        << "  --solver <sparselu|bicgstab|umfpack>\n"
        << "  --iter-maxiter <value>\n"
        << "  --iter-tol <value>\n"
        << "  --output-dir <path>\n"
        << "  --no-exact\n"
        << "  --no-vtu\n"
        << "  --summary-only\n"
        << "  --benchmark\n"
        << "  --benchmark-log <path>\n"
        << "  --help\n";
}

std::string solver_name(nonlinear_esfem_t6::LinearSolverType solver)
{
    switch (solver)
    {
    case nonlinear_esfem_t6::LinearSolverType::SparseLU:
        return "sparselu";
    case nonlinear_esfem_t6::LinearSolverType::BiCGSTAB:
        return "bicgstab";
    case nonlinear_esfem_t6::LinearSolverType::UmfPack:
        return "umfpack";
    }
    throw std::invalid_argument("Unknown solver type");
}

std::string csfem_base_name(const std::string &problem_type, int num_els_x, int num_els_y)
{
    return "CSFEM_" + problem_type + "_" + std::to_string(num_els_x) + "x" + std::to_string(num_els_y);
}

void write_matlab_bundle(const std::filesystem::path &output_dir,
                         const std::string &base_name,
                         const nonlinear_esfem_t6::Mesh &mesh,
                         const nonlinear_esfem_t6::Result &result,
                         bool has_exact)
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

    if (has_exact)
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

void write_matlab_files(const std::filesystem::path &output_dir,
                        const std::string &base_name,
                        const nonlinear_esfem_t6::Mesh &mesh,
                        const nonlinear_esfem_t6::Result &result,
                        const OutputOptions &output_options)
{
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
    for (std::size_t i = 0; i < result.strain_energy_history.size(); ++i)
    {
        summary << "strain_energy_history[" << i + 1 << "] " << result.strain_energy_history[i] << '\n';
    }

    if (!output_options.write_postprocess)
    {
        return;
    }

    const bool has_exact = (result.exact_u.rows() == result.nodal_u.rows() &&
                            result.exact_u.cols() == result.nodal_u.cols());
    write_matlab_bundle(output_dir, base_name, mesh, result, has_exact);

    if (has_exact)
    {
        const std::filesystem::path comparison_file = output_dir / (base_name + "_comparison.csv");
        std::ofstream comparison(comparison_file);
        if (!comparison)
        {
            throw std::runtime_error("Failed to open comparison csv: " + comparison_file.string());
        }
        comparison << std::setprecision(16);
        comparison << "node,ux,uy,ux_exact,uy_exact,ux_error,uy_error\n";
        for (int i = 0; i < result.nodal_u.rows(); ++i)
        {
            comparison << i + 1 << ','
                       << result.nodal_u(i, 0) << ','
                       << result.nodal_u(i, 1) << ','
                       << result.exact_u(i, 0) << ','
                       << result.exact_u(i, 1) << ','
                       << (result.nodal_u(i, 0) - result.exact_u(i, 0)) << ','
                       << (result.nodal_u(i, 1) - result.exact_u(i, 1)) << '\n';
        }
    }

    const std::filesystem::path nodes_file = output_dir / (base_name + "_nodes.csv");
    std::ofstream nodes(nodes_file);
    if (!nodes)
    {
        throw std::runtime_error("Failed to open nodes csv: " + nodes_file.string());
    }
    nodes << std::setprecision(16);
    nodes << "node,x,y\n";
    for (int i = 0; i < mesh.nodes.rows(); ++i)
    {
        nodes << i + 1 << ',' << mesh.nodes(i, 0) << ',' << mesh.nodes(i, 1) << '\n';
    }

    const std::filesystem::path elements_file = output_dir / (base_name + "_elements.csv");
    std::ofstream elements(elements_file);
    if (!elements)
    {
        throw std::runtime_error("Failed to open elements csv: " + elements_file.string());
    }
    elements << "element,n1,n2,n3,n4,n5,n6\n";
    for (std::size_t i = 0; i < mesh.elements.size(); ++i)
    {
        elements << i + 1;
        for (int a = 0; a < 6; ++a)
        {
            elements << ',' << mesh.elements[i][static_cast<std::size_t>(a)] + 1;
        }
        elements << '\n';
    }

    const std::filesystem::path disp_file = output_dir / (base_name + "_displacement.csv");
    std::ofstream disp(disp_file);
    if (!disp)
    {
        throw std::runtime_error("Failed to open displacement csv: " + disp_file.string());
    }
    disp << std::setprecision(16);
    disp << "node,ux,uy\n";
    for (int i = 0; i < result.nodal_u.rows(); ++i)
    {
        disp << i + 1 << ',' << result.nodal_u(i, 0) << ',' << result.nodal_u(i, 1) << '\n';
    }

    if (has_exact)
    {
        const std::filesystem::path exact_file = output_dir / (base_name + "_exact.csv");
        std::ofstream exact(exact_file);
        if (!exact)
        {
            throw std::runtime_error("Failed to open exact csv: " + exact_file.string());
        }
        exact << std::setprecision(16);
        exact << "node,ux,uy\n";
        for (int i = 0; i < result.exact_u.rows(); ++i)
        {
            exact << i + 1 << ',' << result.exact_u(i, 0) << ',' << result.exact_u(i, 1) << '\n';
        }
    }

    const std::filesystem::path initial_shape_file = output_dir / (base_name + "_initial_shape.csv");
    std::ofstream initial_shape(initial_shape_file);
    if (!initial_shape)
    {
        throw std::runtime_error("Failed to open initial shape csv: " + initial_shape_file.string());
    }
    initial_shape << std::setprecision(16);
    initial_shape << "node,x,y\n";
    for (int i = 0; i < mesh.nodes.rows(); ++i)
    {
        initial_shape << i + 1 << ',' << mesh.nodes(i, 0) << ',' << mesh.nodes(i, 1) << '\n';
    }

    const std::filesystem::path deformed_shape_file = output_dir / (base_name + "_deformed_shape.csv");
    std::ofstream deformed_shape(deformed_shape_file);
    if (!deformed_shape)
    {
        throw std::runtime_error("Failed to open deformed shape csv: " + deformed_shape_file.string());
    }
    deformed_shape << std::setprecision(16);
    deformed_shape << "node,x,y\n";
    for (int i = 0; i < mesh.nodes.rows(); ++i)
    {
        deformed_shape << i + 1 << ','
                       << mesh.nodes(i, 0) + result.nodal_u(i, 0) << ','
                       << mesh.nodes(i, 1) + result.nodal_u(i, 1) << '\n';
    }
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

        int num_els = model_parameters::default_num_els[0];
        const auto default_material = model_parameters::nonlinear_material_for("nonlinear_patch");
        Eigen::Vector2d material(default_material.mu, default_material.lambda_like);
        nonlinear_esfem_t6::SolverOptions options;
        options.nstep = model_parameters::default_nstep;
        options.maxiter = model_parameters::default_maxiter;
        options.tolerance = model_parameters::default_tolerance;
        options.residual_tolerance = model_parameters::default_residual_tolerance;
        options.line_search_max_backtracks = model_parameters::default_line_search_max_backtracks;
        options.line_search_reduction = model_parameters::default_line_search_reduction;
        options.line_search_min_alpha = model_parameters::default_line_search_min_alpha;
        options.adaptive_load_stepping = model_parameters::default_adaptive_load_stepping;
        options.max_step_cuts = model_parameters::default_max_step_cuts;
        options.linear_solver = parse_solver_type(model_parameters::default_structural_linear_solver);
        options.iterative_maxiter = model_parameters::default_structural_iterative_maxiter;
        options.iterative_tolerance = model_parameters::default_structural_iterative_tolerance;
        options.compute_exact_solution = model_parameters::default_compute_exact_solution;
        std::filesystem::path output_dir = std::filesystem::path("res/nonlinear_patch/csfem");
        OutputOptions output_options;
        std::optional<std::filesystem::path> benchmark_log;

        for (std::size_t i = 0; i < args.size(); ++i)
        {
            const std::string &arg = args[i];
            if (arg == "nonlinear_patch")
            {
                continue;
            }
            if (!arg.empty() && arg[0] == '-')
            {
                auto require_value = [&](const char *name) -> const std::string & {
                    if (i + 1 >= args.size())
                    {
                        throw std::invalid_argument(std::string("Missing value for ") + name);
                    }
                    ++i;
                    return args[i];
                };

                if (arg == "--num-els")
                {
                    num_els = std::stoi(require_value("--num-els"));
                }
                else if (arg == "--mu")
                {
                    material(0) = std::stod(require_value("--mu"));
                }
                else if (arg == "--lambda_like")
                {
                    material(1) = std::stod(require_value("--lambda_like"));
                }
                else if (arg == "--nstep")
                {
                    options.nstep = std::stoi(require_value("--nstep"));
                }
                else if (arg == "--maxiter")
                {
                    options.maxiter = std::stoi(require_value("--maxiter"));
                }
                else if (arg == "--tol")
                {
                    options.tolerance = std::stod(require_value("--tol"));
                }
                else if (arg == "--solver")
                {
                    options.linear_solver = parse_solver_type(require_value("--solver"));
                }
                else if (arg == "--iter-maxiter")
                {
                    options.iterative_maxiter = std::stoi(require_value("--iter-maxiter"));
                }
                else if (arg == "--iter-tol")
                {
                    options.iterative_tolerance = std::stod(require_value("--iter-tol"));
                }
                else if (arg == "--output-dir")
                {
                    output_dir = require_value("--output-dir");
                }
                else if (arg == "--no-exact")
                {
                    options.compute_exact_solution = false;
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
                else
                {
                    throw std::invalid_argument("Unsupported option: " + arg);
                }
                continue;
            }
            throw std::invalid_argument("Unsupported positional argument: " + arg);
        }

        const auto total_start = std::chrono::steady_clock::now();
        const auto problem =
            nonlinear_esfem_t6::NonlinearEsfemT6Problem::make_not_so_simple_shear(Eigen::Vector2i::Constant(num_els), material);
        const nonlinear_esfem_t6::NonlinearCsfemT6Solver solver;
        const auto solve_start = std::chrono::steady_clock::now();
        const nonlinear_esfem_t6::Result result = solver.solve(problem, options);
        const auto solve_end = std::chrono::steady_clock::now();

        const std::string base_name = csfem_base_name("nonlinear_patch", num_els, num_els);
        if (output_options.write_vtu)
        {
            nonlinear_esfem_t6::write_vtu(output_dir / (base_name + ".vtu"), problem.data().mesh, result.nodal_u, result.exact_u);
        }
        write_matlab_files(output_dir, base_name, problem.data().mesh, result, output_options);
        const auto total_end = std::chrono::steady_clock::now();

        if (output_options.write_benchmark)
        {
            benchmark::Record record;
            record.timestamp_utc = benchmark::utc_timestamp_now();
            record.problem_type = "nonlinear_patch";
            record.method = "csfem";
            record.num_els_x = num_els;
            record.num_els_y = num_els;
            record.nstep = options.nstep;
            record.maxiter = options.maxiter;
            record.solver = solver_name(options.linear_solver);
            record.compute_exact_solution = options.compute_exact_solution;
            record.write_vtu = output_options.write_vtu;
            record.write_postprocess = output_options.write_postprocess;
            record.relative_error = result.relative_error;
            record.strain_energy = result.strain_energy;
            record.solve_wall_time_sec = benchmark::seconds_between(solve_start, solve_end);
            record.total_wall_time_sec = benchmark::seconds_between(total_start, total_end);
            record.peak_rss_mb = benchmark::peak_rss_mb();
            const std::filesystem::path log_path =
                benchmark_log.value_or(output_dir / (base_name + "_benchmark.csv"));
            benchmark::append_csv(log_path, record);
        }

        std::cout << "--------------------------------------------------\n";
        std::cout << std::fixed << std::setprecision(6)
                  << "Final Timing [s] solve="
                  << benchmark::seconds_between(solve_start, solve_end)
                  << " total="
                  << benchmark::seconds_between(total_start, total_end)
                  << '\n';
        std::cout << std::setprecision(6)
                  << "nonlinear_patch / csfem: rel_error=" << result.relative_error
                  << ", strain_energy=" << result.strain_energy
                  << ", output=" << output_dir << '\n';
        return 0;
    }
    catch (const std::exception &ex)
    {
        std::cerr << "Error: " << ex.what() << '\n';
        return 1;
    }
}
