#include "structural_solver.hpp"
#include "model_parameters.hpp"
#include "sfem_shared.hpp"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/LU>
#include <Eigen/SparseLU>
#include <Eigen/SparseQR>
#ifdef USE_EIGEN_UMFPACK
#if __has_include(<Eigen/UmfPackSupport>)
#include <Eigen/UmfPackSupport>
#elif __has_include(<unsupported/Eigen/UmfPackSupport>)
#include <unsupported/Eigen/UmfPackSupport>
#else
#error "USE_EIGEN_UMFPACK is set, but Eigen UmfPackSupport header was not found"
#endif
#endif

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <fstream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace sfem
{
    namespace
    {

        using Matrix23 = Eigen::Matrix<double, 2, 3>;
        using Matrix33 = Eigen::Matrix3d;
        using Matrix34 = Eigen::Matrix<double, 3, 4>;
        using Matrix44 = Eigen::Matrix4d;
        using Matrix46 = Eigen::Matrix<double, 4, 6>;
        using Matrix66 = Eigen::Matrix<double, 6, 6>;
        using SparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
        using Triplet = Eigen::Triplet<double>;
        using Vector6 = Eigen::Matrix<double, 6, 1>;
        namespace shared_sfem = nonlinear_esfem_t6::shared;

        bool is_acceptable_solution(const SparseMatrix &a,
                                    const Eigen::VectorXd &b,
                                    const Eigen::VectorXd &x)
        {
            if (!x.allFinite())
            {
                return false;
            }
            const Eigen::VectorXd residual = a * x - b;
            if (!residual.allFinite())
            {
                return false;
            }
            const double scale = std::max(1.0, b.norm());
            return residual.norm() <= 1.0e-8 * scale;
        }

        bool is_nonlinear_scenario(Scenario scenario)
        {
            return scenario != Scenario::LinearPatch;
        }

        nonlinear_esfem_t6::LinearSolverType to_shared_solver_type(LinearSolverType solver)
        {
            switch (solver)
            {
            case LinearSolverType::SparseLU:
                return nonlinear_esfem_t6::LinearSolverType::SparseLU;
            case LinearSolverType::BiCGSTAB:
                return nonlinear_esfem_t6::LinearSolverType::BiCGSTAB;
            case LinearSolverType::UmfPack:
                return nonlinear_esfem_t6::LinearSolverType::UmfPack;
            }
            throw std::invalid_argument("Unsupported linear solver type");
        }

        nonlinear_esfem_t6::Problem to_shared_problem(const Model &model)
        {
            nonlinear_esfem_t6::Problem problem;
            problem.mesh.nodes = model.mesh.nodes;
            problem.mesh.elements = model.mesh.elements;
            problem.bc.dofs = model.bc.dofs;
            problem.bc.values = model.bc.values;
            problem.force = model.force;
            problem.problem_type = scenario_name(model.scenario);
            problem.material = model.nonlinear_material;
            problem.num_els = model.num_els;
            problem.edge_boundary_ng = model_parameters::esfem_edge_boundary_ng;
            problem.edge_tri_quad_order = model_parameters::esfem_edge_tri_quad_order;
            problem.edge_quad_quad_order = model_parameters::esfem_edge_quad_quad_order;
            problem.edge_reg_param = model_parameters::esfem_edge_reg_param;
            return problem;
        }

        nonlinear_esfem_t6::SolverOptions to_shared_options(const NewtonOptions &options)
        {
            nonlinear_esfem_t6::SolverOptions shared_options;
            shared_options.nstep = options.nstep;
            shared_options.maxiter = options.maxiter;
            shared_options.tolerance = options.tolerance;
            shared_options.residual_tolerance = options.residual_tolerance;
            shared_options.linear_solver = to_shared_solver_type(options.linear_solver);
            shared_options.iterative_maxiter = options.iterative_maxiter;
            shared_options.iterative_tolerance = options.iterative_tolerance;
            shared_options.line_search_max_backtracks = options.line_search_max_backtracks;
            shared_options.line_search_reduction = options.line_search_reduction;
            shared_options.line_search_min_alpha = options.line_search_min_alpha;
            shared_options.adaptive_load_stepping = options.adaptive_load_stepping;
            shared_options.max_step_cuts = options.max_step_cuts;
            shared_options.allow_step_growth = options.allow_step_growth;
            shared_options.aggressive_stagnation_control = options.aggressive_stagnation_control;
            shared_options.compute_exact_solution = false;
            shared_options.debug_csfem_bending = options.debug_csfem_bending;
            shared_options.debug_output_dir = options.debug_output_dir;
            return shared_options;
        }

        struct QuadratureRule
        {
            Eigen::VectorXd weights;
            Eigen::MatrixXd points;
        };

        struct NonlinearMaterial
        {
            Matrix33 cmat = Matrix33::Zero();
            Matrix44 smat = Matrix44::Zero();
            double w0 = 0.0;
        };

        struct SupportData
        {
            std::vector<int> nodes;
            Eigen::MatrixXd gcoord;
            Eigen::MatrixXd fx;
            Eigen::MatrixXd fy;
            std::vector<double> weights;
            Eigen::MatrixXd points;
            Eigen::MatrixXd dx;
            Eigen::MatrixXd dy;
            std::vector<int> edof;
            double sub_area = 0.0;
        };

        void write_csfem_bending_debug_dump(const Model &model,
                                            const NewtonOptions &options,
                                            const SupportData &data,
                                            const Eigen::VectorXd &uu,
                                            const Eigen::VectorXd &rhs,
                                            const SparseMatrix &k,
                                            int step,
                                            int iteration,
                                            double scale)
        {
            if (!options.debug_csfem_bending || options.debug_output_dir.empty())
            {
                return;
            }

            std::filesystem::create_directories(options.debug_output_dir);
            const std::filesystem::path file =
                options.debug_output_dir / ("CSFEM_bending_block_debug_step" + std::to_string(step) +
                                            "_iter" + std::to_string(iteration) + ".m");
            std::ofstream out(file);
            if (!out)
            {
                return;
            }

            out << std::setprecision(16);
            out << "debug.step = " << step << ";\n";
            out << "debug.iteration = " << iteration << ";\n";
            out << "debug.scale = " << scale << ";\n";
            out << "debug.uu_norm = " << uu.norm() << ";\n";
            out << "debug.rhs_norm = " << rhs.norm() << ";\n";
            out << "debug.k_diag_min = " << k.diagonal().minCoeff() << ";\n";
            out << "debug.k_diag_max = " << k.diagonal().maxCoeff() << ";\n";
            out << "debug.num_nodes = " << model.mesh.nodes.rows() << ";\n";
            out << "debug.num_elements = " << model.mesh.elements.size() << ";\n";

            out << "debug.first_element_nodes = [";
            for (std::size_t i = 0; i < data.nodes.size(); ++i)
            {
                out << data.nodes[i] + 1;
                out << (i + 1 == data.nodes.size() ? "" : ", ");
            }
            out << "];\n";

            out << "debug.first_element_coords = [\n";
            for (int i = 0; i < data.gcoord.rows(); ++i)
            {
                out << "  " << data.gcoord(i, 0) << ", " << data.gcoord(i, 1) << ";\n";
            }
            out << "];\n";

            out << "debug.first_element_edof = [";
            for (std::size_t i = 0; i < data.edof.size(); ++i)
            {
                out << data.edof[i] + 1;
                out << (i + 1 == data.edof.size() ? "" : ", ");
            }
            out << "];\n";

            out << "debug.first_element_uu = [\n";
            for (int i = 0; i < static_cast<int>(data.nodes.size()); ++i)
            {
                out << "  " << uu(2 * data.nodes[static_cast<std::size_t>(i)]) << ", "
                    << uu(2 * data.nodes[static_cast<std::size_t>(i)] + 1) << ";\n";
            }
            out << "];\n";

            auto write_matrix = [&](const char *name, const Eigen::MatrixXd &mat)
            {
                out << "debug." << name << " = [\n";
                for (int r = 0; r < mat.rows(); ++r)
                {
                    out << "  ";
                    for (int c = 0; c < mat.cols(); ++c)
                    {
                        out << mat(r, c);
                        if (c + 1 < mat.cols())
                        {
                            out << ", ";
                        }
                    }
                    out << ";\n";
                }
                out << "];\n";
            };

            write_matrix("first_element_fx", data.fx);
            write_matrix("first_element_fy", data.fy);
            write_matrix("first_element_points", data.points);
            write_matrix("first_element_dx", data.dx);
            write_matrix("first_element_dy", data.dy);

            out << "debug.first_element_weights = [";
            for (std::size_t i = 0; i < data.weights.size(); ++i)
            {
                out << data.weights[i];
                out << (i + 1 == data.weights.size() ? "" : ", ");
            }
            out << "];\n";
            out << "debug.first_element_sub_area = " << data.sub_area << ";\n";
        }

        void write_csfem_bending_failure_dump(const Model &model,
                                              const NewtonOptions &options,
                                              const SupportData &data,
                                              const Eigen::MatrixXd &wku,
                                              const Eigen::MatrixXd &dndx,
                                              const Eigen::MatrixXd &bmat,
                                              const Eigen::MatrixXd &bgeo,
                                              const Eigen::Matrix3d &fmat,
                                              const NonlinearMaterial &material,
                                              int step,
                                              int iteration,
                                              int element_index,
                                              int gauss_index,
                                              double scale,
                                              double weight,
                                              const char *reason)
        {
            if (!options.debug_csfem_bending || options.debug_output_dir.empty())
            {
                return;
            }

            std::filesystem::create_directories(options.debug_output_dir);
            const std::filesystem::path file =
                options.debug_output_dir / ("CSFEM_bending_block_failure_step" + std::to_string(step) +
                                            "_iter" + std::to_string(iteration) +
                                            "_elem" + std::to_string(element_index + 1) +
                                            "_gp" + std::to_string(gauss_index + 1) + ".m");
            std::ofstream out(file);
            if (!out)
            {
                return;
            }

            out << std::setprecision(16);
            out << "failure.reason = '" << reason << "';\n";
            out << "failure.step = " << step << ";\n";
            out << "failure.iteration = " << iteration << ";\n";
            out << "failure.element_index = " << element_index + 1 << ";\n";
            out << "failure.gauss_index = " << gauss_index + 1 << ";\n";
            out << "failure.scale = " << scale << ";\n";
            out << "failure.weight = " << weight << ";\n";

            auto write_matrix = [&](const char *name, const Eigen::MatrixXd &mat)
            {
                out << "failure." << name << " = [\n";
                for (int r = 0; r < mat.rows(); ++r)
                {
                    out << "  ";
                    for (int c = 0; c < mat.cols(); ++c)
                    {
                        out << mat(r, c);
                        if (c + 1 < mat.cols())
                        {
                            out << ", ";
                        }
                    }
                    out << ";\n";
                }
                out << "];\n";
            };

            out << "failure.nodes = [";
            for (std::size_t i = 0; i < data.nodes.size(); ++i)
            {
                out << data.nodes[i] + 1;
                out << (i + 1 == data.nodes.size() ? "" : ", ");
            }
            out << "];\n";

            write_matrix("gcoord", data.gcoord);
            write_matrix("wku", wku);
            write_matrix("dndx", dndx);
            write_matrix("bmat", bmat);
            write_matrix("bgeo", bgeo);
            write_matrix("fmat", fmat);
            write_matrix("cmat", material.cmat);
            write_matrix("smat", material.smat);
        }

        Eigen::VectorXd solve_with_fresh_sparse_lu(SparseMatrix a, const Eigen::VectorXd &b);
        Eigen::VectorXd solve_with_umfpack_fallback(const SparseMatrix &a, const Eigen::VectorXd &b);
        Eigen::VectorXd solve_with_sparse_qr_fallback(const SparseMatrix &a, const Eigen::VectorXd &b);
        Eigen::VectorXd solve_with_iterative_fallback(const SparseMatrix &a,
                                                      const Eigen::VectorXd &b,
                                                      int maxiter,
                                                      double tolerance);
        Eigen::VectorXd solve_with_dense_fallback(const SparseMatrix &a, const Eigen::VectorXd &b);
        int regularize_near_zero_rows(SparseMatrix &a, Eigen::VectorXd &b);

        class SparseDirectSolverCache
        {
        public:
            Eigen::VectorXd solve(SparseMatrix a, const Eigen::VectorXd &b, const NewtonOptions &options)
            {
                a.makeCompressed();
                if (options.linear_solver == LinearSolverType::UmfPack)
                {
#ifdef USE_EIGEN_UMFPACK
                    Eigen::UmfPackLU<SparseMatrix> umfpack;
                    umfpack.compute(a);
                    if (umfpack.info() == Eigen::Success)
                    {
                        const Eigen::VectorXd x = umfpack.solve(b);
                        if (umfpack.info() == Eigen::Success && x.allFinite())
                        {
                            return x;
                        }
                    }
#else
                    throw std::runtime_error("UMFPACK solver requested, but Eigen UmfPackSupport is unavailable in this build");
#endif
                }
                if (options.linear_solver == LinearSolverType::BiCGSTAB)
                {
                    Eigen::BiCGSTAB<SparseMatrix, Eigen::IncompleteLUT<double>> iterative_solver;
                    iterative_solver.setMaxIterations(options.iterative_maxiter);
                    iterative_solver.setTolerance(options.iterative_tolerance);
                    iterative_solver.compute(a);
                    if (iterative_solver.info() == Eigen::Success)
                    {
                        const Eigen::VectorXd x = iterative_solver.solve(b);
                        if (iterative_solver.info() == Eigen::Success && x.allFinite())
                        {
                            return x;
                        }
                    }

                    return solve_with_iterative_fallback(a,
                                                         b,
                                                         std::max(options.iterative_maxiter, 6000),
                                                         std::min(options.iterative_tolerance, 1.0e-11));
                }

                analyze_if_needed(a);
                sparse_lu_.factorize(a);
                if (sparse_lu_.info() != Eigen::Success)
                {
                    analyzed_ = false;
                    analyze_if_needed(a);
                    sparse_lu_.factorize(a);
                }
                if (sparse_lu_.info() == Eigen::Success)
                {
                    return sparse_lu_.solve(b);
                }

                a.diagonal().array() += 1.0e-10;
                a.makeCompressed();
                analyze_if_needed(a);
                sparse_lu_.factorize(a);
                if (sparse_lu_.info() != Eigen::Success)
                {
                    analyzed_ = false;
                    analyze_if_needed(a);
                    sparse_lu_.factorize(a);
                }
                if (sparse_lu_.info() != Eigen::Success)
                {
                    analyzed_ = false;
                    try
                    {
                        return solve_with_umfpack_fallback(a, b);
                    }
                    catch (const std::runtime_error &)
                    {
                    }
                    try
                    {
                        return solve_with_fresh_sparse_lu(a, b);
                    }
                    catch (const std::runtime_error &)
                    {
                        try
                        {
                            return solve_with_sparse_qr_fallback(a, b);
                        }
                        catch (const std::runtime_error &)
                        {
                            try
                            {
                                return solve_with_iterative_fallback(a,
                                                                     b,
                                                                     std::max(options.iterative_maxiter, 4000),
                                                                     std::min(options.iterative_tolerance, 1.0e-10));
                            }
                            catch (const std::runtime_error &)
                            {
                                return solve_with_dense_fallback(a, b);
                            }
                        }
                    }
                }
                const Eigen::VectorXd x = sparse_lu_.solve(b);
                if (is_acceptable_solution(a, b, x))
                {
                    return x;
                }
                try
                {
                    return solve_with_umfpack_fallback(a, b);
                }
                catch (const std::runtime_error &)
                {
                }
                return x;
            }

        private:
            void analyze_if_needed(const SparseMatrix &a)
            {
                if (same_pattern(a))
                {
                    return;
                }
                sparse_lu_.analyzePattern(a);
                rows_ = a.rows();
                cols_ = a.cols();
                outer_index_.assign(a.outerIndexPtr(), a.outerIndexPtr() + a.outerSize() + 1);
                inner_index_.assign(a.innerIndexPtr(), a.innerIndexPtr() + a.nonZeros());
                analyzed_ = true;
            }

            bool same_pattern(const SparseMatrix &a) const
            {
                if (!analyzed_ || rows_ != a.rows() || cols_ != a.cols())
                {
                    return false;
                }
                if (static_cast<int>(outer_index_.size()) != a.outerSize() + 1 ||
                    static_cast<int>(inner_index_.size()) != a.nonZeros())
                {
                    return false;
                }
                for (int i = 0; i < a.outerSize() + 1; ++i)
                {
                    if (outer_index_[static_cast<std::size_t>(i)] != a.outerIndexPtr()[i])
                    {
                        return false;
                    }
                }
                for (int i = 0; i < a.nonZeros(); ++i)
                {
                    if (inner_index_[static_cast<std::size_t>(i)] != a.innerIndexPtr()[i])
                    {
                        return false;
                    }
                }
                return true;
            }

            Eigen::SparseLU<SparseMatrix> sparse_lu_;
            int rows_ = -1;
            int cols_ = -1;
            std::vector<int> outer_index_;
            std::vector<int> inner_index_;
            bool analyzed_ = false;
        };

        struct TargetEdge
        {
            int n1 = 0;
            int n2 = 0;
            int elem1 = 0;
            int elem2 = -1;
            int local_edge1 = 0;
            int local_edge2 = -1;
        };

        Eigen::VectorXd solve_with_fresh_sparse_lu(SparseMatrix a, const Eigen::VectorXd &b);

        Eigen::VectorXd solve_with_fresh_sparse_lu(SparseMatrix a, const Eigen::VectorXd &b)
        {
            constexpr double kShifts[] = {0.0, 1.0e-12, 1.0e-10, 1.0e-8, 1.0e-6};
            for (const double shift : kShifts)
            {
                SparseMatrix trial = a;
                if (shift > 0.0)
                {
                    trial.diagonal().array() += shift;
                }
                trial.makeCompressed();
                Eigen::SparseLU<SparseMatrix> lu;
                lu.analyzePattern(trial);
                lu.factorize(trial);
                if (lu.info() == Eigen::Success)
                {
                    const Eigen::VectorXd x = lu.solve(b);
                    if (lu.info() == Eigen::Success && x.allFinite())
                    {
                        return x;
                    }
                }
            }
            throw std::runtime_error("Failed to factorize sparse tangent matrix");
        }

        Eigen::VectorXd solve_with_umfpack_fallback(const SparseMatrix &a, const Eigen::VectorXd &b)
        {
#ifdef USE_EIGEN_UMFPACK
            Eigen::UmfPackLU<SparseMatrix> umfpack;
            umfpack.compute(a);
            if (umfpack.info() == Eigen::Success)
            {
                const Eigen::VectorXd x = umfpack.solve(b);
                if (umfpack.info() == Eigen::Success && x.allFinite())
                {
                    return x;
                }
            }
#endif
            throw std::runtime_error("Failed to factorize sparse tangent matrix");
        }

        Eigen::VectorXd solve_with_sparse_qr_fallback(const SparseMatrix &a, const Eigen::VectorXd &b)
        {
            constexpr double kShifts[] = {0.0, 1.0e-12, 1.0e-10, 1.0e-8, 1.0e-6, 1.0e-4, 1.0e-2};
            Eigen::VectorXd best_x;
            double best_residual = std::numeric_limits<double>::infinity();

            for (const double shift : kShifts)
            {
                SparseMatrix trial = a;
                if (shift > 0.0)
                {
                    trial.diagonal().array() += shift;
                }
                trial.makeCompressed();

                Eigen::SparseQR<SparseMatrix, Eigen::COLAMDOrdering<int>> qr;
                qr.compute(trial);
                if (qr.info() != Eigen::Success)
                {
                    continue;
                }

                const Eigen::VectorXd x = qr.solve(b);
                if (qr.info() != Eigen::Success || !x.allFinite())
                {
                    continue;
                }

                const Eigen::VectorXd residual = a * x - b;
                if (!residual.allFinite())
                {
                    continue;
                }

                const double norm = residual.norm();
                if (norm < best_residual)
                {
                    best_residual = norm;
                    best_x = x;
                }

                if (is_acceptable_solution(a, b, x))
                {
                    return x;
                }
            }

            if (best_x.size() != 0)
            {
                return best_x;
            }

            throw std::runtime_error("Failed to factorize sparse tangent matrix");
        }

        Eigen::VectorXd solve_with_iterative_fallback(const SparseMatrix &a,
                                                      const Eigen::VectorXd &b,
                                                      int maxiter,
                                                      double tolerance)
        {
            auto try_bicgstab_ilut = [&](SparseMatrix mat,
                                         double diag_shift,
                                         int local_maxiter,
                                         double local_tolerance) -> Eigen::VectorXd
            {
                if (diag_shift > 0.0)
                {
                    mat.diagonal().array() += diag_shift;
                    mat.makeCompressed();
                }
                Eigen::BiCGSTAB<SparseMatrix, Eigen::IncompleteLUT<double>> iterative_solver;
                iterative_solver.setMaxIterations(local_maxiter);
                iterative_solver.setTolerance(local_tolerance);
                iterative_solver.preconditioner().setDroptol(1.0e-3);
                iterative_solver.preconditioner().setFillfactor(20);
                iterative_solver.compute(mat);
                if (iterative_solver.info() != Eigen::Success)
                {
                    return Eigen::VectorXd();
                }
                const Eigen::VectorXd x = iterative_solver.solve(b);
                if (iterative_solver.info() != Eigen::Success || !x.allFinite())
                {
                    return Eigen::VectorXd();
                }
                return x;
            };

            auto try_bicgstab_diag = [&](SparseMatrix mat,
                                         double diag_shift,
                                         int local_maxiter,
                                         double local_tolerance) -> Eigen::VectorXd
            {
                if (diag_shift > 0.0)
                {
                    mat.diagonal().array() += diag_shift;
                    mat.makeCompressed();
                }
                Eigen::BiCGSTAB<SparseMatrix, Eigen::DiagonalPreconditioner<double>> iterative_solver;
                iterative_solver.setMaxIterations(local_maxiter);
                iterative_solver.setTolerance(local_tolerance);
                iterative_solver.compute(mat);
                if (iterative_solver.info() != Eigen::Success)
                {
                    return Eigen::VectorXd();
                }
                const Eigen::VectorXd x = iterative_solver.solve(b);
                if (iterative_solver.info() != Eigen::Success || !x.allFinite())
                {
                    return Eigen::VectorXd();
                }
                return x;
            };

            Eigen::VectorXd x = try_bicgstab_ilut(a, 0.0, maxiter, tolerance);
            if (x.size() != 0)
            {
                return x;
            }

            x = try_bicgstab_ilut(a, 1.0e-10, std::max(maxiter, 4000), std::min(tolerance, 1.0e-10));
            if (x.size() != 0)
            {
                return x;
            }

            x = try_bicgstab_ilut(a, 1.0e-6, std::max(maxiter, 8000), std::min(tolerance, 1.0e-10));
            if (x.size() != 0)
            {
                return x;
            }

            x = try_bicgstab_diag(a, 1.0e-8, std::max(maxiter, 6000), std::min(tolerance, 1.0e-11));
            if (x.size() != 0)
            {
                return x;
            }

            x = try_bicgstab_diag(a, 1.0e-4, std::max(maxiter, 10000), std::min(tolerance, 1.0e-10));
            if (x.size() != 0)
            {
                return x;
            }

            throw std::runtime_error("Failed to factorize sparse tangent matrix");
        }

        Eigen::VectorXd solve_with_dense_fallback(const SparseMatrix &a, const Eigen::VectorXd &b)
        {
            constexpr int kDenseFallbackMaxDofs = 6000;
            if (a.rows() > kDenseFallbackMaxDofs || a.cols() > kDenseFallbackMaxDofs)
            {
                throw std::runtime_error("Failed to factorize sparse tangent matrix");
            }

            Eigen::MatrixXd dense = Eigen::MatrixXd(a);
            Eigen::VectorXd best_x;
            double best_residual = std::numeric_limits<double>::infinity();
            auto consider_candidate = [&](const Eigen::VectorXd &x)
            {
                if (!x.allFinite())
                {
                    return false;
                }
                const Eigen::VectorXd residual = a * x - b;
                if (!residual.allFinite())
                {
                    return false;
                }
                const double norm = residual.norm();
                if (norm < best_residual)
                {
                    best_residual = norm;
                    best_x = x;
                }
                return is_acceptable_solution(a, b, x);
            };

            {
                Eigen::FullPivLU<Eigen::MatrixXd> lu(dense);
                if (lu.isInvertible())
                {
                    const Eigen::VectorXd x = lu.solve(b);
                    if (consider_candidate(x))
                    {
                        return x;
                    }
                }
            }

            constexpr double kShifts[] = {1.0e-12, 1.0e-10, 1.0e-8, 1.0e-6};
            for (const double shift : kShifts)
            {
                Eigen::MatrixXd shifted = dense;
                shifted.diagonal().array() += shift;
                Eigen::FullPivLU<Eigen::MatrixXd> lu(shifted);
                if (!lu.isInvertible())
                {
                    continue;
                }
                const Eigen::VectorXd x = lu.solve(b);
                if (consider_candidate(x))
                {
                    return x;
                }
            }

            if (consider_candidate(dense.completeOrthogonalDecomposition().solve(b)))
            {
                return best_x;
            }
            constexpr double kRelaxedShifts[] = {1.0e-4, 1.0e-3, 1.0e-2};
            for (const double shift : kRelaxedShifts)
            {
                Eigen::MatrixXd shifted = dense;
                shifted.diagonal().array() += shift;
                if (consider_candidate(shifted.completeOrthogonalDecomposition().solve(b)))
                {
                    return best_x;
                }
            }
            if (best_x.size() != 0)
            {
                return best_x;
            }
            throw std::runtime_error("Failed to factorize sparse tangent matrix");
        }

        int regularize_near_zero_rows(SparseMatrix &a, Eigen::VectorXd &b)
        {
            constexpr double kRowTol = 1.0e-14;
            std::vector<int> pinned_rows;
            pinned_rows.reserve(32);

            for (int row = 0; row < a.rows(); ++row)
            {
                double row_abs_sum = 0.0;
                for (SparseMatrix::InnerIterator it(a, row); it; ++it)
                {
                    row_abs_sum += std::abs(it.value());
                }
                if (row_abs_sum <= kRowTol)
                {
                    pinned_rows.push_back(row);
                }
            }

            for (const int row : pinned_rows)
            {
                a.coeffRef(row, row) = 1.0;
                b(row) = 0.0;
            }

            if (!pinned_rows.empty())
            {
                a.makeCompressed();
            }
            return static_cast<int>(pinned_rows.size());
        }

        double triangle_area_from_vertices(const Eigen::MatrixXd &coords)
        {
            return 0.5 * std::abs((coords(1, 0) - coords(0, 0)) * (coords(2, 1) - coords(0, 1)) -
                                  (coords(2, 0) - coords(0, 0)) * (coords(1, 1) - coords(0, 1)));
        }

        double polygon_area(const Eigen::MatrixXd &coords)
        {
            double area = 0.0;
            for (int i = 0; i < coords.rows(); ++i)
            {
                const int j = (i + 1) % coords.rows();
                area += coords(i, 0) * coords(j, 1) - coords(j, 0) * coords(i, 1);
            }
            return 0.5 * std::abs(area);
        }

        QuadratureRule gauss_line_rule(int order)
        {
            QuadratureRule rule;
            if (order == 2)
            {
                rule.weights.resize(2);
                rule.points.resize(2, 1);
                const double a = 0.5773502691896257;
                rule.weights << 1.0, 1.0;
                rule.points << a, -a;
            }
            else if (order == 3)
            {
                rule.weights.resize(3);
                rule.points.resize(3, 1);
                const double a = 0.7745966692414834;
                rule.weights << 0.5555555555555556, 0.5555555555555556, 0.8888888888888888;
                rule.points << a, -a, 0.0;
            }
            else
            {
                throw std::invalid_argument("Unsupported line quadrature order");
            }
            return rule;
        }

        QuadratureRule gauss_quad_rule(int order)
        {
            const QuadratureRule line = gauss_line_rule(order);
            QuadratureRule rule;
            rule.weights.resize(line.weights.size() * line.weights.size());
            rule.points.resize(line.weights.size() * line.weights.size(), 2);
            int k = 0;
            for (int i = 0; i < line.weights.size(); ++i)
            {
                for (int j = 0; j < line.weights.size(); ++j)
                {
                    rule.weights(k) = line.weights(i) * line.weights(j);
                    rule.points.row(k) << line.points(i, 0), line.points(j, 0);
                    ++k;
                }
            }
            return rule;
        }

        QuadratureRule triangle_rule(int order)
        {
            QuadratureRule rule;
            if (order == 2)
            {
                rule.weights = Eigen::Vector3d::Constant(1.0 / 6.0);
                rule.points.resize(3, 2);
                rule.points << 1.0 / 6.0, 1.0 / 6.0,
                    2.0 / 3.0, 1.0 / 6.0,
                    1.0 / 6.0, 2.0 / 3.0;
            }
            else if (order == 3)
            {
                rule.weights.resize(4);
                rule.points.resize(4, 2);
                rule.points << 1.0 / 3.0, 1.0 / 3.0,
                    0.2, 0.2,
                    0.2, 0.6,
                    0.6, 0.2;
                rule.weights << -0.28125, 0.2604166666666667, 0.2604166666666667, 0.2604166666666667;
            }
            else if (order == 4)
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
            }
            else
            {
                throw std::invalid_argument("Unsupported triangular quadrature order");
            }
            return rule;
        }

        void lagrange_basis(const std::string &type, const Eigen::Vector2d &coord, Eigen::VectorXd &n, Eigen::MatrixXd &dndxi)
        {
            const double x = coord(0);
            const double y = coord(1);
            if (type == "T3")
            {
                n.resize(3);
                dndxi.resize(3, 2);
                n << 1.0 - x - y, x, y;
                dndxi << -1.0, -1.0,
                    1.0, 0.0,
                    0.0, 1.0;
            }
            else if (type == "T6")
            {
                n.resize(6);
                dndxi.resize(6, 2);
                n << (1.0 - x - y) * (1.0 - 2.0 * x - 2.0 * y),
                    x * (2.0 * x - 1.0),
                    y * (2.0 * y - 1.0),
                    4.0 * x * (1.0 - x - y),
                    4.0 * x * y,
                    4.0 * y * (1.0 - x - y);
                dndxi.col(0) << 4.0 * (x + y) - 3.0, 4.0 * x - 1.0, 0.0, -4.0 * (2.0 * x + y - 1.0), 4.0 * y, -4.0 * y;
                dndxi.col(1) << 4.0 * (x + y) - 3.0, 0.0, 4.0 * y - 1.0, -4.0 * x, 4.0 * x, -4.0 * (x + 2.0 * y - 1.0);
            }
            else if (type == "Q4")
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
            }
            else if (type == "Q9")
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
            }
            else if (type == "L3")
            {
                n.resize(3);
                dndxi.resize(3, 1);
                n << -0.5 * (1.0 - x) * x, 0.5 * (1.0 + x) * x, 1.0 - x * x;
                dndxi << x - 0.5, x + 0.5, -2.0 * x;
            }
            else
            {
                throw std::invalid_argument("Unsupported element type in lagrange_basis");
            }
        }

        Eigen::Vector2d inverse_mapping(const std::string &type, const Eigen::MatrixXd &nodes, const Eigen::Vector2d &pt)
        {
            Eigen::Vector2d coord = Eigen::Vector2d::Zero();
            for (int iter = 0; iter < 10; ++iter)
            {
                Eigen::VectorXd n;
                Eigen::MatrixXd dndxi;
                lagrange_basis(type, coord, n, dndxi);
                const Eigen::Vector2d xy = nodes.transpose() * n;
                const Eigen::Matrix2d f = nodes.transpose() * dndxi;
                coord -= f.inverse() * (xy - pt);
            }
            return coord;
        }

        Eigen::VectorXd serendipity_shape(const std::string &type, const Eigen::MatrixXd &nodes, const Eigen::Vector2d &pt)
        {
            Eigen::VectorXd n;
            Eigen::MatrixXd dndxi;
            const Eigen::Vector2d xi_eta = inverse_mapping(type, nodes, pt);
            lagrange_basis(type, xi_eta, n, dndxi);
            return n;
        }

        Matrix33 linear_constitutive(double e, double nu)
        {
            Matrix33 c = Matrix33::Zero();
            const double c0 = e / ((1.0 + nu) * (1.0 - 2.0 * nu));
            c(0, 0) = 1.0 - nu;
            c(0, 1) = nu;
            c(1, 0) = nu;
            c(1, 1) = 1.0 - nu;
            c(2, 2) = 0.5 * (1.0 - 2.0 * nu);
            return c0 * c;
        }

        NonlinearMaterial nonlinear_constitutive(const Eigen::Vector2d &param, const Eigen::Matrix3d &fmat)
        {
            Eigen::Matrix3d c = fmat.transpose() * fmat;
            double i3 = c.determinant();
            if (!std::isfinite(i3) || i3 <= 1.0e-12)
            {
                c.diagonal().array() += 1.0e-10;
                i3 = c.determinant();
            }
            const Eigen::Matrix3d inv_c = c.inverse();
            const double i1 = c.trace();
            i3 = std::max(i3, 1.0e-12);

            const double lambda = param(1) - (2.0 / 3.0) * param(0);
            const double mu = param(0) - 0.5 * lambda * std::log(i3);

            const Eigen::Matrix3d ident = Eigen::Matrix3d::Identity();
            const Eigen::Matrix3d s = 0.5 * lambda * std::log(i3) * inv_c + param(0) * (ident - inv_c);

            NonlinearMaterial material;
            material.smat.setZero();
            material.smat.block<2, 2>(0, 0) = s.block<2, 2>(0, 0);
            material.smat.block<2, 2>(2, 2) = s.block<2, 2>(0, 0);

            material.cmat(0, 0) = lambda * inv_c(0, 0) * inv_c(0, 0) + mu * (2.0 * inv_c(0, 0) * inv_c(0, 0));
            material.cmat(0, 1) = lambda * inv_c(0, 0) * inv_c(1, 1) + mu * (2.0 * inv_c(0, 1) * inv_c(0, 1));
            material.cmat(0, 2) = lambda * inv_c(0, 0) * inv_c(0, 1) + mu * (inv_c(0, 0) * inv_c(0, 1) + inv_c(0, 1) * inv_c(0, 0));
            material.cmat(1, 0) = material.cmat(0, 1);
            material.cmat(1, 1) = lambda * inv_c(1, 1) * inv_c(1, 1) + mu * (2.0 * inv_c(1, 1) * inv_c(1, 1));
            material.cmat(1, 2) = lambda * inv_c(1, 1) * inv_c(0, 1) + mu * (inv_c(1, 0) * inv_c(1, 1) + inv_c(1, 1) * inv_c(1, 0));
            material.cmat(2, 0) = material.cmat(0, 2);
            material.cmat(2, 1) = material.cmat(1, 2);
            material.cmat(2, 2) = lambda * inv_c(0, 1) * inv_c(0, 1) + mu * (inv_c(0, 0) * inv_c(1, 1) + inv_c(0, 1) * inv_c(1, 0));

            material.w0 = lambda * std::pow(std::log(i3), 2.0) / 8.0 - 0.5 * param(0) * std::log(i3) + 0.5 * param(0) * (i1 - 3.0);
            return material;
        }

        void nonlinear_bmat(const Eigen::MatrixXd &dndx,
                            const Eigen::MatrixXd &wk_u,
                            Eigen::MatrixXd &bmat,
                            Eigen::MatrixXd &bgeo,
                            Eigen::Matrix3d &fmat)
        {
            const Eigen::Matrix2d dxdx0 = wk_u.transpose() * dndx;
            fmat.setIdentity();
            fmat.block<2, 2>(0, 0) += dxdx0;

            const int nnode = static_cast<int>(dndx.rows());
            const int ncol = 2 * nnode;
            if (bmat.rows() != 3 || bmat.cols() != ncol)
            {
                bmat.resize(3, ncol);
            }
            if (bgeo.rows() != 4 || bgeo.cols() != ncol)
            {
                bgeo.resize(4, ncol);
            }
            bmat.setZero();
            bgeo.setZero();

            for (int a = 0; a < nnode; ++a)
            {
                const double dx = dndx(a, 0);
                const double dy = dndx(a, 1);
                bmat(0, 2 * a) = dx * fmat(0, 0);
                bmat(0, 2 * a + 1) = dx * fmat(1, 0);
                bmat(1, 2 * a) = dy * fmat(0, 1);
                bmat(1, 2 * a + 1) = dy * fmat(1, 1);
                bmat(2, 2 * a) = dy * fmat(0, 0) + dx * fmat(0, 1);
                bmat(2, 2 * a + 1) = dx * fmat(1, 1) + dy * fmat(1, 0);

                bgeo(0, 2 * a) = dx;
                bgeo(1, 2 * a) = dy;
                bgeo(2, 2 * a + 1) = dx;
                bgeo(3, 2 * a + 1) = dy;
            }
        }

        std::vector<double> side_lengths(const Eigen::MatrixXd &xy)
        {
            std::vector<double> side(static_cast<std::size_t>(xy.rows()));
            for (int i = 0; i < xy.rows(); ++i)
            {
                const int j = (i + 1) % xy.rows();
                side[static_cast<std::size_t>(i)] = (xy.row(j) - xy.row(i)).norm();
            }
            return side;
        }

        void outward_normals(const Eigen::MatrixXd &xy, const std::vector<double> &side, std::vector<double> &nx, std::vector<double> &ny)
        {
            nx.resize(static_cast<std::size_t>(xy.rows()));
            ny.resize(static_cast<std::size_t>(xy.rows()));
            for (int i = 0; i < xy.rows(); ++i)
            {
                const int j = (i + 1) % xy.rows();
                nx[static_cast<std::size_t>(i)] = (xy(j, 1) - xy(i, 1)) / side[static_cast<std::size_t>(i)];
                ny[static_cast<std::size_t>(i)] = -(xy(j, 0) - xy(i, 0)) / side[static_cast<std::size_t>(i)];
            }
        }

        std::vector<int> element_dofs(const std::array<int, 6> &element)
        {
            std::vector<int> edof(12);
            for (int a = 0; a < 6; ++a)
            {
                edof[static_cast<std::size_t>(2 * a)] = 2 * element[static_cast<std::size_t>(a)];
                edof[static_cast<std::size_t>(2 * a + 1)] = 2 * element[static_cast<std::size_t>(a)] + 1;
            }
            return edof;
        }

        void scatter_add(std::vector<Triplet> &global, const std::vector<int> &edof, const Eigen::MatrixXd &local)
        {
            for (int i = 0; i < static_cast<int>(edof.size()); ++i)
            {
                for (int j = 0; j < static_cast<int>(edof.size()); ++j)
                {
                    const double value = local(i, j);
                    if (value != 0.0)
                    {
                        global.emplace_back(edof[static_cast<std::size_t>(i)], edof[static_cast<std::size_t>(j)], value);
                    }
                }
            }
        }

        void scatter_add(Eigen::VectorXd &global, const std::vector<int> &edof, const Eigen::VectorXd &local)
        {
            for (int i = 0; i < static_cast<int>(edof.size()); ++i)
            {
                global(edof[static_cast<std::size_t>(i)]) += local(i);
            }
        }

        SparseMatrix make_sparse_matrix(int rows, int cols, const std::vector<Triplet> &triplets)
        {
            SparseMatrix a(rows, cols);
            a.setFromTriplets(triplets.begin(), triplets.end(), [](double lhs, double rhs)
                              { return lhs + rhs; });
            a.makeCompressed();
            return a;
        }

        Eigen::VectorXd solve_system(SparseMatrix a, const Eigen::VectorXd &b, const NewtonOptions &options);

        Eigen::VectorXd solve_system(SparseMatrix a, const Eigen::VectorXd &b)
        {
            NewtonOptions options;
            SparseDirectSolverCache solver_cache;
            return solver_cache.solve(std::move(a), b, options);
        }

        Eigen::VectorXd solve_system(SparseMatrix a, const Eigen::VectorXd &b, const NewtonOptions &options)
        {
            SparseDirectSolverCache solver_cache;
            return solver_cache.solve(std::move(a), b, options);
        }

        void apply_linear_bc(const BoundaryCondition &bc, SparseMatrix &k, Eigen::VectorXd &f)
        {
            std::vector<char> constrained(static_cast<std::size_t>(k.cols()), 0);
            Eigen::VectorXd prescribed = Eigen::VectorXd::Zero(k.cols());
            for (int i = 0; i < static_cast<int>(bc.dofs.size()); ++i)
            {
                const int c = bc.dofs[static_cast<std::size_t>(i)];
                constrained[static_cast<std::size_t>(c)] = 1;
                prescribed(c) = bc.values(i);
            }

            for (int row = 0; row < k.outerSize(); ++row)
            {
                for (SparseMatrix::InnerIterator it(k, row); it; ++it)
                {
                    if (constrained[static_cast<std::size_t>(it.col())])
                    {
                        f(row) -= it.value() * prescribed(it.col());
                        if (row != it.col())
                        {
                            it.valueRef() = 0.0;
                        }
                    }
                }
            }

            for (int i = 0; i < static_cast<int>(bc.dofs.size()); ++i)
            {
                const int c = bc.dofs[static_cast<std::size_t>(i)];
                for (SparseMatrix::InnerIterator it(k, c); it; ++it)
                {
                    if (it.col() != c)
                    {
                        it.valueRef() = 0.0;
                    }
                }
                k.coeffRef(c, c) = 1.0;
                f(c) = prescribed(c);
            }
            k.prune(0.0);
        }

        void apply_nonlinear_bc(const BoundaryCondition &bc, const Eigen::VectorXd &uu, double scale, SparseMatrix &k, Eigen::VectorXd &rhs)
        {
            std::vector<char> constrained(static_cast<std::size_t>(k.cols()), 0);
            Eigen::VectorXd prescribed = Eigen::VectorXd::Zero(k.cols());
            for (int i = 0; i < static_cast<int>(bc.dofs.size()); ++i)
            {
                const int c = bc.dofs[static_cast<std::size_t>(i)];
                constrained[static_cast<std::size_t>(c)] = 1;
                prescribed(c) = scale * bc.values(i) - uu(c);
            }

            for (int row = 0; row < k.outerSize(); ++row)
            {
                for (SparseMatrix::InnerIterator it(k, row); it; ++it)
                {
                    if (constrained[static_cast<std::size_t>(it.col())])
                    {
                        rhs(row) -= it.value() * prescribed(it.col());
                        if (row != it.col())
                        {
                            it.valueRef() = 0.0;
                        }
                    }
                }
            }

            for (int i = 0; i < static_cast<int>(bc.dofs.size()); ++i)
            {
                const int c = bc.dofs[static_cast<std::size_t>(i)];
                for (SparseMatrix::InnerIterator it(k, c); it; ++it)
                {
                    if (it.col() != c)
                    {
                        it.valueRef() = 0.0;
                    }
                }
                k.coeffRef(c, c) = 1.0;
                rhs(c) = prescribed(c);
            }
            k.prune(0.0);
        }

        void apply_nonlinear_bc_row_only(const BoundaryCondition &bc,
                                         const Eigen::VectorXd &uu,
                                         double scale,
                                         SparseMatrix &k,
                                         Eigen::VectorXd &rhs)
        {
            for (int i = 0; i < static_cast<int>(bc.dofs.size()); ++i)
            {
                const int c = bc.dofs[static_cast<std::size_t>(i)];
                for (SparseMatrix::InnerIterator it(k, c); it; ++it)
                {
                    it.valueRef() = (it.col() == c) ? 1.0 : 0.0;
                }
                k.coeffRef(c, c) = 1.0;
                rhs(c) = scale * bc.values(i) - uu(c);
            }
            k.prune(0.0);
        }

        Eigen::MatrixXd element_coords(const Mesh &mesh, const std::array<int, 6> &element)
        {
            Eigen::MatrixXd coords(6, 2);
            for (int a = 0; a < 6; ++a)
            {
                coords.row(a) = mesh.nodes.row(element[static_cast<std::size_t>(a)]);
            }
            return coords;
        }

        Eigen::MatrixXd nodal_matrix(const Eigen::VectorXd &uu)
        {
            Eigen::MatrixXd u(uu.size() / 2, 2);
            for (int i = 0; i < u.rows(); ++i)
            {
                u(i, 0) = uu(2 * i);
                u(i, 1) = uu(2 * i + 1);
            }
            return u;
        }

        int local_edge_index_from_corner_nodes(const std::array<int, 6> &element, int node_a, int node_b)
        {
            for (int edge = 0; edge < 3; ++edge)
            {
                const int edge_a = element[static_cast<std::size_t>(edge)];
                const int edge_b = element[static_cast<std::size_t>((edge + 1) % 3)];
                if ((edge_a == node_a && edge_b == node_b) || (edge_a == node_b && edge_b == node_a))
                {
                    return edge;
                }
            }
            throw std::runtime_error("Failed to match a shared T6 edge while building ESFEM support data");
        }

        SupportData build_cell_support_data(const Model &model, int elem_index)
        {
            const auto &element = model.mesh.elements[static_cast<std::size_t>(elem_index)];
            const Eigen::MatrixXd wkx = element_coords(model.mesh, element);
            const QuadratureRule wi = triangle_rule(2);
            const QuadratureRule wb = gauss_line_rule(2);

            Eigen::MatrixXd xy = wkx;
            Eigen::RowVectorXd ni = Eigen::RowVectorXd::Zero(6);
            Eigen::MatrixXd wmat(3, wi.weights.size());
            std::vector<double> detj0(static_cast<std::size_t>(wi.weights.size()));
            Eigen::MatrixXd mr(wi.weights.size(), 2);

            for (int ig = 0; ig < wi.weights.size(); ++ig)
            {
                Eigen::VectorXd n1;
                Eigen::MatrixXd dndxi1;
                lagrange_basis("T3", wi.points.row(ig), n1, dndxi1);
                detj0[static_cast<std::size_t>(ig)] = (dndxi1.transpose() * xy.topRows(3)).determinant();
                mr.row(ig) = (xy.topRows(3).transpose() * n1).transpose();
                const Eigen::VectorXd n = serendipity_shape("T6", wkx, mr.row(ig));
                ni += (n * wi.weights(ig) * detj0[static_cast<std::size_t>(ig)]).transpose();
                wmat.col(ig) << wi.weights(ig) * detj0[static_cast<std::size_t>(ig)],
                    wi.weights(ig) * detj0[static_cast<std::size_t>(ig)] * mr(ig, 0),
                    wi.weights(ig) * detj0[static_cast<std::size_t>(ig)] * mr(ig, 1);
            }

            const int bound[3][3] = {{0, 1, 3}, {1, 2, 4}, {2, 0, 5}};
            std::vector<double> nx;
            std::vector<double> ny;
            outward_normals(xy.topRows(3), side_lengths(xy.topRows(3)), nx, ny);

            Eigen::MatrixXd fx = Eigen::MatrixXd::Zero(3, 6);
            Eigen::MatrixXd fy = Eigen::MatrixXd::Zero(3, 6);
            for (int is = 0; is < 3; ++is)
            {
                Eigen::MatrixXd bxy = Eigen::MatrixXd::Zero(3, 6);
                Eigen::MatrixXd x(3, 2);
                for (int a = 0; a < 3; ++a)
                {
                    x.row(a) = xy.row(bound[is][a]);
                }
                for (int ig = 0; ig < wb.weights.size(); ++ig)
                {
                    Eigen::VectorXd ng;
                    Eigen::MatrixXd dndxi;
                    lagrange_basis("L3", Eigen::Vector2d(wb.points(ig, 0), 0.0), ng, dndxi);
                    const Eigen::Vector2d gp = x.transpose() * ng;
                    const double detj = (x.transpose() * dndxi).norm();
                    const Eigen::VectorXd nt = serendipity_shape("T6", wkx, gp);
                    bxy.row(0) += (nt * detj * wb.weights(ig)).transpose();
                    bxy.row(1) += (nt * detj * wb.weights(ig) * gp(0)).transpose();
                    bxy.row(2) += (nt * detj * wb.weights(ig) * gp(1)).transpose();
                }
                fx += nx[static_cast<std::size_t>(is)] * bxy;
                fy += ny[static_cast<std::size_t>(is)] * bxy;
            }
            fx.row(1) -= ni;
            fy.row(2) -= ni;

            SupportData data;
            data.nodes.assign(element.begin(), element.end());
            data.gcoord = wkx;
            data.fx = fx;
            data.fy = fy;
            const auto wmat_lu = wmat.fullPivLu();
            data.dx = wmat_lu.solve(fx);
            data.dy = wmat_lu.solve(fy);
            data.points = mr;
            data.weights = detj0;
            data.sub_area = triangle_area_from_vertices(wkx.topRows(3));
            data.edof = element_dofs(element);
            return data;
        }

        std::vector<TargetEdge> build_target_edges(const Model &model, std::vector<std::vector<int>> &supp, std::vector<double> &sub_area)
        {
            for (const auto &element : model.mesh.elements)
            {
                supp.emplace_back(element.begin(), element.end());
            }

            std::vector<TargetEdge> edges;
            for (int ie = 0; ie < static_cast<int>(model.mesh.elements.size()); ++ie)
            {
                const auto &element = model.mesh.elements[static_cast<std::size_t>(ie)];
                for (int j = 0; j < 3; ++j)
                {
                    const int n1 = element[static_cast<std::size_t>(j)];
                    const int n2 = element[static_cast<std::size_t>((j + 1) % 3)];
                    bool found = false;
                    for (auto &edge : edges)
                    {
                        if ((edge.n1 == n1 && edge.n2 == n2) || (edge.n1 == n2 && edge.n2 == n1))
                        {
                            edge.elem2 = ie;
                            edge.local_edge2 = j;
                            found = true;
                            break;
                        }
                    }
                    if (!found)
                    {
                        edges.push_back(TargetEdge{n1, n2, ie, -1, j, -1});
                    }
                }
            }

            sub_area.resize(edges.size(), 0.0);
            std::vector<double> tri_area(model.mesh.elements.size(), 0.0);
            for (int ie = 0; ie < static_cast<int>(model.mesh.elements.size()); ++ie)
            {
                tri_area[static_cast<std::size_t>(ie)] = triangle_area_from_vertices(element_coords(model.mesh, model.mesh.elements[static_cast<std::size_t>(ie)]).topRows(3));
            }

            for (int i = 0; i < static_cast<int>(edges.size()); ++i)
            {
                sub_area[static_cast<std::size_t>(i)] = tri_area[static_cast<std::size_t>(edges[static_cast<std::size_t>(i)].elem1)] / 3.0;
                if (edges[static_cast<std::size_t>(i)].elem2 >= 0)
                {
                    sub_area[static_cast<std::size_t>(i)] += tri_area[static_cast<std::size_t>(edges[static_cast<std::size_t>(i)].elem2)] / 3.0;
                }
            }
            return edges;
        }

        SupportData build_edge_support_data(const Model &model,
                                            const std::vector<TargetEdge> &target_edges,
                                            const std::vector<std::vector<int>> &supp,
                                            const std::vector<double> &sub_area,
                                            int edge_index,
                                            bool nonlinear)
        {
            const TargetEdge &target = target_edges[static_cast<std::size_t>(edge_index)];
            std::vector<int> neighbour = {target.elem1};
            if (target.elem2 >= 0)
            {
                neighbour.push_back(target.elem2);
            }
            const int nc = static_cast<int>(neighbour.size());
            const QuadratureRule wb = gauss_line_rule(nonlinear ? 3 : 2);

            const int m = nonlinear ? 6 : (nc == 1 ? 3 : 4);
            Eigen::MatrixXd fx;
            Eigen::MatrixXd fy;
            std::vector<int> nodl;
            int nn = 0;

            const int bound[3][3] = {{0, 1, 3}, {1, 2, 4}, {2, 0, 5}};
            for (int ic = 0; ic < nc; ++ic)
            {
                const auto &element = model.mesh.elements[static_cast<std::size_t>(neighbour[static_cast<std::size_t>(ic)])];
                const Eigen::MatrixXd wkx = element_coords(model.mesh, element);
                std::vector<double> nx;
                std::vector<double> ny;
                outward_normals(wkx.topRows(3), side_lengths(wkx.topRows(3)), nx, ny);

                Eigen::MatrixXd lfx = Eigen::MatrixXd::Zero(m, 6);
                Eigen::MatrixXd lfy = Eigen::MatrixXd::Zero(m, 6);
                for (int is = 0; is < 3; ++is)
                {
                    Eigen::MatrixXd bxy = Eigen::MatrixXd::Zero(m, 6);
                    Eigen::MatrixXd x(3, 2);
                    for (int a = 0; a < 3; ++a)
                    {
                        x.row(a) = wkx.row(bound[is][a]);
                    }
                    for (int ig = 0; ig < wb.weights.size(); ++ig)
                    {
                        Eigen::VectorXd ng;
                        Eigen::MatrixXd dndxi;
                        lagrange_basis("L3", Eigen::Vector2d(wb.points(ig, 0), 0.0), ng, dndxi);
                        const Eigen::Vector2d gp = x.transpose() * ng;
                        const double detj = (x.transpose() * dndxi).norm();
                        const Eigen::VectorXd nt = serendipity_shape("T6", wkx, gp);
                        const Eigen::RowVectorXd n = (nt * detj * wb.weights(ig)).transpose();

                        if (nonlinear)
                        {
                            Eigen::VectorXd p(6);
                            p << 1.0, gp(0), gp(1), gp(0) * gp(0), gp(0) * gp(1), gp(1) * gp(1);
                            bxy += p * n;
                        }
                        else if (nc == 1)
                        {
                            bxy.row(0) += n;
                            bxy.row(1) += n * gp(0);
                            bxy.row(2) += n * gp(1);
                        }
                        else
                        {
                            bxy.row(0) += n;
                            bxy.row(1) += n * gp(0);
                            bxy.row(2) += n * gp(1);
                            bxy.row(3) += n * gp(0) * gp(1);
                        }
                    }
                    lfx += nx[static_cast<std::size_t>(is)] * bxy;
                    lfy += ny[static_cast<std::size_t>(is)] * bxy;
                }

                if (ic == 0)
                {
                    nodl = supp[static_cast<std::size_t>(neighbour[static_cast<std::size_t>(ic)])];
                    nn = static_cast<int>(nodl.size());
                    fx = lfx;
                    fy = lfy;
                }
                else
                {
                    for (int jj = 0; jj < static_cast<int>(supp[static_cast<std::size_t>(neighbour[static_cast<std::size_t>(ic)])].size()); ++jj)
                    {
                        const int node = supp[static_cast<std::size_t>(neighbour[static_cast<std::size_t>(ic)])][static_cast<std::size_t>(jj)];
                        auto it = std::find(nodl.begin(), nodl.end(), node);
                        if (it != nodl.end())
                        {
                            const int pos = static_cast<int>(std::distance(nodl.begin(), it));
                            fx.col(pos) += lfx.col(jj);
                            fy.col(pos) += lfy.col(jj);
                        }
                        else
                        {
                            nodl.push_back(node);
                            fx.conservativeResize(Eigen::NoChange, fx.cols() + 1);
                            fy.conservativeResize(Eigen::NoChange, fy.cols() + 1);
                            fx.col(fx.cols() - 1) = lfx.col(jj);
                            fy.col(fy.cols() - 1) = lfy.col(jj);
                        }
                    }
                    nn = static_cast<int>(nodl.size());
                }
            }

            std::string element_type;
            QuadratureRule wi;
            std::vector<int> reordered_nodes;
            if (nc == 1)
            {
                reordered_nodes.assign(nodl.begin(), nodl.end());
                element_type = "T6";
                wi = triangle_rule(nonlinear ? 4 : 2);
            }
            else
            {
                const auto &elem1 = model.mesh.elements[static_cast<std::size_t>(target.elem1)];
                const auto &elem2 = model.mesh.elements[static_cast<std::size_t>(target.elem2)];
                const int shared_local_a = target.local_edge1;
                const int shared_local_b = (target.local_edge1 + 1) % 3;
                const int opposite_local = (target.local_edge1 + 2) % 3;

                const int shared_a = elem1[static_cast<std::size_t>(shared_local_a)];
                const int shared_b = elem1[static_cast<std::size_t>(shared_local_b)];
                const int opposite1 = elem1[static_cast<std::size_t>(opposite_local)];
                const int mid_shared = elem1[static_cast<std::size_t>(3 + target.local_edge1)];
                const int mid_b_opp1 = elem1[static_cast<std::size_t>(3 + ((target.local_edge1 + 1) % 3))];
                const int mid_opp1_a = elem1[static_cast<std::size_t>(3 + ((target.local_edge1 + 2) % 3))];

                int opposite2 = -1;
                for (int local = 0; local < 3; ++local)
                {
                    const int node = elem2[static_cast<std::size_t>(local)];
                    if (node != shared_a && node != shared_b)
                    {
                        opposite2 = node;
                        break;
                    }
                }
                if (opposite2 < 0)
                {
                    throw std::runtime_error("Failed to locate the opposite node of a shared edge in Cook ESFEM support data");
                }

                const int edge_a_opp2 = local_edge_index_from_corner_nodes(elem2, shared_a, opposite2);
                const int edge_opp2_b = local_edge_index_from_corner_nodes(elem2, opposite2, shared_b);
                const int mid_a_opp2 = elem2[static_cast<std::size_t>(3 + edge_a_opp2)];
                const int mid_opp2_b = elem2[static_cast<std::size_t>(3 + edge_opp2_b)];

                if (target.local_edge1 == 0)
                {
                    reordered_nodes = {shared_a, opposite2, shared_b, opposite1,
                                       mid_a_opp2, mid_opp2_b, mid_b_opp1, mid_opp1_a, mid_shared};
                }
                else if (target.local_edge1 == 1)
                {
                    reordered_nodes = {opposite1, shared_a, opposite2, shared_b,
                                       mid_opp1_a, mid_a_opp2, mid_opp2_b, mid_b_opp1, mid_shared};
                }
                else
                {
                    reordered_nodes = {shared_b, opposite1, shared_a, opposite2,
                                       mid_b_opp1, mid_opp1_a, mid_a_opp2, mid_opp2_b, mid_shared};
                }

                element_type = "Q9";
                wi = gauss_quad_rule(nonlinear ? 3 : 2);
            }

            std::unordered_map<int, int> support_pos;
            support_pos.reserve(nodl.size());
            for (int i = 0; i < static_cast<int>(nodl.size()); ++i)
            {
                support_pos.emplace(nodl[static_cast<std::size_t>(i)], i);
            }

            Eigen::MatrixXd gcoord(reordered_nodes.size(), 2);
            Eigen::MatrixXd rfx(fx.rows(), reordered_nodes.size());
            Eigen::MatrixXd rfy(fy.rows(), reordered_nodes.size());
            for (int i = 0; i < static_cast<int>(reordered_nodes.size()); ++i)
            {
                gcoord.row(i) = model.mesh.nodes.row(reordered_nodes[static_cast<std::size_t>(i)]);
                const int src = support_pos.at(reordered_nodes[static_cast<std::size_t>(i)]);
                rfx.col(i) = fx.col(src);
                rfy.col(i) = fy.col(src);
            }

            Eigen::MatrixXd ni = Eigen::MatrixXd::Zero(nonlinear ? 5 : 3, gcoord.rows());
            Eigen::MatrixXd wmat((nonlinear ? 6 : (nc == 1 ? 3 : 4)), wi.weights.size());
            Eigen::MatrixXd pmat(wi.weights.size(), nonlinear ? 6 : (nc == 1 ? 3 : 4));
            std::vector<double> mw(static_cast<std::size_t>(wi.weights.size()));
            Eigen::MatrixXd mq(wi.weights.size(), 2);
            Eigen::MatrixXd mmat = Eigen::MatrixXd::Zero(wmat.rows(), wmat.rows());

            for (int ig = 0; ig < wi.weights.size(); ++ig)
            {
                Eigen::VectorXd n1;
                Eigen::MatrixXd dndxi1;
                if (nc == 1)
                {
                    lagrange_basis(nonlinear ? "T6" : "T3", wi.points.row(ig), n1, dndxi1);
                    const int nmap = nonlinear ? 6 : 3;
                    const Eigen::MatrixXd base = gcoord.topRows(nmap);
                    const double detj = (dndxi1.transpose() * base).determinant();
                    mq.row(ig) = (base.transpose() * n1).transpose();
                    mw[static_cast<std::size_t>(ig)] = wi.weights(ig) * detj;
                }
                else
                {
                    lagrange_basis(nonlinear ? "Q9" : "Q4", wi.points.row(ig), n1, dndxi1);
                    const int nmap = nonlinear ? 9 : 4;
                    const Eigen::MatrixXd base = gcoord.topRows(nmap);
                    const double detj = (dndxi1.transpose() * base).determinant();
                    mq.row(ig) = (base.transpose() * n1).transpose();
                    mw[static_cast<std::size_t>(ig)] = wi.weights(ig) * detj;
                }

                const Eigen::VectorXd n = serendipity_shape(element_type, gcoord, mq.row(ig));
                if (nonlinear)
                {
                    Eigen::VectorXd p(6);
                    p << 1.0, mq(ig, 0), mq(ig, 1), mq(ig, 0) * mq(ig, 0), mq(ig, 0) * mq(ig, 1), mq(ig, 1) * mq(ig, 1);
                    pmat.row(ig) = p.transpose();
                    mmat += mw[static_cast<std::size_t>(ig)] * (p * p.transpose());
                    ni.row(0) += (n * mw[static_cast<std::size_t>(ig)]).transpose();
                    ni.row(1) += (2.0 * n * mw[static_cast<std::size_t>(ig)] * mq(ig, 0)).transpose();
                    ni.row(2) += (2.0 * n * mw[static_cast<std::size_t>(ig)] * mq(ig, 1)).transpose();
                    ni.row(3) += (n * mw[static_cast<std::size_t>(ig)] * mq(ig, 0)).transpose();
                    ni.row(4) += (n * mw[static_cast<std::size_t>(ig)] * mq(ig, 1)).transpose();
                }
                else if (nc == 1)
                {
                    ni.row(0) += (n * mw[static_cast<std::size_t>(ig)]).transpose();
                    wmat.col(ig) << mw[static_cast<std::size_t>(ig)], mw[static_cast<std::size_t>(ig)] * mq(ig, 0), mw[static_cast<std::size_t>(ig)] * mq(ig, 1);
                }
                else
                {
                    ni.row(0) += (n * mw[static_cast<std::size_t>(ig)]).transpose();
                    ni.row(1) += (n * mw[static_cast<std::size_t>(ig)] * mq(ig, 0)).transpose();
                    ni.row(2) += (n * mw[static_cast<std::size_t>(ig)] * mq(ig, 1)).transpose();
                    wmat.col(ig) << mw[static_cast<std::size_t>(ig)],
                        mw[static_cast<std::size_t>(ig)] * mq(ig, 0),
                        mw[static_cast<std::size_t>(ig)] * mq(ig, 1),
                        mw[static_cast<std::size_t>(ig)] * mq(ig, 0) * mq(ig, 1);
                }
            }

            if (nonlinear)
            {
                rfx.row(1) -= ni.row(0);
                rfx.row(3) -= ni.row(1);
                rfx.row(4) -= ni.row(4);
                rfy.row(2) -= ni.row(0);
                rfy.row(4) -= ni.row(3);
                rfy.row(5) -= ni.row(2);
            }
            else if (nc == 1)
            {
                rfx.row(1) -= ni.row(0);
                rfy.row(2) -= ni.row(0);
            }
            else
            {
                rfx.row(1) -= ni.row(0);
                rfx.row(3) -= ni.row(2);
                rfy.row(2) -= ni.row(0);
                rfy.row(3) -= ni.row(1);
            }

            Eigen::MatrixXd dx;
            Eigen::MatrixXd dy;
            if (nonlinear)
            {
                const double alpha = 1.0e-10 * std::max(mmat.trace() / 6.0, 1.0);
                const Eigen::MatrixXd mreg = mmat + alpha * Eigen::MatrixXd::Identity(6, 6);
                const auto mreg_lu = mreg.fullPivLu();
                dx = (pmat * mreg_lu.solve(rfx)).eval();
                dy = (pmat * mreg_lu.solve(rfy)).eval();
            }
            else
            {
                const auto wmat_lu = wmat.fullPivLu();
                dx = wmat_lu.solve(rfx);
                dy = wmat_lu.solve(rfy);
            }

            SupportData data;
            data.nodes = reordered_nodes;
            data.gcoord = gcoord;
            data.fx = rfx;
            data.fy = rfy;
            data.points = mq;
            data.weights = mw;
            data.dx = dx;
            data.dy = dy;
            data.sub_area = sub_area[static_cast<std::size_t>(edge_index)];
            data.edof.resize(2 * reordered_nodes.size());
            for (int i = 0; i < static_cast<int>(reordered_nodes.size()); ++i)
            {
                data.edof[static_cast<std::size_t>(2 * i)] = 2 * reordered_nodes[static_cast<std::size_t>(i)];
                data.edof[static_cast<std::size_t>(2 * i + 1)] = 2 * reordered_nodes[static_cast<std::size_t>(i)] + 1;
            }
            return data;
        }

        std::vector<SupportData> build_all_cell_support_data(const Model &model)
        {
            std::vector<SupportData> supports;
            supports.reserve(model.mesh.elements.size());
            for (int ie = 0; ie < static_cast<int>(model.mesh.elements.size()); ++ie)
            {
                supports.push_back(build_cell_support_data(model, ie));
            }
            return supports;
        }

        std::vector<SupportData> build_all_edge_support_data(const Model &model, bool nonlinear)
        {
            std::vector<std::vector<int>> supp;
            std::vector<double> sub_area;
            const std::vector<TargetEdge> target_edges = build_target_edges(model, supp, sub_area);

            std::vector<SupportData> supports;
            supports.reserve(target_edges.size());
            for (int ie = 0; ie < static_cast<int>(target_edges.size()); ++ie)
            {
                supports.push_back(build_edge_support_data(model, target_edges, supp, sub_area, ie, nonlinear));
            }
            return supports;
        }

        Result finalize_result(const Model &model,
                               const Eigen::VectorXd &uu,
                               double strain_energy,
                               const std::vector<double> &history,
                               bool compute_exact_solution)
        {
            Result result;
            result.uu = uu;
            result.nodal_u = nodal_matrix(uu);
            result.strain_energy_history = history;
            result.strain_energy = strain_energy;

            if (model.has_exact_solution && compute_exact_solution)
            {
                result.exact_u = exact_field(model.scenario, model.mesh.nodes);
                const Eigen::MatrixXd diff = result.exact_u - result.nodal_u;
                const double err = diff.array().square().sum();
                const double den = result.exact_u.array().square().sum();
                result.relative_error = den > 0.0 ? std::sqrt(err / den) : 0.0;
            }
            else
            {
                result.relative_error = std::numeric_limits<double>::quiet_NaN();
            }
            return result;
        }

    } // namespace

    std::string method_name(Method method)
    {
        switch (method)
        {
        case Method::FEM:
            return "fem";
        case Method::CSFEM:
            return "csfem";
        case Method::ESFEM:
            return "esfem";
        }
        throw std::invalid_argument("Unknown method");
    }

    Method parse_method(const std::string &value)
    {
        const std::string lower = [&]()
        {
            std::string s = value;
            std::transform(s.begin(), s.end(), s.begin(), [](unsigned char c)
                           { return static_cast<char>(std::tolower(c)); });
            return s;
        }();
        if (lower == "fem")
        {
            return Method::FEM;
        }
        if (lower == "csfem")
        {
            return Method::CSFEM;
        }
        if (lower == "esfem")
        {
            return Method::ESFEM;
        }
        throw std::invalid_argument("Unsupported method: " + value);
    }

    Result solve_linear_patch(const Model &model, const NewtonOptions &options)
    {
        const int ndof = 2 * model.mesh.nodes.rows();
        std::vector<Triplet> triplets;
        std::vector<SupportData> cell_supports;
        std::vector<SupportData> edge_supports;

        if (model.method == Method::CSFEM)
        {
            cell_supports = build_all_cell_support_data(model);
        }
        else if (model.method == Method::ESFEM)
        {
            edge_supports = build_all_edge_support_data(model, false);
        }

        if (model.method == Method::FEM)
        {
            const QuadratureRule qr = triangle_rule(2);
#pragma omp parallel for if (model.mesh.elements.size() > 1) default(none) shared(model, qr, triplets)
            for (int ie = 0; ie < static_cast<int>(model.mesh.elements.size()); ++ie)
            {
                const auto &element = model.mesh.elements[static_cast<std::size_t>(ie)];
                const Eigen::MatrixXd wkx = element_coords(model.mesh, element);
                const std::vector<int> edof = element_dofs(element);
                Eigen::MatrixXd ke = Eigen::MatrixXd::Zero(12, 12);
                for (int ig = 0; ig < qr.weights.size(); ++ig)
                {
                    Eigen::VectorXd n;
                    Eigen::MatrixXd dndxi;
                    lagrange_basis("T6", qr.points.row(ig), n, dndxi);
                    const Eigen::Matrix2d j0 = wkx.transpose() * dndxi;
                    const Eigen::MatrixXd dndx = dndxi * j0.inverse();
                    Eigen::MatrixXd b = Eigen::MatrixXd::Zero(3, 12);
                    for (int a = 0; a < 6; ++a)
                    {
                        b(0, 2 * a) = dndx(a, 0);
                        b(1, 2 * a + 1) = dndx(a, 1);
                        b(2, 2 * a) = dndx(a, 1);
                        b(2, 2 * a + 1) = dndx(a, 0);
                    }
                    ke += b.transpose() * model.cmat * b * qr.weights(ig) * j0.determinant();
                }
#pragma omp critical(structural_linear_fem_scatter)
                scatter_add(triplets, edof, ke);
            }
        }
        else if (model.method == Method::CSFEM)
        {
#pragma omp parallel for if (cell_supports.size() > 1) default(none) shared(model, triplets, cell_supports)
            for (int ie = 0; ie < static_cast<int>(cell_supports.size()); ++ie)
            {
                const SupportData &data = cell_supports[static_cast<std::size_t>(ie)];
                Eigen::MatrixXd ke = Eigen::MatrixXd::Zero(12, 12);
                for (int ig = 0; ig < data.points.rows(); ++ig)
                {
                    Eigen::MatrixXd b = Eigen::MatrixXd::Zero(3, 12);
                    for (int a = 0; a < 6; ++a)
                    {
                        b(0, 2 * a) = data.dx(ig, a);
                        b(1, 2 * a + 1) = data.dy(ig, a);
                        b(2, 2 * a) = data.dy(ig, a);
                        b(2, 2 * a + 1) = data.dx(ig, a);
                    }
                    ke += b.transpose() * model.cmat * b * triangle_rule(2).weights(ig) * data.weights[static_cast<std::size_t>(ig)];
                }
#pragma omp critical(structural_linear_csfem_scatter)
                scatter_add(triplets, data.edof, ke);
            }
        }
        else
        {
#pragma omp parallel for if (edge_supports.size() > 1) default(none) shared(model, edge_supports, triplets)
            for (int ie = 0; ie < static_cast<int>(edge_supports.size()); ++ie)
            {
                const SupportData &data = edge_supports[static_cast<std::size_t>(ie)];
                Eigen::MatrixXd ke = Eigen::MatrixXd::Zero(data.edof.size(), data.edof.size());
                for (int ig = 0; ig < data.points.rows(); ++ig)
                {
                    Eigen::MatrixXd b = Eigen::MatrixXd::Zero(3, data.edof.size());
                    for (int a = 0; a < static_cast<int>(data.nodes.size()); ++a)
                    {
                        b(0, 2 * a) = data.dx(ig, a);
                        b(1, 2 * a + 1) = data.dy(ig, a);
                        b(2, 2 * a) = data.dy(ig, a);
                        b(2, 2 * a + 1) = data.dx(ig, a);
                    }
                    ke += b.transpose() * model.cmat * b * data.weights[static_cast<std::size_t>(ig)];
                }
#pragma omp critical(structural_linear_esfem_scatter)
                scatter_add(triplets, data.edof, ke);
            }
        }

        SparseMatrix k = make_sparse_matrix(ndof, ndof, triplets);
        Eigen::VectorXd f = model.force;
        apply_linear_bc(model.bc, k, f);
        const Eigen::VectorXd uu = solve_system(k, f);
        return finalize_result(model, uu, 0.5 * uu.dot(k * uu), {}, options.compute_exact_solution);
    }

    Result solve_nonlinear_patch(const Model &model, const NewtonOptions &options)
    {
        const bool use_shared_solver =
            model.method == Method::ESFEM ||
            (model.method == Method::CSFEM && model.scenario != Scenario::BendingBlock);
        if (use_shared_solver)
        {
            const nonlinear_esfem_t6::Problem shared_problem = to_shared_problem(model);
            const nonlinear_esfem_t6::SolverOptions shared_options = to_shared_options(options);
            const shared_sfem::NonlinearSmoothedFemSolver solver;

            if (model.method == Method::CSFEM)
            {
                const shared_sfem::CellSmoothedDomainAssembler assembler;
                const nonlinear_esfem_t6::Result shared_result = solver.solve(shared_problem, shared_options, assembler);
                return finalize_result(model,
                                       shared_result.uu,
                                       shared_result.strain_energy,
                                       shared_result.strain_energy_history,
                                       options.compute_exact_solution);
            }

            const shared_sfem::EdgeSmoothedDomainAssembler assembler;
            const nonlinear_esfem_t6::Result shared_result = solver.solve(shared_problem, shared_options, assembler);
            return finalize_result(model,
                                   shared_result.uu,
                                   shared_result.strain_energy,
                                   shared_result.strain_energy_history,
                                   options.compute_exact_solution);
        }

        const int ndof = 2 * model.mesh.nodes.rows();
        Eigen::VectorXd uu = Eigen::VectorXd::Zero(ndof);
        SparseDirectSolverCache solver_cache;
        std::vector<double> history;
        history.reserve(options.nstep);

        std::vector<SupportData> cell_supports;
        std::vector<SupportData> edge_supports;
        if (model.method == Method::CSFEM)
        {
            cell_supports = build_all_cell_support_data(model);
        }
        else if (model.method == Method::ESFEM)
        {
            edge_supports = build_all_edge_support_data(model, true);
        }

        double last_energy = 0.0;
        for (int istp = 1; istp <= options.nstep; ++istp)
        {
            const double scale = static_cast<double>(istp) / static_cast<double>(options.nstep);
            double condition = 1.0;
            int niter = 0;
            bool unstable = false;

            std::cout << "\n Step " << istp << "\t Scale " << (scale * 100.0) << "%\n";

            while (condition > options.tolerance && niter < options.maxiter)
            {
                ++niter;
                const auto iter_begin = std::chrono::steady_clock::now();
                double energy = 0.0;
                int nthreads = 1;
#ifdef _OPENMP
                nthreads = omp_get_max_threads();
#endif
                std::vector<std::vector<Triplet>> triplet_buffers(static_cast<std::size_t>(nthreads));
                std::vector<Eigen::VectorXd> residual_buffers(static_cast<std::size_t>(nthreads), Eigen::VectorXd::Zero(ndof));

                if (model.method == Method::FEM)
                {
                    const QuadratureRule qr = triangle_rule(3);
                    const std::size_t triplet_estimate = model.mesh.elements.size() * 144;
                    for (auto &buffer : triplet_buffers)
                    {
                        buffer.reserve(triplet_estimate / static_cast<std::size_t>(nthreads) + 256);
                    }
#pragma omp parallel for if (model.mesh.elements.size() > 1) default(none) shared(model, qr, uu, triplet_buffers, residual_buffers, nthreads) reduction(+ : energy)
                    for (int ie = 0; ie < static_cast<int>(model.mesh.elements.size()); ++ie)
                    {
#ifdef _OPENMP
                        const int tid = omp_get_thread_num();
#else
                        const int tid = 0;
#endif
                        const auto &element = model.mesh.elements[static_cast<std::size_t>(ie)];
                        const Eigen::MatrixXd wkx = element_coords(model.mesh, element);
                        const std::vector<int> edof = element_dofs(element);
                        Eigen::MatrixXd wku(6, 2);
                        Eigen::MatrixXd bmat;
                        Eigen::MatrixXd bgeo;
                        for (int a = 0; a < 6; ++a)
                        {
                            wku.row(a) << uu(2 * element[static_cast<std::size_t>(a)]), uu(2 * element[static_cast<std::size_t>(a)] + 1);
                        }
                        Eigen::MatrixXd ke = Eigen::MatrixXd::Zero(12, 12);
                        Eigen::VectorXd re = Eigen::VectorXd::Zero(12);
                        for (int ig = 0; ig < qr.weights.size(); ++ig)
                        {
                            Eigen::VectorXd n;
                            Eigen::MatrixXd dndxi;
                            lagrange_basis("T6", qr.points.row(ig), n, dndxi);
                            const Eigen::Matrix2d j0 = wkx.transpose() * dndxi;
                            const Eigen::MatrixXd dndx = dndxi * j0.inverse();
                            Eigen::Matrix3d fmat;
                            nonlinear_bmat(dndx, wku, bmat, bgeo, fmat);
                            const NonlinearMaterial material = nonlinear_constitutive(model.nonlinear_material, fmat);
                            ke += (bmat.transpose() * material.cmat * bmat + bgeo.transpose() * material.smat * bgeo) * qr.weights(ig) * j0.determinant();
                            Eigen::Vector3d svec(material.smat(0, 0), material.smat(1, 1), material.smat(0, 1));
                            re += bmat.transpose() * svec * qr.weights(ig) * j0.determinant();
                            energy += material.w0 * triangle_area_from_vertices(wkx.topRows(3));
                        }
                        scatter_add(triplet_buffers[static_cast<std::size_t>(tid)], edof, ke);
                        scatter_add(residual_buffers[static_cast<std::size_t>(tid)], edof, re);
                    }
                    energy /= static_cast<double>(qr.weights.size());
                }
                else if (model.method == Method::CSFEM)
                {
                    const QuadratureRule qr = triangle_rule(3);
                    const std::size_t triplet_estimate = cell_supports.size() * 144;
                    for (auto &buffer : triplet_buffers)
                    {
                        buffer.reserve(triplet_estimate / static_cast<std::size_t>(nthreads) + 256);
                    }
#pragma omp parallel for if (cell_supports.size() > 1) default(none) shared(model, options, uu, triplet_buffers, residual_buffers, nthreads, qr, cell_supports, istp, niter, scale) reduction(+ : energy)
                    for (int ie = 0; ie < static_cast<int>(cell_supports.size()); ++ie)
                    {
#ifdef _OPENMP
                        const int tid = omp_get_thread_num();
#else
                        const int tid = 0;
#endif
                        const SupportData &data = cell_supports[static_cast<std::size_t>(ie)];
                        Eigen::MatrixXd ke = Eigen::MatrixXd::Zero(12, 12);
                        Eigen::VectorXd re = Eigen::VectorXd::Zero(12);
                        Eigen::MatrixXd wku(6, 2);
                        Eigen::MatrixXd dndx(6, 2);
                        Eigen::MatrixXd bmat;
                        Eigen::MatrixXd bgeo;
                        for (int a = 0; a < 6; ++a)
                        {
                            wku.row(a) << uu(2 * data.nodes[static_cast<std::size_t>(a)]), uu(2 * data.nodes[static_cast<std::size_t>(a)] + 1);
                        }
                        for (int ig = 0; ig < data.points.rows(); ++ig)
                        {
                            dndx.col(0) = data.dx.row(ig).transpose();
                            dndx.col(1) = data.dy.row(ig).transpose();
                            Eigen::Matrix3d fmat;
                            nonlinear_bmat(dndx, wku, bmat, bgeo, fmat);
                            const NonlinearMaterial material = nonlinear_constitutive(model.nonlinear_material, fmat);
                            const double w = qr.weights(ig) * data.weights[static_cast<std::size_t>(ig)];
                            const bool bad_material =
                                !material.cmat.allFinite() || !material.smat.allFinite() || !std::isfinite(material.w0) ||
                                !fmat.allFinite() || !bmat.allFinite() || !bgeo.allFinite() || !dndx.allFinite() || !std::isfinite(w);
                            if (bad_material)
                            {
                                write_csfem_bending_failure_dump(model, options, data, wku, dndx, bmat, bgeo, fmat,
                                                                 material, istp, niter, ie, ig, scale, w, "nonfinite_material_or_kinematics");
                                continue;
                            }
                            ke += (bmat.transpose() * material.cmat * bmat + bgeo.transpose() * material.smat * bgeo) * w;
                            Eigen::Vector3d svec(material.smat(0, 0), material.smat(1, 1), material.smat(0, 1));
                            re += bmat.transpose() * svec * w;
                            if (!ke.allFinite() || !re.allFinite())
                            {
                                write_csfem_bending_failure_dump(model, options, data, wku, dndx, bmat, bgeo, fmat,
                                                                 material, istp, niter, ie, ig, scale, w, "nonfinite_element_contribution");
                                continue;
                            }
                            energy += material.w0 * data.sub_area;
                        }
                        scatter_add(triplet_buffers[static_cast<std::size_t>(tid)], data.edof, ke);
                        scatter_add(residual_buffers[static_cast<std::size_t>(tid)], data.edof, re);
                    }
                    energy /= 3.0;
                }
                else
                {
                    const std::size_t triplet_estimate = edge_supports.size() * 400;
                    for (auto &buffer : triplet_buffers)
                    {
                        buffer.reserve(triplet_estimate / static_cast<std::size_t>(nthreads) + 256);
                    }
#pragma omp parallel for if (edge_supports.size() > 1) default(none) shared(model, uu, triplet_buffers, residual_buffers, nthreads, edge_supports) reduction(+ : energy)
                    for (int ie = 0; ie < static_cast<int>(edge_supports.size()); ++ie)
                    {
#ifdef _OPENMP
                        const int tid = omp_get_thread_num();
#else
                        const int tid = 0;
#endif
                        const SupportData &data = edge_supports[static_cast<std::size_t>(ie)];
                        Eigen::MatrixXd ke = Eigen::MatrixXd::Zero(data.edof.size(), data.edof.size());
                        Eigen::VectorXd re = Eigen::VectorXd::Zero(data.edof.size());
                        Eigen::MatrixXd wku(data.nodes.size(), 2);
                        Eigen::MatrixXd dndx(data.nodes.size(), 2);
                        Eigen::MatrixXd bmat;
                        Eigen::MatrixXd bgeo;
                        for (int a = 0; a < static_cast<int>(data.nodes.size()); ++a)
                        {
                            wku.row(a) << uu(2 * data.nodes[static_cast<std::size_t>(a)]), uu(2 * data.nodes[static_cast<std::size_t>(a)] + 1);
                        }
                        double w0_acc = 0.0;
                        double w_acc = 0.0;
                        for (int ig = 0; ig < data.points.rows(); ++ig)
                        {
                            dndx.col(0) = data.dx.row(ig).transpose();
                            dndx.col(1) = data.dy.row(ig).transpose();
                            Eigen::Matrix3d fmat;
                            nonlinear_bmat(dndx, wku, bmat, bgeo, fmat);
                            const NonlinearMaterial material = nonlinear_constitutive(model.nonlinear_material, fmat);
                            const double w = data.weights[static_cast<std::size_t>(ig)];
                            ke += (bmat.transpose() * material.cmat * bmat + bgeo.transpose() * material.smat * bgeo) * w;
                            Eigen::Vector3d svec(material.smat(0, 0), material.smat(1, 1), material.smat(0, 1));
                            re += bmat.transpose() * svec * w;
                            w0_acc += material.w0 * w;
                            w_acc += w;
                        }
                        if (w_acc > 0.0)
                        {
                            energy += (w0_acc / w_acc) * data.sub_area;
                        }
                        scatter_add(triplet_buffers[static_cast<std::size_t>(tid)], data.edof, ke);
                        scatter_add(residual_buffers[static_cast<std::size_t>(tid)], data.edof, re);
                    }
                }

                std::vector<Triplet> triplets;
                std::size_t merged_triplets = 0;
                for (const auto &buffer : triplet_buffers)
                {
                    merged_triplets += buffer.size();
                }
                triplets.reserve(merged_triplets);
                Eigen::VectorXd r = Eigen::VectorXd::Zero(ndof);
                for (int tid = 0; tid < nthreads; ++tid)
                {
                    auto &buffer = triplet_buffers[static_cast<std::size_t>(tid)];
                    triplets.insert(triplets.end(), buffer.begin(), buffer.end());
                    r += residual_buffers[static_cast<std::size_t>(tid)];
                }
                const auto assembly_end = std::chrono::steady_clock::now();

                SparseMatrix k = make_sparse_matrix(ndof, ndof, triplets);
                const auto matrix_end = std::chrono::steady_clock::now();
                Eigen::VectorXd rhs = scale * model.force - r;
                if (model.method == Method::CSFEM && model.scenario == Scenario::BendingBlock)
                {
                    apply_nonlinear_bc_row_only(model.bc, uu, scale, k, rhs);
                }
                else
                {
                    apply_nonlinear_bc(model.bc, uu, scale, k, rhs);
                }
                const int pinned_rows = regularize_near_zero_rows(k, rhs);
                if (pinned_rows > 0 && model.scenario == Scenario::Cook && istp == 1 && niter == 1)
                {
                    std::cout << "Pinned " << pinned_rows << " near-zero tangent rows for cook at the first Newton step.\n";
                }
                if (model.method == Method::CSFEM && model.scenario == Scenario::BendingBlock && !cell_supports.empty())
                {
                    write_csfem_bending_debug_dump(model, options, cell_supports.front(), uu, rhs, k, istp, niter, scale);
                }
                const auto bc_end = std::chrono::steady_clock::now();
                const double residual = rhs.norm() / static_cast<double>(ndof);
                const Eigen::VectorXd duu = solver_cache.solve(k, rhs, options);
                const auto solve_end = std::chrono::steady_clock::now();
                if (!duu.allFinite())
                {
                    unstable = true;
                    break;
                }
                uu += duu;
                if (!uu.allFinite() || !std::isfinite(energy))
                {
                    uu -= duu;
                    unstable = true;
                    break;
                }
                condition = duu.norm() / std::max(uu.norm(), std::numeric_limits<double>::epsilon());
                if (!std::isfinite(condition))
                {
                    uu -= duu;
                    unstable = true;
                    break;
                }
                last_energy = energy;

                std::cout << std::fixed << std::setprecision(6)
                          << "Iteration Number " << niter
                          << " Condition " << condition
                          << " Residual " << residual
                          << " Tolerance " << options.tolerance << '\n';
                std::cout << "Timing [ms] assembly="
                          << std::chrono::duration<double, std::milli>(assembly_end - iter_begin).count()
                          << " sparse_build="
                          << std::chrono::duration<double, std::milli>(matrix_end - assembly_end).count()
                          << " bc="
                          << std::chrono::duration<double, std::milli>(bc_end - matrix_end).count()
                          << " solve="
                          << std::chrono::duration<double, std::milli>(solve_end - bc_end).count()
                          << " total="
                          << std::chrono::duration<double, std::milli>(solve_end - iter_begin).count()
                          << '\n';
            }

            history.push_back(last_energy);
            if (unstable)
            {
                break;
            }
        }

        return finalize_result(model, uu, last_energy, history, options.compute_exact_solution);
    }

    StructuralProblem::StructuralProblem(Model model)
        : model_(std::move(model))
    {
    }

    const Model &StructuralProblem::data() const
    {
        return model_;
    }

    Model &StructuralProblem::data()
    {
        return model_;
    }

    Method StructuralProblem::method() const
    {
        return model_.method;
    }

    Scenario StructuralProblem::scenario() const
    {
        return model_.scenario;
    }

    const std::string &StructuralProblem::name() const
    {
        return model_.name;
    }

    StructuralProblem StructuralProblem::make_patch_test(Method method,
                                                         Scenario scenario,
                                                         const Eigen::Vector2i &num_els)
    {
        return StructuralProblem(build_model(method, scenario, num_els));
    }

    Method FemSolver::method() const
    {
        return Method::FEM;
    }

    Result FemSolver::solve(const StructuralProblem &problem, const NewtonOptions &options) const
    {
        if (!is_nonlinear_scenario(problem.scenario()))
        {
            return solve_linear_patch(problem.data(), options);
        }
        return solve_nonlinear_patch(problem.data(), options);
    }

    Method CellSmoothedFemSolver::method() const
    {
        return Method::CSFEM;
    }

    Result CellSmoothedFemSolver::solve(const StructuralProblem &problem, const NewtonOptions &options) const
    {
        if (!is_nonlinear_scenario(problem.scenario()))
        {
            return solve_linear_patch(problem.data(), options);
        }
        return solve_nonlinear_patch(problem.data(), options);
    }

    Method EdgeSmoothedFemSolver::method() const
    {
        return Method::ESFEM;
    }

    Result EdgeSmoothedFemSolver::solve(const StructuralProblem &problem, const NewtonOptions &options) const
    {
        if (!is_nonlinear_scenario(problem.scenario()))
        {
            return solve_linear_patch(problem.data(), options);
        }
        return solve_nonlinear_patch(problem.data(), options);
    }

    std::unique_ptr<StructuralSolver> make_solver(Method method)
    {
        switch (method)
        {
        case Method::FEM:
            return std::make_unique<FemSolver>();
        case Method::CSFEM:
            return std::make_unique<CellSmoothedFemSolver>();
        case Method::ESFEM:
            return std::make_unique<EdgeSmoothedFemSolver>();
        }
        throw std::invalid_argument("Unknown method");
    }

} // namespace fem
