#include "sfem_solver_utils.hpp"
#include "sfem_mesh.hpp"

#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseLU>
#ifdef USE_EIGEN_UMFPACK
#include <unsupported/Eigen/UmfPackSupport>
#endif

#include <cmath>
#include <limits>
#include <stdexcept>
#include <vector>

namespace nonlinear_esfem_t6::shared
{

namespace
{

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

Eigen::VectorXd solve_with_dense_fallback(const SparseMatrix &a, const Eigen::VectorXd &b)
{
    constexpr int kDenseFallbackMaxDofs = 6000;
    if (a.rows() > kDenseFallbackMaxDofs || a.cols() > kDenseFallbackMaxDofs)
    {
        throw std::runtime_error("Failed to factorize sparse tangent matrix");
    }

    const Eigen::MatrixXd dense = Eigen::MatrixXd(a);
    Eigen::VectorXd best_x;
    double best_residual = std::numeric_limits<double>::infinity();
    auto consider_candidate = [&](const Eigen::VectorXd &x) {
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
        const Eigen::FullPivLU<Eigen::MatrixXd> lu(dense);
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
        const Eigen::FullPivLU<Eigen::MatrixXd> lu(shifted);
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

Eigen::VectorXd solve_with_iterative_fallback(const SparseMatrix &a,
                                              const Eigen::VectorXd &b,
                                              int maxiter,
                                              double tolerance)
{
    auto try_bicgstab_ilut = [&](SparseMatrix mat,
                                 double diag_shift,
                                 int local_maxiter,
                                 double local_tolerance,
                                 double droptol,
                                 int fillfactor) -> Eigen::VectorXd {
        if (diag_shift > 0.0)
        {
            mat.diagonal().array() += diag_shift;
            mat.makeCompressed();
        }
        Eigen::BiCGSTAB<SparseMatrix, Eigen::IncompleteLUT<double>> iterative_solver;
        iterative_solver.setMaxIterations(local_maxiter);
        iterative_solver.setTolerance(local_tolerance);
        iterative_solver.preconditioner().setDroptol(droptol);
        iterative_solver.preconditioner().setFillfactor(fillfactor);
        iterative_solver.compute(mat);
        if (iterative_solver.info() != Eigen::Success)
        {
            return Eigen::VectorXd();
        }
        const Eigen::VectorXd x = iterative_solver.solve(b);
        if (iterative_solver.info() != Eigen::Success || !is_acceptable_solution(mat, b, x))
        {
            return Eigen::VectorXd();
        }
        return x;
    };

    auto try_bicgstab_diag = [&](SparseMatrix mat,
                                 double diag_shift,
                                 int local_maxiter,
                                 double local_tolerance) -> Eigen::VectorXd {
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
        if (iterative_solver.info() != Eigen::Success || !is_acceptable_solution(mat, b, x))
        {
            return Eigen::VectorXd();
        }
        return x;
    };

    Eigen::VectorXd x = try_bicgstab_ilut(a, 0.0, maxiter, tolerance, 1.0e-3, 20);
    if (x.size() != 0)
    {
        return x;
    }

    x = try_bicgstab_ilut(a, 1.0e-10, std::max(maxiter, 4000), std::min(tolerance, 1.0e-10), 1.0e-3, 20);
    if (x.size() != 0)
    {
        return x;
    }

    x = try_bicgstab_ilut(a, 1.0e-8, std::max(maxiter, 6000), std::min(tolerance, 1.0e-10), 1.0e-4, 30);
    if (x.size() != 0)
    {
        return x;
    }

    x = try_bicgstab_ilut(a, 1.0e-6, std::max(maxiter, 8000), std::min(tolerance, 1.0e-11), 1.0e-4, 40);
    if (x.size() != 0)
    {
        return x;
    }

    x = try_bicgstab_diag(a, 1.0e-8, std::max(maxiter, 6000), std::min(tolerance, 1.0e-10));
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
            if (is_acceptable_solution(trial, b, x))
            {
                return x;
            }
        }
    }
    if (a.rows() <= 2500)
    {
        return solve_with_dense_fallback(a, b);
    }
    throw std::runtime_error("Failed to factorize sparse tangent matrix");
}

} // namespace

struct LinearSystemSolver::Impl
{
    Eigen::SparseLU<SparseMatrix> sparse_lu;
    int rows = -1;
    int cols = -1;
    std::vector<int> outer_index;
    std::vector<int> inner_index;
    bool analyzed = false;

    bool same_pattern(const SparseMatrix &a) const
    {
        if (!analyzed || rows != a.rows() || cols != a.cols())
        {
            return false;
        }
        const int outer_size = a.outerSize() + 1;
        if (static_cast<int>(outer_index.size()) != outer_size || static_cast<int>(inner_index.size()) != a.nonZeros())
        {
            return false;
        }
        for (int i = 0; i < outer_size; ++i)
        {
            if (outer_index[static_cast<std::size_t>(i)] != a.outerIndexPtr()[i])
            {
                return false;
            }
        }
        for (int i = 0; i < a.nonZeros(); ++i)
        {
            if (inner_index[static_cast<std::size_t>(i)] != a.innerIndexPtr()[i])
            {
                return false;
            }
        }
        return true;
    }

    void analyze_if_needed(const SparseMatrix &a)
    {
        if (same_pattern(a))
        {
            return;
        }
        sparse_lu.analyzePattern(a);
        rows = a.rows();
        cols = a.cols();
        outer_index.assign(a.outerIndexPtr(), a.outerIndexPtr() + a.outerSize() + 1);
        inner_index.assign(a.innerIndexPtr(), a.innerIndexPtr() + a.nonZeros());
        analyzed = true;
    }
};

LinearSystemSolver::LinearSystemSolver()
    : impl_(std::make_unique<Impl>())
{
}

LinearSystemSolver::~LinearSystemSolver() = default;

void AssemblyAccumulator::scatter_add(std::vector<Triplet> &global, const std::vector<int> &edof, const Eigen::MatrixXd &local) const
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

void AssemblyAccumulator::scatter_add(Eigen::VectorXd &global, const std::vector<int> &edof, const Eigen::VectorXd &local) const
{
    for (int i = 0; i < static_cast<int>(edof.size()); ++i)
    {
        global(edof[static_cast<std::size_t>(i)]) += local(i);
    }
}

Eigen::VectorXd LinearSystemSolver::solve(SparseMatrix a,
                                          const Eigen::VectorXd &b,
                                          const SolverOptions &options)
{
    a.makeCompressed();
    constexpr int kDensePreferredMaxDofs = 500;
    if (a.rows() <= kDensePreferredMaxDofs && a.cols() <= kDensePreferredMaxDofs)
    {
        return solve_with_dense_fallback(a, b);
    }

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
        try
        {
            return solve_with_iterative_fallback(a, b, options.iterative_maxiter, options.iterative_tolerance);
        }
        catch (const std::runtime_error &)
        {
        }
    }

    impl_->analyze_if_needed(a);
    impl_->sparse_lu.factorize(a);
    if (impl_->sparse_lu.info() != Eigen::Success)
    {
        impl_->analyzed = false;
        impl_->analyze_if_needed(a);
        impl_->sparse_lu.factorize(a);
    }
    if (impl_->sparse_lu.info() == Eigen::Success)
    {
        const Eigen::VectorXd x = impl_->sparse_lu.solve(b);
        if (is_acceptable_solution(a, b, x))
        {
            return x;
        }
    }

    a.diagonal().array() += 1.0e-10;
    a.makeCompressed();
    impl_->analyze_if_needed(a);
    impl_->sparse_lu.factorize(a);
    if (impl_->sparse_lu.info() != Eigen::Success)
    {
        impl_->analyzed = false;
        impl_->analyze_if_needed(a);
        impl_->sparse_lu.factorize(a);
    }
    if (impl_->sparse_lu.info() != Eigen::Success)
    {
        impl_->analyzed = false;
        try
        {
            return solve_with_fresh_sparse_lu(a, b);
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
                if (a.rows() <= 2500)
                {
                    return solve_with_dense_fallback(a, b);
                }
                throw;
            }
        }
    }
    {
        const Eigen::VectorXd x = impl_->sparse_lu.solve(b);
        if (is_acceptable_solution(a, b, x))
        {
            return x;
        }
    }
    if (a.rows() <= 2500)
    {
        return solve_with_dense_fallback(a, b);
    }
    return solve_with_iterative_fallback(a,
                                         b,
                                         std::max(options.iterative_maxiter, 4000),
                                         std::min(options.iterative_tolerance, 1.0e-10));
}

void ConstraintHandler::apply_dirichlet(const BoundaryCondition &bc,
                                        const Eigen::VectorXd &uu,
                                        double scale,
                                        SparseMatrix &k,
                                        Eigen::VectorXd &rhs) const
{
    std::vector<char> constrained(static_cast<std::size_t>(k.cols()), 0);
    Eigen::VectorXd prescribed_vals = Eigen::VectorXd::Zero(k.cols());
    for (int i = 0; i < static_cast<int>(bc.dofs.size()); ++i)
    {
        const int dof = bc.dofs[static_cast<std::size_t>(i)];
        constrained[static_cast<std::size_t>(dof)] = 1;
        prescribed_vals(dof) = scale * bc.values(i) - uu(dof);
    }

    for (int row = 0; row < k.outerSize(); ++row)
    {
        for (SparseMatrix::InnerIterator it(k, row); it; ++it)
        {
            if (constrained[static_cast<std::size_t>(it.col())])
            {
                rhs(row) -= it.value() * prescribed_vals(it.col());
                if (row != it.col())
                {
                    it.valueRef() = 0.0;
                }
            }
        }
    }

    for (int i = 0; i < static_cast<int>(bc.dofs.size()); ++i)
    {
        const int dof = bc.dofs[static_cast<std::size_t>(i)];
        for (SparseMatrix::InnerIterator it(k, dof); it; ++it)
        {
            if (it.col() != dof)
            {
                it.valueRef() = 0.0;
            }
        }
        k.coeffRef(dof, dof) = 1.0;
        rhs(dof) = prescribed_vals(dof);
    }
    k.prune(0.0);
}

Result ResultBuilder::finalize(const Problem &problem,
                               const Eigen::VectorXd &uu,
                               double strain_energy,
                               const std::vector<double> &history,
                               bool compute_exact_solution) const
{
    Result result;
    result.uu = uu;
    result.nodal_u = GeometryUtils::nodal_displacements(uu);
    result.strain_energy = strain_energy;
    result.strain_energy_history = history;
    if (compute_exact_solution)
    {
        PatchTestField field;
        if (problem.problem_type == "bending_block")
        {
            result.exact_u = field.exact_bending_block(problem.mesh.nodes);
        }
        else
        {
            result.exact_u = field.exact_not_so_simple_shear(problem.mesh.nodes);
        }
        const Eigen::MatrixXd diff = result.exact_u - result.nodal_u;
        result.relative_error = std::sqrt(diff.array().square().sum() / result.exact_u.array().square().sum());
    }
    else
    {
        result.relative_error = std::numeric_limits<double>::quiet_NaN();
    }
    return result;
}

} // namespace nonlinear_esfem_t6::shared
