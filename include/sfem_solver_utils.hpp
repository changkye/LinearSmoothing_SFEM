#pragma once

#include "nonlinear_sfem_t6.hpp"

#include <Eigen/Dense>
#include <Eigen/Sparse>

#include <memory>
#include <vector>

namespace nonlinear_esfem_t6::shared
{

using SparseMatrix = Eigen::SparseMatrix<double, Eigen::RowMajor>;
using Triplet = Eigen::Triplet<double>;

class AssemblyAccumulator
{
public:
    void scatter_add(std::vector<Triplet> &global, const std::vector<int> &edof, const Eigen::MatrixXd &local) const;
    void scatter_add(Eigen::VectorXd &global, const std::vector<int> &edof, const Eigen::VectorXd &local) const;
};

class LinearSystemSolver
{
public:
    LinearSystemSolver();
    ~LinearSystemSolver();

    Eigen::VectorXd solve(SparseMatrix a,
                          const Eigen::VectorXd &b,
                          const SolverOptions &options);

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

class ConstraintHandler
{
public:
    void apply_dirichlet(const BoundaryCondition &bc,
                         const Eigen::VectorXd &uu,
                         double scale,
                         SparseMatrix &k,
                         Eigen::VectorXd &rhs) const;
};

int regularize_near_zero_rows(SparseMatrix &a, Eigen::VectorXd &b);

class ResultBuilder
{
public:
    Result finalize(const Problem &problem,
                    const Eigen::VectorXd &uu,
                    double strain_energy,
                    const std::vector<double> &history,
                    bool compute_exact_solution) const;
};

} // namespace nonlinear_esfem_t6::shared
