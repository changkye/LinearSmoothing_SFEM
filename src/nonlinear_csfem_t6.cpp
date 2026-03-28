#include "nonlinear_esfem_t6.hpp"
#include "sfem_shared.hpp"

namespace nonlinear_esfem_t6
{

Result NonlinearCsfemT6Solver::solve(const NonlinearEsfemT6Problem &wrapped_problem,
                                     const SolverOptions &options) const
{
    const Problem &problem = wrapped_problem.data();
    const shared::CellSmoothedDomainAssembler assembler;
    const shared::NonlinearSmoothedFemSolver solver;
    return solver.solve(problem, options, assembler);
}

} // namespace nonlinear_esfem_t6
