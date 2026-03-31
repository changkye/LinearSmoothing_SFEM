#include "nonlinear_sfem_t6.hpp"
#include "sfem_shared.hpp"

namespace nonlinear_esfem_t6
{
Result NonlinearEsfemT6Solver::solve(const NonlinearSfemT6Problem &wrapped_problem,
                                     const SolverOptions &options) const
{
    const Problem &problem = wrapped_problem.data();
    const shared::EdgeSmoothedDomainAssembler assembler;
    const shared::NonlinearSmoothedFemSolver solver;
    return solver.solve(problem, options, assembler);
}

} // namespace nonlinear_esfem_t6
