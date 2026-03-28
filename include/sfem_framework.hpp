#pragma once

#include "nonlinear_esfem_t6.hpp"
#include "sfem_material.hpp"

#include <memory>
#include <vector>

namespace nonlinear_esfem_t6::shared
{

struct SupportDomain
{
    std::vector<int> nodes;
    std::vector<int> edof;
    Eigen::MatrixXd dx;
    Eigen::MatrixXd dy;
    std::vector<double> weights;
    double sub_area = 0.0;
};

class ISmoothingDomainAssembler
{
public:
    virtual ~ISmoothingDomainAssembler() = default;
    virtual std::vector<SupportDomain> build_domains(const Problem &problem) const = 0;
};

class NonlinearSmoothedFemSolver
{
public:
    Result solve(const Problem &problem,
                 const SolverOptions &options,
                 const ISmoothingDomainAssembler &assembler) const;
};

} // namespace nonlinear_esfem_t6::shared
