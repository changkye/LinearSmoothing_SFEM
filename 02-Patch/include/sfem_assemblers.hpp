#pragma once

#include "sfem_framework.hpp"

namespace nonlinear_esfem_t6::shared
{

class EdgeSmoothedDomainAssembler final : public ISmoothingDomainAssembler
{
public:
    std::vector<SupportDomain> build_domains(const Problem &problem) const override;
};

class CellSmoothedDomainAssembler final : public ISmoothingDomainAssembler
{
public:
    std::vector<SupportDomain> build_domains(const Problem &problem) const override;
};

} // namespace nonlinear_esfem_t6::shared
