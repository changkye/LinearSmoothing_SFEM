#pragma once

#include "nonlinear_esfem_t6.hpp"

#include <filesystem>

namespace nonlinear_esfem_t6::shared
{

class VtuWriter
{
public:
    void write(const std::filesystem::path &file_path,
               const Mesh &mesh,
               const Eigen::MatrixXd &displacement,
               const Eigen::MatrixXd &exact_displacement) const;
};

} // namespace nonlinear_esfem_t6::shared
