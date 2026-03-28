#pragma once

#include "structural_solver.hpp"

namespace fem
{

class ProblemLibrary
{
public:
    static std::vector<Scenario> available();
    static std::string name(Scenario scenario);
    static Scenario parse(const std::string &value);
    static Eigen::MatrixXd exact_field(Scenario scenario, const Eigen::MatrixXd &nodes);
    static BoundaryCondition make_boundary_condition(Scenario scenario, const Eigen::MatrixXd &nodes);
    static Model build_model(Method method, Scenario scenario, const Eigen::Vector2i &num_els);
};

} // namespace fem
