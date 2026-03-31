#include "model_parameters.hpp"
#include "nonlinear_sfem_t6.hpp"
#include "sfem_shared.hpp"

namespace nonlinear_esfem_t6
{

    NonlinearSfemT6Problem::NonlinearSfemT6Problem(Problem problem)
        : problem_(std::move(problem))
    {
    }

    NonlinearSfemT6Problem NonlinearSfemT6Problem::make_not_so_simple_shear(
        const Eigen::Vector2i &num_els,
        const Eigen::Vector2d &material)
    {
        Problem problem;
        const auto geometry = model_parameters::geometry_for("nonlinear_patch");
        const Eigen::Matrix<double, 2, 2> limits =
            (Eigen::Matrix<double, 2, 2>() << 0.0, geometry.L, 0.0, geometry.H)
                .finished();
        const shared::T6MeshFactory mesh_factory;
        const shared::PatchTestField field;
        const shared::BoundaryConditionBuilder bc_builder;
        problem.mesh = mesh_factory.make_uniform(limits, num_els);
        problem.bc = bc_builder.make_from_exact_field(problem.mesh.nodes,
                                                      field.exact_not_so_simple_shear(problem.mesh.nodes));
        problem.force = Eigen::VectorXd::Zero(2 * problem.mesh.nodes.rows());
        problem.problem_type = "nonlinear_patch";
        problem.material = material;
        problem.num_els = num_els;
        problem.edge_boundary_ng = model_parameters::esfem_edge_boundary_ng;
        problem.edge_tri_quad_order = model_parameters::esfem_edge_tri_quad_order;
        problem.edge_quad_quad_order = model_parameters::esfem_edge_quad_quad_order;
        problem.edge_reg_param = model_parameters::esfem_edge_reg_param;
        return NonlinearSfemT6Problem(std::move(problem));
    }

    const Problem &NonlinearSfemT6Problem::data() const
    {
        return problem_;
    }

    Problem &NonlinearSfemT6Problem::data()
    {
        return problem_;
    }

    void write_vtu(const std::filesystem::path &file_path,
                   const Mesh &mesh,
                   const Eigen::MatrixXd &displacement,
                   const Eigen::MatrixXd &exact_displacement)
    {
        const shared::VtuWriter writer;
        writer.write(file_path, mesh, displacement, exact_displacement);
    }

} // namespace nonlinear_esfem_t6
