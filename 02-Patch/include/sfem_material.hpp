#pragma once

#include <Eigen/Dense>

namespace nonlinear_esfem_t6::shared
{

struct NonlinearMaterial
{
    Eigen::Matrix3d cmat = Eigen::Matrix3d::Zero();
    Eigen::Matrix4d smat = Eigen::Matrix4d::Zero();
    double w0 = 0.0;
};

class NeoHookeanMaterialModel
{
public:
    explicit NeoHookeanMaterialModel(const Eigen::Vector2d &parameters);

    void initialize_bmat_workspace(int nnode, Eigen::MatrixXd &bmat, Eigen::MatrixXd &bgeo) const;
    void nonlinear_bmat(const Eigen::MatrixXd &dndx,
                        const Eigen::MatrixXd &wk_u,
                        Eigen::MatrixXd &bmat,
                        Eigen::MatrixXd &bgeo,
                        Eigen::Matrix3d &fmat) const;
    NonlinearMaterial constitutive(const Eigen::Matrix3d &fmat) const;

private:
    Eigen::Vector2d parameters_;
};

} // namespace nonlinear_esfem_t6::shared
