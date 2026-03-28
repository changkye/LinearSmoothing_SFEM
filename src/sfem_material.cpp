#include "sfem_material.hpp"

#include <cmath>

namespace nonlinear_esfem_t6::shared
{

NeoHookeanMaterialModel::NeoHookeanMaterialModel(const Eigen::Vector2d &parameters)
    : parameters_(parameters)
{
}

void NeoHookeanMaterialModel::initialize_bmat_workspace(int nnode,
                                                        Eigen::MatrixXd &bmat,
                                                        Eigen::MatrixXd &bgeo) const
{
    const int ncol = 2 * nnode;
    if (bmat.rows() != 3 || bmat.cols() != ncol)
    {
        bmat.resize(3, ncol);
    }
    if (bgeo.rows() != 4 || bgeo.cols() != ncol)
    {
        bgeo.resize(4, ncol);
    }
}

void NeoHookeanMaterialModel::nonlinear_bmat(const Eigen::MatrixXd &dndx,
                                             const Eigen::MatrixXd &wk_u,
                                             Eigen::MatrixXd &bmat,
                                             Eigen::MatrixXd &bgeo,
                                             Eigen::Matrix3d &fmat) const
{
    const Eigen::Matrix2d dxdx0 = wk_u.transpose() * dndx;
    fmat.setIdentity();
    fmat.block<2, 2>(0, 0) += dxdx0;

    const int nnode = static_cast<int>(dndx.rows());
    initialize_bmat_workspace(nnode, bmat, bgeo);
    bmat.setZero();
    bgeo.setZero();

    for (int a = 0; a < nnode; ++a)
    {
        const double dx = dndx(a, 0);
        const double dy = dndx(a, 1);
        bmat(0, 2 * a) = dx * fmat(0, 0);
        bmat(0, 2 * a + 1) = dx * fmat(1, 0);
        bmat(1, 2 * a) = dy * fmat(0, 1);
        bmat(1, 2 * a + 1) = dy * fmat(1, 1);
        bmat(2, 2 * a) = dy * fmat(0, 0) + dx * fmat(0, 1);
        bmat(2, 2 * a + 1) = dx * fmat(1, 1) + dy * fmat(1, 0);
        bgeo(0, 2 * a) = dx;
        bgeo(1, 2 * a) = dy;
        bgeo(2, 2 * a + 1) = dx;
        bgeo(3, 2 * a + 1) = dy;
    }
}

NonlinearMaterial NeoHookeanMaterialModel::constitutive(const Eigen::Matrix3d &fmat) const
{
    const Eigen::Matrix3d c = fmat.transpose() * fmat;
    const Eigen::Matrix3d inv_c = c.inverse();
    const double i1 = c.trace();
    const double i3 = c.determinant();

    const double lambda = parameters_(1) - (2.0 / 3.0) * parameters_(0);
    const double mu = parameters_(0) - 0.5 * lambda * std::log(i3);

    const Eigen::Matrix3d ident = Eigen::Matrix3d::Identity();
    const Eigen::Matrix3d s = 0.5 * lambda * std::log(i3) * inv_c + parameters_(0) * (ident - inv_c);

    NonlinearMaterial result;
    result.smat.block<2, 2>(0, 0) = s.block<2, 2>(0, 0);
    result.smat.block<2, 2>(2, 2) = s.block<2, 2>(0, 0);

    result.cmat(0, 0) = lambda * inv_c(0, 0) * inv_c(0, 0) + 2.0 * mu * inv_c(0, 0) * inv_c(0, 0);
    result.cmat(0, 1) = lambda * inv_c(0, 0) * inv_c(1, 1) + 2.0 * mu * inv_c(0, 1) * inv_c(0, 1);
    result.cmat(0, 2) = lambda * inv_c(0, 0) * inv_c(0, 1) + mu * (inv_c(0, 0) * inv_c(0, 1) + inv_c(0, 1) * inv_c(0, 0));
    result.cmat(1, 0) = result.cmat(0, 1);
    result.cmat(1, 1) = lambda * inv_c(1, 1) * inv_c(1, 1) + 2.0 * mu * inv_c(1, 1) * inv_c(1, 1);
    result.cmat(1, 2) = lambda * inv_c(1, 1) * inv_c(0, 1) + mu * (inv_c(1, 0) * inv_c(1, 1) + inv_c(1, 1) * inv_c(1, 0));
    result.cmat(2, 0) = result.cmat(0, 2);
    result.cmat(2, 1) = result.cmat(1, 2);
    result.cmat(2, 2) = lambda * inv_c(0, 1) * inv_c(0, 1) + mu * (inv_c(0, 0) * inv_c(1, 1) + inv_c(0, 1) * inv_c(1, 0));

    result.w0 = lambda * std::pow(std::log(i3), 2.0) / 8.0 - 0.5 * parameters_(0) * std::log(i3) + 0.5 * parameters_(0) * (i1 - 3.0);
    return result;
}

} // namespace nonlinear_esfem_t6::shared
