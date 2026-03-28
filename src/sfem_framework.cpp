#include "sfem_framework.hpp"
#include "sfem_assemblers.hpp"
#include "sfem_material.hpp"
#include "sfem_solver_utils.hpp"

#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace nonlinear_esfem_t6::shared
{

namespace
{

bool sparse_matrix_all_finite(const SparseMatrix &k)
{
    for (int outer = 0; outer < k.outerSize(); ++outer)
    {
        for (SparseMatrix::InnerIterator it(k, outer); it; ++it)
        {
            if (!std::isfinite(it.value()))
            {
                return false;
            }
        }
    }
    return true;
}

void assemble_internal_force_and_energy(const std::vector<SupportDomain> &domains,
                                        const NeoHookeanMaterialModel &material_model,
                                        const Eigen::VectorXd &uu,
                                        Eigen::VectorXd &internal_force,
                                        double &energy)
{
    internal_force.setZero();
    energy = 0.0;

    for (const SupportDomain &domain : domains)
    {
        Eigen::MatrixXd wk_u(domain.nodes.size(), 2);
        Eigen::MatrixXd dndx(domain.nodes.size(), 2);
        Eigen::MatrixXd bmat;
        Eigen::MatrixXd bgeo;
        material_model.initialize_bmat_workspace(static_cast<int>(domain.nodes.size()), bmat, bgeo);
        for (int a = 0; a < static_cast<int>(domain.nodes.size()); ++a)
        {
            wk_u.row(a) << uu(2 * domain.nodes[static_cast<std::size_t>(a)]),
                uu(2 * domain.nodes[static_cast<std::size_t>(a)] + 1);
        }

        Eigen::VectorXd re = Eigen::VectorXd::Zero(domain.edof.size());
        double w0_acc = 0.0;
        double w_acc = 0.0;
        for (int ig = 0; ig < domain.dx.rows(); ++ig)
        {
            dndx.col(0) = domain.dx.row(ig).transpose();
            dndx.col(1) = domain.dy.row(ig).transpose();

            Eigen::Matrix3d fmat;
            material_model.nonlinear_bmat(dndx, wk_u, bmat, bgeo, fmat);
            const NonlinearMaterial material = material_model.constitutive(fmat);
            const double weight = domain.weights[static_cast<std::size_t>(ig)];
            const Eigen::Vector3d stress(material.smat(0, 0), material.smat(1, 1), material.smat(0, 1));
            re += bmat.transpose() * stress * weight;
            w0_acc += material.w0 * weight;
            w_acc += weight;
        }

        for (int i = 0; i < static_cast<int>(domain.edof.size()); ++i)
        {
            internal_force(domain.edof[static_cast<std::size_t>(i)]) += re(i);
        }
        if (w_acc > 0.0)
        {
            energy += (w0_acc / w_acc) * domain.sub_area;
        }
    }
}

void apply_dirichlet_to_rhs_only(const BoundaryCondition &bc,
                                 const Eigen::VectorXd &uu,
                                 double scale,
                                 Eigen::VectorXd &rhs)
{
    for (int i = 0; i < static_cast<int>(bc.dofs.size()); ++i)
    {
        const int dof = bc.dofs[static_cast<std::size_t>(i)];
        rhs(dof) = scale * bc.values(i) - uu(dof);
    }
}

void write_matrix_m(std::ofstream &out, const char *name, const Eigen::MatrixXd &mat)
{
    out << name << " = [\n";
    for (int r = 0; r < mat.rows(); ++r)
    {
        out << "  ";
        for (int c = 0; c < mat.cols(); ++c)
        {
            out << mat(r, c);
            if (c + 1 < mat.cols())
            {
                out << ", ";
            }
        }
        out << ";\n";
    }
    out << "];\n";
}

void write_worst_cell_domain_bundle(const Problem &problem,
                                    const std::vector<SupportDomain> &domains,
                                    const NeoHookeanMaterialModel &material_model,
                                    const Eigen::VectorXd &uu,
                                    const SolverOptions &options,
                                    int step,
                                    int iter)
{
    if (!options.debug_csfem_bending || options.debug_output_dir.empty())
    {
        return;
    }
    std::filesystem::create_directories(options.debug_output_dir);
    const std::filesystem::path file =
        options.debug_output_dir /
        ("csfem_worst_domain_" + std::to_string(problem.num_els(0)) + "x" +
         std::to_string(problem.num_els(1)) + "_step" + std::to_string(step) +
         "_iter" + std::to_string(iter) + ".m");
    std::ofstream out(file);
    if (!out)
    {
        return;
    }

    double worst_eig = std::numeric_limits<double>::infinity();
    int worst_domain = -1;
    Eigen::MatrixXd worst_ke;
    Eigen::MatrixXd worst_re;
    Eigen::MatrixXd worst_wku;
    Eigen::MatrixXd worst_detc;
    Eigen::MatrixXd worst_fnorm;

    for (int idomain = 0; idomain < static_cast<int>(domains.size()); ++idomain)
    {
        const SupportDomain &domain = domains[static_cast<std::size_t>(idomain)];
        Eigen::MatrixXd wk_u(domain.nodes.size(), 2);
        for (int a = 0; a < static_cast<int>(domain.nodes.size()); ++a)
        {
            wk_u.row(a) << uu(2 * domain.nodes[static_cast<std::size_t>(a)]),
                uu(2 * domain.nodes[static_cast<std::size_t>(a)] + 1);
        }

        Eigen::MatrixXd ke = Eigen::MatrixXd::Zero(domain.edof.size(), domain.edof.size());
        Eigen::VectorXd re = Eigen::VectorXd::Zero(domain.edof.size());
        Eigen::MatrixXd dndx(domain.nodes.size(), 2);
        Eigen::MatrixXd detc(domain.dx.rows(), 1);
        Eigen::MatrixXd fnorm(domain.dx.rows(), 1);
        Eigen::MatrixXd bmat;
        Eigen::MatrixXd bgeo;
        material_model.initialize_bmat_workspace(static_cast<int>(domain.nodes.size()), bmat, bgeo);
        for (int ig = 0; ig < domain.dx.rows(); ++ig)
        {
            dndx.col(0) = domain.dx.row(ig).transpose();
            dndx.col(1) = domain.dy.row(ig).transpose();
            Eigen::Matrix3d fmat;
            material_model.nonlinear_bmat(dndx, wk_u, bmat, bgeo, fmat);
            const Eigen::Matrix2d cmat2 = fmat.topLeftCorner<2, 2>().transpose() * fmat.topLeftCorner<2, 2>();
            detc(ig, 0) = cmat2.determinant();
            fnorm(ig, 0) = fmat.norm();
            const NonlinearMaterial material = material_model.constitutive(fmat);
            const double weight = domain.weights[static_cast<std::size_t>(ig)];
            ke += (bmat.transpose() * material.cmat * bmat + bgeo.transpose() * material.smat * bgeo) * weight;
            const Eigen::Vector3d stress(material.smat(0, 0), material.smat(1, 1), material.smat(0, 1));
            re += bmat.transpose() * stress * weight;
        }

        const Eigen::MatrixXd ke_sym = 0.5 * (ke + ke.transpose());
        const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(ke_sym);
        const double eig_min = eig.info() == Eigen::Success ? eig.eigenvalues().minCoeff() : std::numeric_limits<double>::quiet_NaN();
        if (eig_min < worst_eig)
        {
            worst_eig = eig_min;
            worst_domain = idomain;
            worst_ke = ke;
            worst_re = re;
            worst_wku = wk_u;
            worst_detc = detc;
            worst_fnorm = fnorm;
        }
    }

    if (worst_domain < 0)
    {
        return;
    }

    const SupportDomain &domain = domains[static_cast<std::size_t>(worst_domain)];
    Eigen::MatrixXd coords(domain.nodes.size(), 2);
    for (int i = 0; i < static_cast<int>(domain.nodes.size()); ++i)
    {
        coords.row(i) = problem.mesh.nodes.row(domain.nodes[static_cast<std::size_t>(i)]);
    }

    out << std::setprecision(16);
    out << "debug.step = " << step << ";\n";
    out << "debug.iter = " << iter << ";\n";
    out << "debug.worst_domain = " << worst_domain + 1 << ";\n";
    out << "debug.worst_eig_min = " << worst_eig << ";\n";
    out << "debug.nodes = [";
    for (std::size_t i = 0; i < domain.nodes.size(); ++i)
    {
        out << domain.nodes[i] + 1;
        out << (i + 1 == domain.nodes.size() ? "" : ", ");
    }
    out << "];\n";
    out << "debug.edof = [";
    for (std::size_t i = 0; i < domain.edof.size(); ++i)
    {
        out << domain.edof[i] + 1;
        out << (i + 1 == domain.edof.size() ? "" : ", ");
    }
    out << "];\n";
    write_matrix_m(out, "debug.coords", coords);
    write_matrix_m(out, "debug.wk_u", worst_wku);
    write_matrix_m(out, "debug.dx", domain.dx);
    write_matrix_m(out, "debug.dy", domain.dy);
    write_matrix_m(out, "debug.ke", worst_ke);
    write_matrix_m(out, "debug.re", worst_re);
    write_matrix_m(out, "debug.detc", worst_detc);
    write_matrix_m(out, "debug.fnorm", worst_fnorm);
    out << "debug.weights = [";
    for (std::size_t i = 0; i < domain.weights.size(); ++i)
    {
        out << domain.weights[i];
        out << (i + 1 == domain.weights.size() ? "" : ", ");
    }
    out << "];\n";
}

void write_cell_domain_trace(const Problem &problem,
                             const std::vector<SupportDomain> &domains,
                             const NeoHookeanMaterialModel &material_model,
                             const Eigen::VectorXd &uu,
                             const SolverOptions &options,
                             int step,
                             int iter)
{
    if (!options.debug_csfem_bending || options.debug_output_dir.empty())
    {
        return;
    }
    std::filesystem::create_directories(options.debug_output_dir);
    const std::filesystem::path file =
        options.debug_output_dir /
        ("csfem_cell_trace_" + std::to_string(problem.num_els(0)) + "x" +
         std::to_string(problem.num_els(1)) + "_step" + std::to_string(step) +
         "_iter" + std::to_string(iter) + ".txt");
    std::ofstream out(file);
    if (!out)
    {
        return;
    }

    out << std::setprecision(16);
    out << "num_domains " << domains.size() << '\n';
    out << "num_nodes " << problem.mesh.nodes.rows() << '\n';
    out << "num_elements_x " << problem.num_els(0) << '\n';
    out << "num_elements_y " << problem.num_els(1) << '\n';

    std::vector<int> node_counts(static_cast<std::size_t>(problem.mesh.nodes.rows()), 0);
    for (const SupportDomain &domain : domains)
    {
        for (const int node : domain.nodes)
        {
            ++node_counts[static_cast<std::size_t>(node)];
        }
    }

    out << "node_participation\n";
    for (int node = 0; node < static_cast<int>(node_counts.size()); ++node)
    {
        out << node + 1 << ' ' << node_counts[static_cast<std::size_t>(node)] << '\n';
    }

    out << "domain_stats\n";
    for (int idomain = 0; idomain < static_cast<int>(domains.size()); ++idomain)
    {
        const SupportDomain &domain = domains[static_cast<std::size_t>(idomain)];
        Eigen::MatrixXd wk_u(domain.nodes.size(), 2);
        for (int a = 0; a < static_cast<int>(domain.nodes.size()); ++a)
        {
            wk_u.row(a) << uu(2 * domain.nodes[static_cast<std::size_t>(a)]),
                uu(2 * domain.nodes[static_cast<std::size_t>(a)] + 1);
        }

        Eigen::MatrixXd ke = Eigen::MatrixXd::Zero(domain.edof.size(), domain.edof.size());
        Eigen::MatrixXd dndx(domain.nodes.size(), 2);
        Eigen::MatrixXd bmat;
        Eigen::MatrixXd bgeo;
        material_model.initialize_bmat_workspace(static_cast<int>(domain.nodes.size()), bmat, bgeo);
        double detc_min = std::numeric_limits<double>::infinity();
        double detc_max = -std::numeric_limits<double>::infinity();
        double f_norm_max = 0.0;
        bool has_nonfinite_state = false;
        for (int ig = 0; ig < domain.dx.rows(); ++ig)
        {
            dndx.col(0) = domain.dx.row(ig).transpose();
            dndx.col(1) = domain.dy.row(ig).transpose();
            Eigen::Matrix3d fmat;
            material_model.nonlinear_bmat(dndx, wk_u, bmat, bgeo, fmat);
            const Eigen::Matrix2d cmat2 = fmat.topLeftCorner<2, 2>().transpose() * fmat.topLeftCorner<2, 2>();
            const double detc = cmat2.determinant();
            detc_min = std::min(detc_min, detc);
            detc_max = std::max(detc_max, detc);
            f_norm_max = std::max(f_norm_max, fmat.norm());
            if (!fmat.allFinite() || !std::isfinite(detc))
            {
                has_nonfinite_state = true;
            }
            const NonlinearMaterial material = material_model.constitutive(fmat);
            const double weight = domain.weights[static_cast<std::size_t>(ig)];
            ke += (bmat.transpose() * material.cmat * bmat + bgeo.transpose() * material.smat * bgeo) * weight;
        }

        const Eigen::MatrixXd ke_sym = 0.5 * (ke + ke.transpose());
        const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(ke_sym);
        const double eig_min = eig.info() == Eigen::Success ? eig.eigenvalues().minCoeff() : std::numeric_limits<double>::quiet_NaN();
        const double eig_max = eig.info() == Eigen::Success ? eig.eigenvalues().maxCoeff() : std::numeric_limits<double>::quiet_NaN();

        out << idomain + 1
            << " dx_max " << domain.dx.cwiseAbs().maxCoeff()
            << " dy_max " << domain.dy.cwiseAbs().maxCoeff()
            << " uu_norm " << wk_u.norm()
            << " detc_min " << detc_min
            << " detc_max " << detc_max
            << " f_norm_max " << f_norm_max
            << " has_nonfinite_state " << (has_nonfinite_state ? 1 : 0)
            << " ke_diag_min " << ke.diagonal().minCoeff()
            << " ke_diag_max " << ke.diagonal().maxCoeff()
            << " eig_min " << eig_min
            << " eig_max " << eig_max
            << '\n';
    }
}

void append_global_matrix_trace(const Problem &problem,
                                const BoundaryCondition &bc,
                                const SparseMatrix &k,
                                const Eigen::VectorXd &rhs,
                                const SolverOptions &options,
                                int step,
                                int iter,
                                const char *label)
{
    if (!options.debug_csfem_bending || options.debug_output_dir.empty())
    {
        return;
    }
    const std::filesystem::path file =
        options.debug_output_dir /
        ("csfem_cell_trace_" + std::to_string(problem.num_els(0)) + "x" +
         std::to_string(problem.num_els(1)) + "_step" + std::to_string(step) +
         "_iter" + std::to_string(iter) + ".txt");
    std::ofstream out(file, std::ios::app);
    if (!out)
    {
        return;
    }

    const Eigen::VectorXd diag = k.diagonal();
    out << "global_matrix " << label << '\n';
    out << "rhs_norm " << rhs.norm() << '\n';
    out << "diag_abs_min " << diag.cwiseAbs().minCoeff() << '\n';
    out << "diag_abs_max " << diag.cwiseAbs().maxCoeff() << '\n';
    if (k.rows() <= 500)
    {
        Eigen::MatrixXd dense = Eigen::MatrixXd(k);
        const Eigen::MatrixXd sym = 0.5 * (dense + dense.transpose());
        const Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> eig(sym);
        if (eig.info() == Eigen::Success)
        {
            out << "global_eig_min " << eig.eigenvalues().minCoeff() << '\n';
            out << "global_eig_max " << eig.eigenvalues().maxCoeff() << '\n';
            out << "global_abs_eig_min " << eig.eigenvalues().cwiseAbs().minCoeff() << '\n';
        }
    }
    out << "near_zero_diag_rows\n";
    for (int i = 0; i < diag.size(); ++i)
    {
        if (std::abs(diag(i)) < 1.0e-10)
        {
            int nnz = 0;
            for (SparseMatrix::InnerIterator it(k, i); it; ++it)
            {
                ++nnz;
            }
            out << i + 1 << ' ' << diag(i) << ' ' << nnz << '\n';
        }
    }
    out << "constrained_rows\n";
    for (int dof : bc.dofs)
    {
        int nnz = 0;
        double max_offdiag = 0.0;
        for (SparseMatrix::InnerIterator it(k, dof); it; ++it)
        {
            ++nnz;
            if (it.col() != dof)
            {
                max_offdiag = std::max(max_offdiag, std::abs(it.value()));
            }
        }
        out << dof + 1 << ' ' << diag(dof) << ' ' << nnz << ' ' << max_offdiag << '\n';
    }
}

void append_update_trace(const Problem &problem,
                         const BoundaryCondition &bc,
                         const Eigen::VectorXd &uu,
                         const Eigen::VectorXd &du,
                         const SolverOptions &options,
                         int step,
                         int iter)
{
    if (!options.debug_csfem_bending || options.debug_output_dir.empty())
    {
        return;
    }
    const std::filesystem::path file =
        options.debug_output_dir /
        ("csfem_cell_trace_" + std::to_string(problem.num_els(0)) + "x" +
         std::to_string(problem.num_els(1)) + "_step" + std::to_string(step) +
         "_iter" + std::to_string(iter) + ".txt");
    std::ofstream out(file, std::ios::app);
    if (!out)
    {
        return;
    }

    Eigen::Index max_idx = 0;
    const double max_du = du.cwiseAbs().maxCoeff(&max_idx);
    bool constrained = false;
    double prescribed_value = 0.0;
    for (int i = 0; i < static_cast<int>(bc.dofs.size()); ++i)
    {
        if (bc.dofs[static_cast<std::size_t>(i)] == static_cast<int>(max_idx))
        {
            constrained = true;
            prescribed_value = bc.values(i);
            break;
        }
    }

    out << "update_trace\n";
    out << "uu_norm_after " << (uu + du).norm() << '\n';
    out << "du_norm " << du.norm() << '\n';
    out << "du_abs_max " << max_du << '\n';
    out << "du_abs_max_dof " << max_idx + 1 << '\n';
    out << "du_abs_max_node " << (max_idx / 2) + 1 << '\n';
    out << "du_abs_max_component " << ((max_idx % 2) == 0 ? "ux" : "uy") << '\n';
    out << "du_abs_max_constrained " << (constrained ? 1 : 0) << '\n';
    out << "du_abs_max_bc_value " << prescribed_value << '\n';
    if ((max_idx / 2) < problem.mesh.nodes.rows())
    {
        out << "du_abs_max_node_x " << problem.mesh.nodes(max_idx / 2, 0) << '\n';
        out << "du_abs_max_node_y " << problem.mesh.nodes(max_idx / 2, 1) << '\n';
    }
}

void append_solver_comparison_trace(const Problem &problem,
                                    const SparseMatrix &k,
                                    const Eigen::VectorXd &rhs,
                                    const Eigen::VectorXd &du_sparse,
                                    const SolverOptions &options,
                                    int step,
                                    int iter)
{
    if (!options.debug_csfem_bending || options.debug_output_dir.empty() || k.rows() > 500)
    {
        return;
    }

    const std::filesystem::path file =
        options.debug_output_dir /
        ("csfem_cell_trace_" + std::to_string(problem.num_els(0)) + "x" +
         std::to_string(problem.num_els(1)) + "_step" + std::to_string(step) +
         "_iter" + std::to_string(iter) + ".txt");
    std::ofstream out(file, std::ios::app);
    if (!out)
    {
        return;
    }

    const Eigen::MatrixXd dense = Eigen::MatrixXd(k);
    const Eigen::VectorXd du_dense = dense.fullPivLu().solve(rhs);
    const Eigen::VectorXd res_sparse = k * du_sparse - rhs;
    const Eigen::VectorXd res_dense = k * du_dense - rhs;
    Eigen::Index max_idx = 0;
    du_sparse.cwiseAbs().maxCoeff(&max_idx);

    out << "solver_compare\n";
    out << "sparse_res_norm " << res_sparse.norm() << '\n';
    out << "dense_res_norm " << res_dense.norm() << '\n';
    out << "sparse_dense_max_diff " << (du_sparse - du_dense).cwiseAbs().maxCoeff() << '\n';
    out << "tracked_dof " << max_idx + 1 << '\n';
    out << "tracked_sparse " << du_sparse(max_idx) << '\n';
    out << "tracked_dense " << du_dense(max_idx) << '\n';
}

} // namespace

Result NonlinearSmoothedFemSolver::solve(const Problem &problem,
                                         const SolverOptions &options,
                                         const ISmoothingDomainAssembler &assembler) const
{
    const NeoHookeanMaterialModel material_model(problem.material);
    const AssemblyAccumulator accumulator;
    const ConstraintHandler constraints;
    LinearSystemSolver linear_solver;
    const ResultBuilder result_builder;
    const std::vector<SupportDomain> domains = assembler.build_domains(problem);

    const int ndof = 2 * problem.mesh.nodes.rows();
    const std::size_t triplet_count_estimate = [&]() {
        std::size_t count = 0;
        for (const SupportDomain &domain : domains)
        {
            count += domain.edof.size() * domain.edof.size();
        }
        return count;
    }();
    Eigen::VectorXd uu = Eigen::VectorXd::Zero(ndof);
    std::vector<double> history;
    history.reserve(options.nstep);
    bool late_instability_dumped = false;

    double last_energy = 0.0;
    const double initial_increment = 1.0 / static_cast<double>(std::max(1, options.nstep));
    const double min_increment = initial_increment / std::pow(2.0, std::max(0, options.max_step_cuts));
    double current_load_factor = 0.0;
    double load_increment = initial_increment;
    int step = 0;
    int step_cut_count = 0;

    while (current_load_factor < 1.0 - 1.0e-12)
    {
        const double target_load_factor = std::min(1.0, current_load_factor + load_increment);
        const Eigen::VectorXd uu_step_start = uu;
        double condition = 1.0;
        double residual = std::numeric_limits<double>::infinity();
        double best_residual = std::numeric_limits<double>::infinity();
        int plateau_count = 0;
        int iter = 0;
        bool step_converged = false;
        bool step_stagnated = false;

        std::cout << "\n Step " << step + 1 << "\t Scale " << (target_load_factor * 100.0) << "%\n";

        while (iter < options.maxiter)
        {
            ++iter;
            const auto iter_begin = std::chrono::steady_clock::now();
            double energy = 0.0;

            if (step == 1 && iter <= 2 &&
                dynamic_cast<const CellSmoothedDomainAssembler *>(&assembler) != nullptr)
            {
                write_cell_domain_trace(problem, domains, material_model, uu, options, step, iter);
                if (iter == 2)
                {
                    write_worst_cell_domain_bundle(problem, domains, material_model, uu, options, step, iter);
                }
            }

            int nthreads = 1;
#ifdef _OPENMP
            nthreads = omp_get_max_threads();
#endif
            std::vector<std::vector<Triplet>> triplet_buffers(static_cast<std::size_t>(nthreads));
            std::vector<Eigen::VectorXd> residual_buffers(static_cast<std::size_t>(nthreads), Eigen::VectorXd::Zero(ndof));
            for (auto &buffer : triplet_buffers)
            {
                buffer.reserve(triplet_count_estimate / static_cast<std::size_t>(nthreads) + 256);
            }

            #pragma omp parallel for if(domains.size() > 1) default(none) shared(domains, uu, material_model, accumulator, triplet_buffers, residual_buffers, nthreads) reduction(+:energy)
            for (int idomain = 0; idomain < static_cast<int>(domains.size()); ++idomain)
            {
#ifdef _OPENMP
                const int tid = omp_get_thread_num();
#else
                const int tid = 0;
#endif
                const SupportDomain &domain = domains[static_cast<std::size_t>(idomain)];
                Eigen::MatrixXd wk_u(domain.nodes.size(), 2);
                Eigen::MatrixXd dndx(domain.nodes.size(), 2);
                Eigen::MatrixXd bmat;
                Eigen::MatrixXd bgeo;
                material_model.initialize_bmat_workspace(static_cast<int>(domain.nodes.size()), bmat, bgeo);
                for (int a = 0; a < static_cast<int>(domain.nodes.size()); ++a)
                {
                    wk_u.row(a) << uu(2 * domain.nodes[static_cast<std::size_t>(a)]), uu(2 * domain.nodes[static_cast<std::size_t>(a)] + 1);
                }

                Eigen::MatrixXd ke = Eigen::MatrixXd::Zero(domain.edof.size(), domain.edof.size());
                Eigen::VectorXd re = Eigen::VectorXd::Zero(domain.edof.size());
                double w0_acc = 0.0;
                double w_acc = 0.0;

                for (int ig = 0; ig < domain.dx.rows(); ++ig)
                {
                    dndx.col(0) = domain.dx.row(ig).transpose();
                    dndx.col(1) = domain.dy.row(ig).transpose();

                    Eigen::Matrix3d fmat;
                    material_model.nonlinear_bmat(dndx, wk_u, bmat, bgeo, fmat);
                    const NonlinearMaterial material = material_model.constitutive(fmat);
                    const double weight = domain.weights[static_cast<std::size_t>(ig)];

                    ke += (bmat.transpose() * material.cmat * bmat + bgeo.transpose() * material.smat * bgeo) * weight;
                    const Eigen::Vector3d stress(material.smat(0, 0), material.smat(1, 1), material.smat(0, 1));
                    re += bmat.transpose() * stress * weight;
                    w0_acc += material.w0 * weight;
                    w_acc += weight;
                }

                if (w_acc > 0.0)
                {
                    energy += (w0_acc / w_acc) * domain.sub_area;
                }

                accumulator.scatter_add(triplet_buffers[static_cast<std::size_t>(tid)], domain.edof, ke);
                accumulator.scatter_add(residual_buffers[static_cast<std::size_t>(tid)], domain.edof, re);
            }

            std::vector<Triplet> triplets;
            triplets.reserve(triplet_count_estimate);
            Eigen::VectorXd r = Eigen::VectorXd::Zero(ndof);
            for (int tid = 0; tid < nthreads; ++tid)
            {
                auto &buffer = triplet_buffers[static_cast<std::size_t>(tid)];
                triplets.insert(triplets.end(), buffer.begin(), buffer.end());
                r += residual_buffers[static_cast<std::size_t>(tid)];
            }
            const auto assembly_end = std::chrono::steady_clock::now();

            SparseMatrix k(ndof, ndof);
            k.setFromTriplets(triplets.begin(), triplets.end(), [](double lhs, double rhs) { return lhs + rhs; });
            k.makeCompressed();
            const auto matrix_end = std::chrono::steady_clock::now();
            Eigen::VectorXd rhs = target_load_factor * problem.force - r;
            constraints.apply_dirichlet(problem.bc, uu, target_load_factor, k, rhs);
            if (step == 1 && iter <= 2 &&
                dynamic_cast<const CellSmoothedDomainAssembler *>(&assembler) != nullptr)
            {
                append_global_matrix_trace(problem, problem.bc, k, rhs, options, step, iter, "before_solve");
            }
            const auto bc_end = std::chrono::steady_clock::now();
            residual = rhs.norm() / static_cast<double>(ndof);
            const double previous_best_residual = best_residual;
            best_residual = std::min(best_residual, residual);
            if (options.aggressive_stagnation_control &&
                std::isfinite(previous_best_residual) &&
                residual <= std::max(options.residual_tolerance * 500.0, 1.0e-4) &&
                residual >= previous_best_residual * 0.98)
            {
                ++plateau_count;
            }
            else
            {
                plateau_count = 0;
            }
            if (dynamic_cast<const CellSmoothedDomainAssembler *>(&assembler) != nullptr &&
                (!rhs.allFinite() || !std::isfinite(residual) || !sparse_matrix_all_finite(k)))
            {
                write_cell_domain_trace(problem, domains, material_model, uu, options, step, iter);
                write_worst_cell_domain_bundle(problem, domains, material_model, uu, options, step, iter);
                append_global_matrix_trace(problem, problem.bc, k, rhs, options, step, iter, "nonfinite_before_solve");
                throw std::runtime_error("Non-finite CS-FEM system detected before solve");
            }
            Eigen::VectorXd du;
            try
            {
                du = linear_solver.solve(k, rhs, options);
            }
            catch (const std::runtime_error &)
            {
                if (step == 1 && iter <= 2 &&
                    dynamic_cast<const CellSmoothedDomainAssembler *>(&assembler) != nullptr)
                {
                    append_global_matrix_trace(problem, problem.bc, k, rhs, options, step, iter, "solve_failure");
                }
                throw;
            }
            const auto solve_end = std::chrono::steady_clock::now();
            if (!du.allFinite())
            {
                if (dynamic_cast<const CellSmoothedDomainAssembler *>(&assembler) != nullptr)
                {
                    write_cell_domain_trace(problem, domains, material_model, uu, options, step, iter);
                    write_worst_cell_domain_bundle(problem, domains, material_model, uu, options, step, iter);
                }
                break;
            }
            if (step == 1 && iter == 1 &&
                dynamic_cast<const CellSmoothedDomainAssembler *>(&assembler) != nullptr)
            {
                append_solver_comparison_trace(problem, k, rhs, du, options, step, iter);
            }

            double accepted_alpha = 1.0;
            Eigen::VectorXd accepted_du = du;
            double accepted_energy = energy;
            if (options.line_search_max_backtracks > 0)
            {
                const double current_rhs_norm = rhs.norm();
                Eigen::VectorXd trial_internal = Eigen::VectorXd::Zero(ndof);
                double trial_energy = energy;
                double alpha = 1.0;
                bool accepted = false;
                for (int iback = 0; iback <= options.line_search_max_backtracks; ++iback)
                {
                    const Eigen::VectorXd trial_uu = uu + alpha * du;
                    assemble_internal_force_and_energy(domains, material_model, trial_uu, trial_internal, trial_energy);
                    Eigen::VectorXd trial_rhs = target_load_factor * problem.force - trial_internal;
                    apply_dirichlet_to_rhs_only(problem.bc, trial_uu, target_load_factor, trial_rhs);
                    const double trial_rhs_norm = trial_rhs.norm();
                    if (trial_rhs.allFinite() && std::isfinite(trial_rhs_norm) &&
                        trial_rhs_norm <= current_rhs_norm)
                    {
                        accepted_alpha = alpha;
                        accepted_du = alpha * du;
                        accepted_energy = trial_energy;
                        accepted = true;
                        break;
                    }
                    alpha *= options.line_search_reduction;
                    if (alpha < options.line_search_min_alpha)
                    {
                        break;
                    }
                }
                if (!accepted)
                {
                    accepted_alpha = std::max(options.line_search_min_alpha, std::pow(options.line_search_reduction, options.line_search_max_backtracks));
                    accepted_du = accepted_alpha * du;
                    Eigen::VectorXd trial_internal = Eigen::VectorXd::Zero(ndof);
                    assemble_internal_force_and_energy(domains, material_model, uu + accepted_du, trial_internal, accepted_energy);
                }
            }

            uu += accepted_du;
            if (step == 1 && iter <= 2 &&
                dynamic_cast<const CellSmoothedDomainAssembler *>(&assembler) != nullptr)
            {
                append_update_trace(problem, problem.bc, uu - accepted_du, accepted_du, options, step, iter);
            }
            condition = accepted_du.norm() / std::max(uu.norm(), std::numeric_limits<double>::epsilon());
            if (!std::isfinite(condition))
            {
                if (dynamic_cast<const CellSmoothedDomainAssembler *>(&assembler) != nullptr)
                {
                    write_cell_domain_trace(problem, domains, material_model, uu, options, step, iter);
                    write_worst_cell_domain_bundle(problem, domains, material_model, uu, options, step, iter);
                }
                break;
            }
            if (dynamic_cast<const CellSmoothedDomainAssembler *>(&assembler) != nullptr &&
                !late_instability_dumped &&
                step > 1 &&
                (condition > 5.0 || residual > 10.0 || iter == options.maxiter))
            {
                write_cell_domain_trace(problem, domains, material_model, uu, options, step, iter);
                write_worst_cell_domain_bundle(problem, domains, material_model, uu, options, step, iter);
                append_global_matrix_trace(problem, problem.bc, k, rhs, options, step, iter, "late_instability");
                late_instability_dumped = true;
            }
            last_energy = accepted_energy;

            std::cout << std::fixed << std::setprecision(6)
                      << "Iteration Number " << iter
                      << " Condition " << condition
                      << " Residual " << residual
                      << " Tolerance " << options.tolerance;
            if (accepted_alpha != 1.0)
            {
                std::cout << " Alpha " << accepted_alpha;
            }
            std::cout << '\n';
            std::cout << "Timing [ms] assembly="
                      << std::chrono::duration<double, std::milli>(assembly_end - iter_begin).count()
                      << " sparse_build="
                      << std::chrono::duration<double, std::milli>(matrix_end - assembly_end).count()
                      << " bc="
                      << std::chrono::duration<double, std::milli>(bc_end - matrix_end).count()
                      << " solve="
                      << std::chrono::duration<double, std::milli>(solve_end - bc_end).count()
                      << " total="
                      << std::chrono::duration<double, std::milli>(solve_end - iter_begin).count()
                      << '\n';

            if (condition <= options.tolerance && residual <= options.residual_tolerance)
            {
                step_converged = true;
                break;
            }

            const bool tiny_alpha = accepted_alpha <= options.line_search_min_alpha * 1.0001;
            const bool residual_stagnating = residual > std::max(options.residual_tolerance * 50.0, best_residual * 1.05);
            const bool many_iterations = iter >= std::min(options.maxiter, 20);
            const bool residual_plateau =
                options.aggressive_stagnation_control &&
                iter >= std::min(options.maxiter, 12) &&
                plateau_count >= 6 &&
                accepted_alpha <= 0.125;
            const bool hard_stall =
                iter >= std::min(options.maxiter, 12) &&
                residual > std::max(options.residual_tolerance * 100.0, 1.0e-5) &&
                tiny_alpha;
            if (options.adaptive_load_stepping &&
                (hard_stall || residual_plateau || (many_iterations && residual_stagnating && tiny_alpha)))
            {
                step_stagnated = true;
                std::cout << "Step stagnation detected. Cutting load increment.\n";
                break;
            }
        }

        if (step_converged)
        {
            current_load_factor = target_load_factor;
            history.push_back(last_energy);
            ++step;
            step_cut_count = 0;
            if (options.adaptive_load_stepping &&
                options.allow_step_growth &&
                load_increment < initial_increment)
            {
                load_increment = std::min(initial_increment, 2.0 * load_increment);
            }
            continue;
        }

        uu = uu_step_start;
        if (!options.adaptive_load_stepping || load_increment <= min_increment + 1.0e-15)
        {
            history.push_back(last_energy);
            break;
        }

        load_increment *= 0.5;
        ++step_cut_count;
        std::cout << (step_stagnated ? "Stagnated step cut applied. " : "Step cut applied. ")
                  << "Retrying with smaller increment (" << (load_increment * 100.0) << "% load).\n";
    }

    if (current_load_factor < 1.0 - 1.0e-12)
    {
        std::cout << "Warning: nonlinear solve stopped at "
                  << std::fixed << std::setprecision(6) << (current_load_factor * 100.0)
                  << "% load. Final result corresponds to the last converged increment only.\n";
    }

    return result_builder.finalize(problem, uu, last_energy, history, options.compute_exact_solution);
}

} // namespace nonlinear_esfem_t6::shared
