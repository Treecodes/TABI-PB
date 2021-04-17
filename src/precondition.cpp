#include <cmath>

#include "constants.h"
#include "boundary_element.h"

static int lu_decomp(double* A, int N, int* pivot);
static void lu_solve(double* A, int N, int* pivot, double* rhs);


void BoundaryElement::precondition_diagonal(double *z, double *r)
{
    timers_.precondition.start();

    double potential_coeff_1 = 0.5 * (1. +      params_.phys_eps_);
    double potential_coeff_2 = 0.5 * (1. + 1. / params_.phys_eps_);
    
    for (std::size_t i = 0;                i <     elements_.num(); ++i) z[i] = r[i] / potential_coeff_1;
    for (std::size_t i = elements_.num(); i < 2 * elements_.num(); ++i) z[i] = r[i] / potential_coeff_2;

    timers_.precondition.stop();
}


void BoundaryElement::precondition_block(double *z, double *r)
{
    timers_.precondition.start();

    double eps    = params_.phys_eps_;
    double kappa  = params_.phys_kappa_;
    double kappa2 = params_.phys_kappa2_;

    const std::size_t num_total_elements         = elements_.num();
    const double* __restrict elements_x_ptr    = elements_.x_ptr();
    const double* __restrict elements_y_ptr    = elements_.y_ptr();
    const double* __restrict elements_z_ptr    = elements_.z_ptr();

    const double* __restrict elements_nx_ptr   = elements_.nx_ptr();
    const double* __restrict elements_ny_ptr   = elements_.ny_ptr();
    const double* __restrict elements_nz_ptr   = elements_.nz_ptr();
    const double* __restrict elements_area_ptr = elements_.area_ptr();

    double potential_coeff_1 = 0.5 * (1. +      params_.phys_eps_);
    double potential_coeff_2 = 0.5 * (1. + 1. / params_.phys_eps_);

    for (auto leaf_idx : tree_.leaves()) {

        auto element_idxs = tree_.node_particle_idxs(leaf_idx);
        std::size_t element_begin = element_idxs[0];
        std::size_t element_end   = element_idxs[1];
        std::size_t num_elements = element_end - element_begin;
        std::size_t num_cols = 2 * num_elements;

        std::vector<double> A(num_cols * num_cols, 0.);
        std::vector<double> rhs(num_cols, 0.);
        std::vector<int> pivot(num_cols, 0);

        for (std::size_t j = element_begin; j < element_end; ++j) {

            std::size_t row = j - element_begin;

            double target_x = elements_x_ptr[j];
            double target_y = elements_y_ptr[j];
            double target_z = elements_z_ptr[j];

            double target_nx = elements_nx_ptr[j];
            double target_ny = elements_ny_ptr[j];
            double target_nz = elements_nz_ptr[j];

            for (std::size_t k = element_begin; k < j; ++k) {

                std::size_t col = k - element_begin;

                double source_x = elements_x_ptr[k];
                double source_y = elements_y_ptr[k];
                double source_z = elements_z_ptr[k];

                double source_nx = elements_nx_ptr[k];
                double source_ny = elements_ny_ptr[k];
                double source_nz = elements_nz_ptr[k];
                double source_area = elements_area_ptr[k];

                double dist_x = source_x - target_x;
                double dist_y = source_y - target_y;
                double dist_z = source_z - target_z;
                double r = std::sqrt(dist_x * dist_x + dist_y * dist_y + dist_z * dist_z);

                if (r > 0) {
                    double one_over_r = 1. / r;
                    double G0 = constants::ONE_OVER_4PI * one_over_r;
                    double kappa_r = kappa * r;
                    double exp_kappa_r = std::exp(-kappa_r);
                    double Gk = exp_kappa_r * G0;
                    double source_cos = (source_nx * dist_x + source_ny * dist_y + source_nz * dist_z) * one_over_r;
                    double target_cos = (target_nx * dist_x + target_ny * dist_y + target_nz * dist_z) * one_over_r;
                    double tp1 = G0 * one_over_r;
                    double tp2 = (1. + kappa_r) * exp_kappa_r;

                    double dot_tqsq = source_nx * target_nx + source_ny * target_ny + source_nz * target_nz;
                    double G3 = (dot_tqsq - 3. * target_cos * source_cos) * one_over_r * tp1;
                    double G4 = tp2 * G3 - kappa2 * target_cos * source_cos * Gk;

                    double L1 = source_cos * tp1 * (1. - tp2 * eps);
                    double L2 = G0 - Gk;
                    double L3 = G4 - G3;
                    double L4 = target_cos * tp1 * (1. - tp2 / eps);

                    A[(row               ) * num_cols + (col               )] = -L1 * source_area;
                    A[(row               ) * num_cols + (col + num_elements)] = -L2 * source_area;
                    A[(row + num_elements) * num_cols + (col               )] = -L3 * source_area;
                    A[(row + num_elements) * num_cols + (col + num_elements)] = -L4 * source_area;

                    L1 = target_cos * tp1 * (1. - tp2 * eps);
                    L4 = source_cos * tp1 * (1. - tp2 / eps);

                    A[(col               ) * num_cols + (row               )] = -L1 * source_area;
                    A[(col               ) * num_cols + (row + num_elements)] = -L2 * source_area;
                    A[(col + num_elements) * num_cols + (row               )] = -L3 * source_area;
                    A[(col + num_elements) * num_cols + (row + num_elements)] = -L4 * source_area;
                }
            }

            A[(row               ) * num_cols + (row               )] = potential_coeff_1;
            A[(row + num_elements) * num_cols + (row + num_elements)] = potential_coeff_2;

            rhs[row]                = r[j];
            rhs[row + num_elements] = r[j + num_total_elements];
        }

        int nrhs = 1, info;
        int num_cols_int = (int)num_cols;
        
        //CLAPACK style call:
        //dgesv_(&num_cols_int, &nrhs, column_major_A.data(), &num_cols_int,
        //       pivot.data(), rhs.data(), &num_cols_int, &info);

        lu_decomp(A.data(), num_cols_int, pivot.data());
        lu_solve(A.data(), num_cols_int, pivot.data(), rhs.data());

        for (std::size_t j = element_begin; j < element_end; ++j) {
            z[j]                      = rhs[j - element_begin];
            z[j + num_total_elements] = rhs[j - element_begin + num_elements];
        }

    }

    timers_.precondition.stop();
}


static int lu_decomp(double* A, int N, int* pivot)
{
    // record pivoting number
    for (int i = 0; i < N; ++i) pivot[i] = i;

    for (int i = 0; i < N; ++i) {
        double A_max = 0.0;
        int idx_max = i;

        for (int k = i; k < N; ++k) {
            double A_abs = std::abs(A[k*N+i]);
            if (A_abs > A_max) {
                A_max = A_abs;
                idx_max = k;
            }
        }

        //failure, matrix is degenerate
        if (A_max < 1.e-14) return 1;

        if (idx_max != i) {
            //pivoting
            int idx = pivot[i];
            pivot[i] = pivot[idx_max];
            pivot[idx_max] = idx;

            //counting pivots starting from N (for determinant)
            pivot[N]++;

            for (int kk = 0; kk < N; ++kk) {
                double temp    = A[i*N + kk];
                A[i*N + kk]    = A[idx_max*N + kk];
                A[idx_max*N + kk] = temp;
            }
        }

        for (int j = i + 1; j < N; ++j) {
            A[j*N + i] /= A[i*N + i];

            for (int k = i + 1; k < N; ++k) {
                A[j*N + k] -= A[j*N + i] * A[i*N + k];
            }
        }
    }

    return 0;
}


static void lu_solve(double* A, int N, int* pivot, double* rhs)
{
    std::vector<double> xtemp(N, 0.);

    for (int i = 0; i < N; ++i) {
        xtemp[i] = rhs[pivot[i]];

        for (int k = 0; k < i; ++k) {
            xtemp[i] -= A[i*N +k] * xtemp[k];
        }
    }

    for (int i = N - 1; i >= 0; --i) {
        for (int k = i + 1; k < N; ++k) {
            xtemp[i] -= A[i*N + k] * xtemp[k];
        }

        xtemp[i] = xtemp[i] / A[i*N + i];
    }

    for (int i = 0; i < N; ++i) {
        rhs[i] = xtemp[i];
    }
}
