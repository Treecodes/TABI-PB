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
    
    for (std::size_t i = 0;                i <     particles_.num(); ++i) z[i] = r[i] / potential_coeff_1;
    for (std::size_t i = particles_.num(); i < 2 * particles_.num(); ++i) z[i] = r[i] / potential_coeff_2;

    timers_.precondition.stop();
}


void BoundaryElement::precondition_block(double *z, double *r)
{
    timers_.precondition.start();

    double eps    = params_.phys_eps_;
    double kappa  = params_.phys_kappa_;
    double kappa2 = params_.phys_kappa2_;

    const std::size_t num_total_particles         = particles_.num();
    const double* __restrict particles_x_ptr    = particles_.x_ptr();
    const double* __restrict particles_y_ptr    = particles_.y_ptr();
    const double* __restrict particles_z_ptr    = particles_.z_ptr();

    const double* __restrict particles_nx_ptr   = particles_.nx_ptr();
    const double* __restrict particles_ny_ptr   = particles_.ny_ptr();
    const double* __restrict particles_nz_ptr   = particles_.nz_ptr();
    const double* __restrict particles_area_ptr = particles_.area_ptr();

    double potential_coeff_1 = 0.5 * (1. +      params_.phys_eps_);
    double potential_coeff_2 = 0.5 * (1. + 1. / params_.phys_eps_);

    for (auto leaf_idx : tree_.leaves()) {

        auto particle_idxs = tree_.node_particle_idxs(leaf_idx);
        std::size_t particle_begin = particle_idxs[0];
        std::size_t particle_end   = particle_idxs[1];
        std::size_t num_particles = particle_end - particle_begin;
        std::size_t num_cols = 2 * num_particles;

        std::vector<double> A(num_cols * num_cols, 0.);
        std::vector<double> rhs(num_cols, 0.);
        std::vector<int> pivot(num_cols, 0);

        for (std::size_t j = particle_begin; j < particle_end; ++j) {

            std::size_t row = j - particle_begin;

            double target_x = particles_x_ptr[j];
            double target_y = particles_y_ptr[j];
            double target_z = particles_z_ptr[j];

            double target_nx = particles_nx_ptr[j];
            double target_ny = particles_ny_ptr[j];
            double target_nz = particles_nz_ptr[j];

            for (std::size_t k = particle_begin; k < j; ++k) {

                std::size_t col = k - particle_begin;

                double source_x = particles_x_ptr[k];
                double source_y = particles_y_ptr[k];
                double source_z = particles_z_ptr[k];

                double source_nx = particles_nx_ptr[k];
                double source_ny = particles_ny_ptr[k];
                double source_nz = particles_nz_ptr[k];
                double source_area = particles_area_ptr[k];

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

                    A[(row                ) * num_cols + (col                )] = -L1 * source_area;
                    A[(row                ) * num_cols + (col + num_particles)] = -L2 * source_area;
                    A[(row + num_particles) * num_cols + (col                )] = -L3 * source_area;
                    A[(row + num_particles) * num_cols + (col + num_particles)] = -L4 * source_area;

                    L1 = target_cos * tp1 * (1. - tp2 * eps);
                    L4 = source_cos * tp1 * (1. - tp2 / eps);

                    A[(col                ) * num_cols + (row                )] = -L1 * source_area;
                    A[(col                ) * num_cols + (row + num_particles)] = -L2 * source_area;
                    A[(col + num_particles) * num_cols + (row                )] = -L3 * source_area;
                    A[(col + num_particles) * num_cols + (row + num_particles)] = -L4 * source_area;
                }
            }

            A[(row                ) * num_cols + (row                )] = potential_coeff_1;
            A[(row + num_particles) * num_cols + (row + num_particles)] = potential_coeff_2;

            rhs[row]                 = r[j];
            rhs[row + num_particles] = r[j + num_total_particles];
        }

        int nrhs = 1, info;
        int num_cols_int = (int)num_cols;
        
        //CLAPACK style call:
        //dgesv_(&num_cols_int, &nrhs, column_major_A.data(), &num_cols_int,
        //       pivot.data(), rhs.data(), &num_cols_int, &info);

        lu_decomp(A.data(), num_cols_int, pivot.data());
        lu_solve(A.data(), num_cols_int, pivot.data(), rhs.data());

        for (std::size_t j = particle_begin; j < particle_end; ++j) {
            z[j]                       = rhs[j - particle_begin];
            z[j + num_total_particles] = rhs[j - particle_begin + num_particles];
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
