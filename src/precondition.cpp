#include "constants.h"
#include "boundary_element.h"

#include <Accelerate/Accelerate.h>

#include <iostream>

static int lu_decomp(double* A, int N, int* ipiv);
static void lu_solve(double* matrixA, int N, int* ipiv, double* rhs);

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
    const double* __restrict__ particles_x_ptr    = particles_.x_ptr();
    const double* __restrict__ particles_y_ptr    = particles_.y_ptr();
    const double* __restrict__ particles_z_ptr    = particles_.z_ptr();

    const double* __restrict__ particles_nx_ptr   = particles_.nx_ptr();
    const double* __restrict__ particles_ny_ptr   = particles_.ny_ptr();
    const double* __restrict__ particles_nz_ptr   = particles_.nz_ptr();
    const double* __restrict__ particles_area_ptr = particles_.area_ptr();

    double potential_coeff_1 = 0.5 * (1. +      params_.phys_eps_);
    double potential_coeff_2 = 0.5 * (1. + 1. / params_.phys_eps_);

    for (auto leaf_idx : tree_.leaves()) {

        auto particle_idxs = tree_.node_particle_idxs(leaf_idx);
        std::size_t particle_begin = particle_idxs[0];
        std::size_t particle_end   = particle_idxs[1];
        std::size_t num_particles = particle_end - particle_begin;
        std::size_t num_rows = 2 * num_particles;

        std::vector<double> column_major_A(num_rows * num_rows, 0.);
        std::vector<double> rhs(num_rows, 0.);
        std::vector<int> ipiv(num_rows, 0);

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

                    column_major_A[(row                ) * num_rows + (col                )] = -L1 * source_area;
                    column_major_A[(row                ) * num_rows + (col + num_particles)] = -L2 * source_area;
                    column_major_A[(row + num_particles) * num_rows + (col                )] = -L3 * source_area;
                    column_major_A[(row + num_particles) * num_rows + (col + num_particles)] = -L4 * source_area;

                    L1 = target_cos * tp1 * (1. - tp2 * eps);
                    L4 = source_cos * tp1 * (1. - tp2 / eps);

                    column_major_A[(col                ) * num_rows + (row                )] = -L1 * source_area;
                    column_major_A[(col                ) * num_rows + (row + num_particles)] = -L2 * source_area;
                    column_major_A[(col + num_particles) * num_rows + (row                )] = -L3 * source_area;
                    column_major_A[(col + num_particles) * num_rows + (row + num_particles)] = -L4 * source_area;
                }
            }

            column_major_A[(row                ) * num_rows + (row                )] = potential_coeff_1;
            column_major_A[(row + num_particles) * num_rows + (row + num_particles)] = potential_coeff_2;

            rhs[row]                 = r[j];
            rhs[row + num_particles] = r[j + num_total_particles];
        }

        int nrhs = 1, info;
        int num_rows_int = (int)num_rows;
        
        //CLAPACK style call:
        //dgesv_(&num_rows_int, &nrhs, column_major_A.data(), &num_rows_int,
        //      ipiv.data(), rhs.data(), &num_rows_int, &info);

        lu_decomp(column_major_A.data(), num_rows_int, ipiv.data());
        lu_solve(column_major_A.data(), num_rows_int, ipiv.data(), rhs.data());

        for (std::size_t j = particle_begin; j < particle_end; ++j) {
            z[j]                       = rhs[j - particle_begin];
            z[j + num_total_particles] = rhs[j - particle_begin + num_particles];
        }

    }

//    int i, j, idx = 0, nrow = 0, nrow2, ibeg = 0, iend = 0;
//    int *ipiv, inc, nrhs, info;
//    //double **matrixA;
//    double *columnMajorA, *rhs;
//    double L1, L2, L3, L4, area;
//    double tp[3], tq[3], sp[3], sq[3];
//    double r_s[3], rs, irs, sumrs;
//    double G0, kappa_rs, exp_kappa_rs, Gk;
//    double cos_theta, cos_theta0, tp1, tp2, dot_tqsq;
//    double G10, G20, G1, G2, G3, G4;
//    double pre1, pre2;
//
//    pre1 = 0.5*(1.0+s_eps);
//    pre2 = 0.5*(1.0+1.0/s_eps);
//
//    //make_matrix(matrixA, 2*s_max_per_leaf, 2*s_max_per_leaf);
//    make_vector(columnMajorA, 4*s_max_per_leaf*s_max_per_leaf);
//    make_vector(ipiv, 2*s_max_per_leaf);
//    make_vector(rhs, 2*s_max_per_leaf);
//
//    while (idx < s_numpars) {
//        leaflength(s_tree_root, idx, &nrow);
//        nrow2 = nrow*2;
//        ibeg  = idx;
//        iend  = idx + nrow - 1;
//
//        memset(columnMajorA, 0, nrow2*nrow2*sizeof(double));
//        memset(ipiv, 0, nrow2*sizeof(int));
//        memset(rhs, 0, nrow2*sizeof(double));
//
//        for (i = ibeg; i <= iend; i++) {
//            tp[0] = s_particle_position[0][i];
//            tp[1] = s_particle_position[1][i];
//            tp[2] = s_particle_position[2][i];
//            tq[0] = s_particle_normal[0][i];
//            tq[1] = s_particle_normal[1][i];
//            tq[2] = s_particle_normal[2][i];
//
//            for (j = ibeg; j < i; j++) {
//                sp[0] = s_particle_position[0][j];
//                sp[1] = s_particle_position[1][j];
//                sp[2] = s_particle_position[2][j];
//                sq[0] = s_particle_normal[0][j];
//                sq[1] = s_particle_normal[1][j];
//                sq[2] = s_particle_normal[2][j];
//
//                r_s[0] = sp[0]-tp[0]; r_s[1] = sp[1]-tp[1]; r_s[2] = sp[2]-tp[2];
//                sumrs = r_s[0]*r_s[0] + r_s[1]*r_s[1] + r_s[2]*r_s[2];
//                rs = sqrt(sumrs);
//                irs = 1.0/rs;
//                G0 = ONE_OVER_4PI * irs;
//                kappa_rs = s_kappa * rs;
//                exp_kappa_rs = exp(-kappa_rs);
//                Gk = exp_kappa_rs * G0;
//
//                cos_theta  = (sq[0]*r_s[0] + sq[1]*r_s[1] + sq[2]*r_s[2]) * irs;
//                cos_theta0 = (tq[0]*r_s[0] + tq[1]*r_s[1] + tq[2]*r_s[2]) * irs;
//                tp1 = G0* irs;
//                tp2 = (1.0 + kappa_rs) * exp_kappa_rs;
//
//                G10 = cos_theta0 * tp1;
//                G20 = tp2 * G10;
//
//                G1 = cos_theta * tp1;
//                G2 = tp2 * G1;
//
//                dot_tqsq = sq[0]*tq[0] + sq[1]*tq[1] + sq[2]*tq[2];
//                G3 = (dot_tqsq - 3.0*cos_theta0*cos_theta) * irs*tp1;
//                G4 = tp2*G3 - s_kappa2*cos_theta0*cos_theta*Gk;
//
//                area = s_particle_area[j];
//
//                L1 = G1 - s_eps*G2;
//                L2 = G0 - Gk;
//                L3 = G4 - G3;
//                L4 = G10 - G20/s_eps;
//
//                //matrixA[i-ibeg][j-ibeg] = -L1*area;
//                //matrixA[i-ibeg][j+nrow-ibeg] = -L2*area;
//                //matrixA[i+nrow-ibeg][j-ibeg] = -L3*area;
//                //matrixA[i+nrow-ibeg][j+nrow-ibeg] = -L4*area;
//                columnMajorA[(j-ibeg)*nrow2 + i-ibeg] = -L1*area;
//                columnMajorA[(j+nrow-ibeg)*nrow2 + i-ibeg] = -L2*area;
//                columnMajorA[(j-ibeg)*nrow2 + i+nrow-ibeg] = -L3*area;
//                columnMajorA[(j+nrow-ibeg)*nrow2 + i+nrow-ibeg] = -L4*area;
//            }
//
//            //matrixA[i-ibeg][i-ibeg] = pre1;
//            //matrixA[i+nrow-ibeg][i+nrow-ibeg] = pre2;
//            columnMajorA[(i-ibeg)*nrow2 + i-ibeg] = pre1;
//            columnMajorA[(i+nrow-ibeg)*nrow2 + i+nrow-ibeg] = pre2;
//
//            for (j = i+1; j <= iend; j++) {
//                sp[0] = s_particle_position[0][j];
//                sp[1] = s_particle_position[1][j];
//                sp[2] = s_particle_position[2][j];
//                sq[0] = s_particle_normal[0][j];
//                sq[1] = s_particle_normal[1][j];
//                sq[2] = s_particle_normal[2][j];
//
//                r_s[0] = sp[0]-tp[0]; r_s[1] = sp[1]-tp[1]; r_s[2] = sp[2]-tp[2];
//                sumrs = r_s[0]*r_s[0] + r_s[1]*r_s[1] + r_s[2]*r_s[2];
//                rs = sqrt(sumrs);
//                irs = 1.0/rs;
//                G0 = ONE_OVER_4PI * irs;
//                kappa_rs = s_kappa * rs;
//                exp_kappa_rs = exp(-kappa_rs);
//                Gk = exp_kappa_rs * G0;
//
//                cos_theta  = (sq[0]*r_s[0] + sq[1]*r_s[1] + sq[2]*r_s[2]) * irs;
//                cos_theta0 = (tq[0]*r_s[0] + tq[1]*r_s[1] + tq[2]*r_s[2]) * irs;
//                tp1 = G0 * irs;
//                tp2 = (1.0 + kappa_rs) * exp_kappa_rs;
//
//                G10 = cos_theta0 * tp1;
//                G20 = tp2 * G10;
//
//                G1 = cos_theta * tp1;
//                G2 = tp2 * G1;
//
//                dot_tqsq = sq[0]*tq[0] + sq[1]*tq[1] + sq[2]*tq[2];
//                G3 = (dot_tqsq - 3.0*cos_theta0*cos_theta) * irs*tp1;
//                G4 = tp2*G3 - s_kappa2*cos_theta0*cos_theta*Gk;
//
//                area = s_particle_area[j];
//
//                L1 = G1 - s_eps*G2;
//                L2 = G0 - Gk;
//                L3 = G4 - G3;
//                L4 = G10 - G20/s_eps;
//
//                //matrixA[i-ibeg][j-ibeg] = -L1*area;
//                //matrixA[i-ibeg][j+nrow-ibeg] = -L2*area;
//                //matrixA[i+nrow-ibeg][j-ibeg] = -L3*area;
//                //matrixA[i+nrow-ibeg][j+nrow-ibeg] = -L4*area;
//                columnMajorA[(j-ibeg)*nrow2 + i-ibeg] = -L1*area;
//                columnMajorA[(j+nrow-ibeg)*nrow2 + i-ibeg] = -L2*area;
//                columnMajorA[(j-ibeg)*nrow2 + i+nrow-ibeg] = -L3*area;
//                columnMajorA[(j+nrow-ibeg)*nrow2 + i+nrow-ibeg] = -L4*area;
//            }
//        }
//
//        for (i = 0; i < nrow; i++) {
//            rhs[i] = r[i+ibeg];
//            rhs[i+nrow] = r[i+ibeg+s_numpars];
//        }
//
//        // Jiahui's implementation of LU decomposition
//        //inc = lu_decomp(matrixA, nrow2, ipiv);
//        //lu_solve(matrixA, nrow2, ipiv, rhs);
//
//        // Apple Accelerate implementation of LAPACK LU decomposition
//        //nrhs = 1;
//        //dgesv_(&nrow2, &nrhs, columnMajorA, &nrow2, ipiv, rhs, &nrow2, &info);
//
//        // LAPACKE implementation of LAPACK LU decomposition
//        LAPACKE_dgesv(LAPACK_COL_MAJOR, nrow2, 1, columnMajorA, nrow2, ipiv, rhs, nrow2);
//
//        for (i = 0; i < nrow; i++) {
//            z[i+ibeg] = rhs[i];
//            z[i+ibeg+s_numpars] = rhs[i+nrow];
//        }
//
//        idx += nrow;
//    }
//
//    //free_matrix(matrixA);
//    free_vector(columnMajorA);
//    free_vector(rhs);
//    free_vector(ipiv);


    timers_.precondition.stop();
}



static int lu_decomp(double* A, int N, int* ipiv)
{
    int imax;
    double maxA, absA, Tol = 1.0e-14;

    for (int i = 0; i < N; i++)
        ipiv[i] = i; // record pivoting number

    for (int i = 0; i < N; i++) {
        maxA = 0.0;
        imax = i;

        for (int k = i; k < N; k++) {
            if ((absA = std::abs(A[k*N+i])) > maxA) {
                maxA = absA;
                imax = k;
            }
        }

        if (maxA < Tol) return 0; //failure, matrix is degenerate

        if (imax != i) {
            //pivoting P
            int idx = ipiv[i];
            ipiv[i] = ipiv[imax];
            ipiv[imax] = idx;

            //counting pivots starting from N (for determinant)
            ipiv[N]++;

            for (int kk = 0; kk < N; ++kk) {
                double temp    = A[i*N + kk];
                A[i*N + kk]    = A[imax*N + kk];
                A[imax*N + kk] = temp;
            }
        }

        for (int j = i + 1; j < N; j++) {
            A[j*N + i] /= A[i*N + i];

            for (int k = i + 1; k < N; k++) {
                A[j*N + k] -= A[j*N + i] * A[i*N + k];
            }
        }
  }

  return 1;
}

static void lu_solve(double* matrixA, int N, int* ipiv, double* rhs)
{
  /* b will contain the solution */
  
  int i, k;

  std::vector<double> xtemp(N, 0.);

  for (int i = 0; i < N; i++) {
    xtemp[i] = rhs[ipiv[i]];

    for (int k = 0; k < i; k++)
      xtemp[i] -= matrixA[i*N +k] * xtemp[k];
  }

  for (int i = N - 1; i >= 0; i--) {
    for (int k = i + 1; k < N; k++) {
      xtemp[i] -= matrixA[i*N + k] * xtemp[k];
    }

    xtemp[i] = xtemp[i] / matrixA[i*N + i];
  }

  for (int i = 0; i < N; i++) {
    rhs[i] = xtemp[i];
  }
  
}
