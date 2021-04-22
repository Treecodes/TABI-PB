#include <iomanip>
#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>
#include <iostream>
#include <numeric>
#include <array>
#include <cmath>
#include <cstdlib>
#include <cstdio>

#include "constants.h"
#include "elements.h"


static double triangle_area(std::array<std::array<double, 3>, 3> v);

 
Elements::Elements(const class Molecule& mol, const struct Params& params, struct Timers_Elements& timers)
    : Particles(params), molecule_(mol), timers_(timers)
{
    timers_.ctor.start();

    Elements::generate_elements(params_.mesh_, params_.mesh_density_, params_.mesh_probe_radius_);
    
    source_charge_.assign(num_, 0.);
    source_charge_dx_.assign(num_, 0.);
    source_charge_dy_.assign(num_, 0.);
    source_charge_dz_.assign(num_, 0.);
    
    target_charge_.assign(num_, 0.);
    target_charge_dx_.assign(num_, 0.);
    target_charge_dy_.assign(num_, 0.);
    target_charge_dz_.assign(num_, 0.);
    
    source_term_.assign(num_ * 2, 0.);
    
    order_.resize(num_);
    std::iota(order_.begin(), order_.end(), 0);

    timers_.ctor.stop();
}


void Elements::generate_elements(Params::Mesh mesh, double mesh_density, double probe_radius)
{
    std::ofstream NS_param_file("surfaceConfiguration.prm");
    
    NS_param_file << "Grid_scale = "                                 << mesh_density     << std::endl;
    NS_param_file << "Grid_perfil = "                                << 90.0             << std::endl;
    NS_param_file << "XYZR_FileName = "                              << "molecule.xyzr"  << std::endl;
    NS_param_file << "Build_epsilon_maps = "                         << "false"          << std::endl;
    NS_param_file << "Build_status_map = "                           << "false"          << std::endl;
    NS_param_file << "Save_Mesh_MSMS_Format = "                      << "true"           << std::endl;
    NS_param_file << "Compute_Vertex_Normals = "                     << "true"           << std::endl;

    if (mesh == Params::Mesh::SES)  NS_param_file                    << "Surface = ses"  << std::endl;
    if (mesh == Params::Mesh::SKIN) NS_param_file                    << "Surface = skin" << std::endl;
    
    NS_param_file << "Smooth_Mesh = "                                << "true"           << std::endl;
    NS_param_file << "Skin_Surface_Parameter = "                     << 0.45             << std::endl;
    NS_param_file << "Cavity_Detection_Filling = "                   << "false"          << std::endl;
    NS_param_file << "Conditional_Volume_Filling_Value = "           << 11.4             << std::endl;
    NS_param_file << "Keep_Water_Shaped_Cavities = "                 << "false"          << std::endl;
    NS_param_file << "Probe_Radius = "                               << probe_radius     << std::endl;
    NS_param_file << "Accurate_Triangulation = "                     << "true"           << std::endl;
    NS_param_file << "Triangulation = "                              << "true"           << std::endl;
    NS_param_file << "Check_duplicated_vertices = "                  << "true"           << std::endl;
    NS_param_file << "Save_Status_map = "                            << "false"          << std::endl;
    NS_param_file << "Save_PovRay = "                                << "false"          << std::endl;
    NS_param_file << "Max_ses_patches_per_auxiliary_grid_2d_cell = " << 1600             << std::endl;
    NS_param_file << "Max_ses_patches_auxiliary_grid_2d_size = "     << 50               << std::endl;

    NS_param_file.close();
    
#ifdef _WIN32
    std::system("NanoShaper.exe");
#else
    std::system("NanoShaper");
#endif

    std::remove("stderror.txt");
    std::remove("surfaceConfiguration.prm");
    std::remove("triangleAreas.txt");
    std::remove("exposed.xyz");
    std::remove("exposedIndices.txt");
    
    std::string line;
    
    // Read in the vert file
    std::ifstream vert_file("triangulatedSurf.vert");
    // Throw away first two lines
    std::getline(vert_file, line);
    std::getline(vert_file, line);
    
    std::getline(vert_file, line);
    num_ = std::stoul(line);

    x_.reserve(num_);
    y_.reserve(num_);
    z_.reserve(num_);
    nx_.reserve(num_);
    ny_.reserve(num_);
    nz_.reserve(num_);
    
    while (std::getline(vert_file, line)) {
        
        std::istringstream iss(line);
        std::vector<std::string> tokenized_line{std::istream_iterator<std::string> {iss},
                                                std::istream_iterator<std::string> {} };
                                                
        x_.push_back(std::stod(tokenized_line[0]));
        y_.push_back(std::stod(tokenized_line[1]));
        z_.push_back(std::stod(tokenized_line[2]));
        nx_.push_back(std::stod(tokenized_line[3]));
        ny_.push_back(std::stod(tokenized_line[4]));
        nz_.push_back(std::stod(tokenized_line[5]));
    }
    
    vert_file.close();
    
    // Read in the face file
    std::ifstream face_file("triangulatedSurf.face");
    // Throw away first two lines
    std::getline(face_file, line);
    std::getline(face_file, line);
    
    std::getline(face_file, line);
    num_faces_ = std::stoul(line);
    
    face_x_.reserve(num_faces_);
    face_y_.reserve(num_faces_);
    face_z_.reserve(num_faces_);
    
    while (std::getline(face_file, line)) {
        
        std::istringstream iss(line);
        std::vector<std::string> tokenized_line{std::istream_iterator<std::string> {iss},
                                                std::istream_iterator<std::string> {} };
                                                
        face_x_.push_back(std::stoul(tokenized_line[0]));
        face_y_.push_back(std::stoul(tokenized_line[1]));
        face_z_.push_back(std::stoul(tokenized_line[2]));
    }
    
    face_file.close();
    
    std::remove("molecule.xyzr");
    std::remove("triangulatedSurf.vert");
    std::remove("triangulatedSurf.face");

    area_.assign(num_, 0.);
    
    for (std::size_t i = 0; i < num_faces_; ++i) {
    
        std::array<std::size_t, 3> iface {face_x_[i], face_y_[i], face_z_[i]};
        std::array<std::array<double, 3>, 3> r; 
        
        for (int ii = 0; ii < 3; ++ii) {
            r[0][ii] = x_[iface[ii]-1];
            r[1][ii] = y_[iface[ii]-1];
            r[2][ii] = z_[iface[ii]-1];
        }

        for (int j = 0; j < 3; ++j) {
            area_[iface[j]-1] += triangle_area(r);
        }
    }
    
    std::transform(area_.begin(), area_.end(), area_.begin(),
                   [=](double x) { return x / 3.; } );
    surface_area_ = std::accumulate(area_.begin(), area_.end(), decltype(area_)::value_type(0));
    std::cout << "Surface area of triangulated mesh is " << surface_area_ 
              << ". " << std::endl << std::endl;
}



void Elements::compute_source_term()
{
/* this computes the source term where
 * S1=sum(qk*G0)/e1 S2=sim(qk*G0')/e1 */
    timers_.compute_source_term.start();

    double eps_solute = params_.phys_eps_solute_;

    std::size_t num_atoms = molecule_.num();
    std::size_t num = num_;
    
    const double* __restrict elements_x_ptr = x_.data();
    const double* __restrict elements_y_ptr = y_.data();
    const double* __restrict elements_z_ptr = z_.data();
    
    const double* __restrict elements_nx_ptr = nx_.data();
    const double* __restrict elements_ny_ptr = ny_.data();
    const double* __restrict elements_nz_ptr = nz_.data();
    
    const double* __restrict molecule_x_ptr = molecule_.x_ptr();
    const double* __restrict molecule_y_ptr = molecule_.y_ptr();
    const double* __restrict molecule_z_ptr = molecule_.z_ptr();
    const double* __restrict molecule_charge_ptr = molecule_.charge_ptr();
    
    double* __restrict elements_source_term_ptr = source_term_.data();

#ifdef OPENACC_ENABLED
    #pragma acc parallel loop gang present(molecule_x_ptr, molecule_y_ptr, \
                                      molecule_z_ptr, molecule_charge_ptr, \
                                      elements_x_ptr, elements_y_ptr, elements_z_ptr, \
                                      elements_nx_ptr, elements_ny_ptr, elements_nz_ptr, \
                                      elements_source_term_ptr)
#elif OPENMP_ENABLED
    #pragma omp parallel for
#endif
    for (std::size_t i = 0; i < num; ++i) {

        double source_term_1 = 0.;
        double source_term_2 = 0.;

#ifdef OPENACC_ENABLED
        #pragma acc loop vector reduction(+:source_term_1,source_term_2)
#endif
        for (std::size_t j = 0; j < num_atoms; ++j) {

  /* r_s = distance of charge position to triangular */
            double x_dist = molecule_x_ptr[j] - elements_x_ptr[i];
            double y_dist = molecule_y_ptr[j] - elements_y_ptr[i];
            double z_dist = molecule_z_ptr[j] - elements_z_ptr[i];
            double dist   = std::sqrt(x_dist*x_dist + y_dist*y_dist + z_dist*z_dist);

  /* cos_theta = <tr_q,r_s>/||r_s||_2 */
            double cos_theta = (elements_nx_ptr[i] * x_dist
                              + elements_ny_ptr[i] * y_dist
                              + elements_nz_ptr[i] * z_dist) / dist;

  /* G0 = 1/(4pi*||r_s||_2) */
            double G0 = constants::ONE_OVER_4PI / dist;

  /* G1 = cos_theta*G0/||r_s||_2 */
            double G1 = cos_theta * G0 / dist;

  /* update source term */
            source_term_1 += molecule_charge_ptr[j] * G0 / eps_solute;
            source_term_2 += molecule_charge_ptr[j] * G1 / eps_solute;
        }

        elements_source_term_ptr[i]       += source_term_1;
        elements_source_term_ptr[num + i] += source_term_2;
    }
    
    Elements::update_source_term_on_host();

    timers_.compute_source_term.stop();
}


void Elements::compute_source_term(const class InterpolationPoints& elem_interp_pts,
                                   const class Tree& elem_tree,
                                   const class Molecule& molecule,
                                   const class InterpolationPoints& mol_interp_pts,
                                   const class Tree& mol_tree,
                                   const class InteractionList& interaction_list)
{
/* this computes the source term where
 * S1=sum(qk*G0)/e1 S2=sim(qk*G0')/e1 */
    timers_.compute_source_term.start();
    
    class SourceTermCompute source_term(source_term_,
                                        *this,    elem_interp_pts, elem_tree,
                                        molecule, mol_interp_pts,  mol_tree,
                                        interaction_list, params_.phys_eps_solute_);
    
    source_term.compute();
    Elements::update_source_term_on_host();

    timers_.compute_source_term.stop();
}


void Elements::reorder()
{
    apply_order(order_.begin(), order_.end(), nx_.begin());
    apply_order(order_.begin(), order_.end(), ny_.begin());
    apply_order(order_.begin(), order_.end(), nz_.begin());
    
    apply_order(order_.begin(), order_.end(), area_.begin());
    apply_order(order_.begin(), order_.end(), source_term_.begin());
    apply_order(order_.begin(), order_.end(), source_term_.begin() + num_);
}


void Elements::unorder()
{
    apply_unorder(order_.begin(), order_.end(), x_.begin());
    apply_unorder(order_.begin(), order_.end(), y_.begin());
    apply_unorder(order_.begin(), order_.end(), z_.begin());
    
    apply_unorder(order_.begin(), order_.end(), nx_.begin());
    apply_unorder(order_.begin(), order_.end(), ny_.begin());
    apply_unorder(order_.begin(), order_.end(), nz_.begin());
    
    apply_unorder(order_.begin(), order_.end(), area_.begin());
    apply_unorder(order_.begin(), order_.end(), source_term_.begin());
    apply_unorder(order_.begin(), order_.end(), source_term_.begin() + num_);
}


void Elements::unorder(std::vector<double>& potential)
{
    apply_unorder(order_.begin(), order_.end(), x_.begin());
    apply_unorder(order_.begin(), order_.end(), y_.begin());
    apply_unorder(order_.begin(), order_.end(), z_.begin());
    
    apply_unorder(order_.begin(), order_.end(), nx_.begin());
    apply_unorder(order_.begin(), order_.end(), ny_.begin());
    apply_unorder(order_.begin(), order_.end(), nz_.begin());
    
    apply_unorder(order_.begin(), order_.end(), area_.begin());
    apply_unorder(order_.begin(), order_.end(), source_term_.begin());
    apply_unorder(order_.begin(), order_.end(), source_term_.begin() + num_);
    
    apply_unorder(order_.begin(), order_.end(), potential.begin());
    apply_unorder(order_.begin(), order_.end(), potential.begin() + num_);
}


void Elements::compute_charges(const double* __restrict potential_ptr)
{
    timers_.compute_charges.start();

    std::size_t num = num_;

    const double* __restrict nx_ptr = nx_.data();
    const double* __restrict ny_ptr = ny_.data();
    const double* __restrict nz_ptr = nz_.data();
    const double* __restrict area_ptr = area_.data();

    double* __restrict target_q_ptr    = target_charge_.data();
    double* __restrict target_q_dx_ptr = target_charge_dx_.data();
    double* __restrict target_q_dy_ptr = target_charge_dy_.data();
    double* __restrict target_q_dz_ptr = target_charge_dz_.data();

    double* __restrict source_q_ptr    = source_charge_.data();
    double* __restrict source_q_dx_ptr = source_charge_dx_.data();
    double* __restrict source_q_dy_ptr = source_charge_dy_.data();
    double* __restrict source_q_dz_ptr = source_charge_dz_.data();

#ifdef OPENACC_ENABLED
    #pragma acc parallel loop present(nx_ptr, ny_ptr, nz_ptr, area_ptr, potential_ptr, \
                target_q_ptr, target_q_dx_ptr, target_q_dy_ptr, target_q_dz_ptr, \
                source_q_ptr, source_q_dx_ptr, source_q_dy_ptr, source_q_dz_ptr)

#elif OPENMP_ENABLED
    #pragma omp parallel for
#endif
    for (std::size_t i = 0; i < num; ++i) {
        target_q_ptr   [i] = constants::ONE_OVER_4PI;
        target_q_dx_ptr[i] = constants::ONE_OVER_4PI * nx_ptr[i];
        target_q_dy_ptr[i] = constants::ONE_OVER_4PI * ny_ptr[i];
        target_q_dz_ptr[i] = constants::ONE_OVER_4PI * nz_ptr[i];
        
        source_q_ptr   [i] =             area_ptr[i] * potential_ptr[num + i];
        source_q_dx_ptr[i] = nx_ptr[i] * area_ptr[i] * potential_ptr[i];
        source_q_dy_ptr[i] = ny_ptr[i] * area_ptr[i] * potential_ptr[i];
        source_q_dz_ptr[i] = nz_ptr[i] * area_ptr[i] * potential_ptr[i];
    }

    timers_.compute_charges.stop();
}


void Elements::copyin_to_device() const
{
    timers_.copyin_to_device.start();

#ifdef OPENACC_ENABLED
    const double* x_ptr = x_.data();
    const double* y_ptr = y_.data();
    const double* z_ptr = z_.data();

    std::size_t x_num = x_.size();
    std::size_t y_num = y_.size();
    std::size_t z_num = z_.size();

    const double* nx_ptr = nx_.data();
    const double* ny_ptr = ny_.data();
    const double* nz_ptr = nz_.data();

    std::size_t nx_num = nx_.size();
    std::size_t ny_num = ny_.size();
    std::size_t nz_num = nz_.size();

    const double* area_ptr = area_.data();
    std::size_t area_num = area_.size();

    const double* source_term_ptr = source_term_.data();
    std::size_t source_term_num = source_term_.size();

    const double* tq_ptr    = target_charge_.data();
    const double* tq_dx_ptr = target_charge_dx_.data();
    const double* tq_dy_ptr = target_charge_dy_.data();
    const double* tq_dz_ptr = target_charge_dz_.data();

    std::size_t tq_num    = target_charge_.size();
    std::size_t tq_dx_num = target_charge_dx_.size();
    std::size_t tq_dy_num = target_charge_dy_.size();
    std::size_t tq_dz_num = target_charge_dz_.size();

    const double* sq_ptr    = source_charge_.data();
    const double* sq_dx_ptr = source_charge_dx_.data();
    const double* sq_dy_ptr = source_charge_dy_.data();
    const double* sq_dz_ptr = source_charge_dz_.data();

    std::size_t sq_num    = source_charge_.size();
    std::size_t sq_dx_num = source_charge_dx_.size();
    std::size_t sq_dy_num = source_charge_dy_.size();
    std::size_t sq_dz_num = source_charge_dz_.size();

    #pragma acc enter data copyin(x_ptr[0:x_num],  y_ptr[0:y_num], z_ptr[0:z_num], \
                nx_ptr[0:nx_num], ny_ptr[0:ny_num], nz_ptr[0:nz_num], area_ptr[0:area_num])
    #pragma acc enter data create(source_term_ptr[0:source_term_num], \
                tq_ptr[0:tq_num], tq_dx_ptr[0:tq_dx_num], tq_dy_ptr[0:tq_dy_num], tq_dz_ptr[0:tq_dz_num], \
                sq_ptr[0:sq_num], sq_dx_ptr[0:sq_dx_num], sq_dy_ptr[0:sq_dy_num], sq_dz_ptr[0:sq_dz_num])
#endif

    timers_.copyin_to_device.stop();
}


void Elements::update_source_term_on_host() const
{
#ifdef OPENACC_ENABLED
    const double* source_term_ptr = source_term_.data();
    std::size_t source_term_num = source_term_.size();

    #pragma acc update self(source_term_ptr[0:source_term_num])
#endif
}


void Elements::delete_from_device() const
{
    timers_.delete_from_device.start();

#ifdef OPENACC_ENABLED
    const double* x_ptr = x_.data();
    const double* y_ptr = y_.data();
    const double* z_ptr = z_.data();

    std::size_t x_num = x_.size();
    std::size_t y_num = y_.size();
    std::size_t z_num = z_.size();

    const double* nx_ptr = nx_.data();
    const double* ny_ptr = ny_.data();
    const double* nz_ptr = nz_.data();

    std::size_t nx_num = nx_.size();
    std::size_t ny_num = ny_.size();
    std::size_t nz_num = nz_.size();

    const double* area_ptr = area_.data();
    std::size_t area_num = area_.size();

    const double* source_term_ptr = source_term_.data();
    std::size_t source_term_num = source_term_.size();

    const double* tq_ptr    = target_charge_.data();
    const double* tq_dx_ptr = target_charge_dx_.data();
    const double* tq_dy_ptr = target_charge_dy_.data();
    const double* tq_dz_ptr = target_charge_dz_.data();

    std::size_t tq_num    = target_charge_.size();
    std::size_t tq_dx_num = target_charge_dx_.size();
    std::size_t tq_dy_num = target_charge_dy_.size();
    std::size_t tq_dz_num = target_charge_dz_.size();

    const double* sq_ptr    = source_charge_.data();
    const double* sq_dx_ptr = source_charge_dx_.data();
    const double* sq_dy_ptr = source_charge_dy_.data();
    const double* sq_dz_ptr = source_charge_dz_.data();

    std::size_t sq_num    = source_charge_.size();
    std::size_t sq_dx_num = source_charge_dx_.size();
    std::size_t sq_dy_num = source_charge_dy_.size();
    std::size_t sq_dz_num = source_charge_dz_.size();

    #pragma acc exit data delete(x_ptr[0:x_num],  y_ptr[0:y_num], z_ptr[0:z_num], \
                nx_ptr[0:nx_num], ny_ptr[0:ny_num], nz_ptr[0:nz_num], area_ptr[0:area_num])
    #pragma acc exit data delete(source_term_ptr[0:source_term_num], \
                tq_ptr[0:tq_num], tq_dx_ptr[0:tq_dx_num], tq_dy_ptr[0:tq_dy_num], tq_dz_ptr[0:tq_dz_num], \
                sq_ptr[0:sq_num], sq_dx_ptr[0:sq_dx_num], sq_dy_ptr[0:sq_dy_num], sq_dz_ptr[0:sq_dz_num])
#endif

    timers_.delete_from_device.stop();
}


static double triangle_area(std::array<std::array<double, 3>, 3> v)
{
    std::array<double, 3> a, b, c;

    for (int i = 0; i < 3; ++i) {
        a[i] = v[i][0] - v[i][1];
        b[i] = v[i][0] - v[i][2];
        c[i] = v[i][1] - v[i][2];
    }

    double aa = std::sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    double bb = std::sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
    double cc = std::sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);

    double ss = 0.5 * (aa + bb + cc);
    return std::sqrt(ss * (ss-aa) * (ss-bb) * (ss-cc));
}


void Timers_Elements::print() const
{
    std::cout.setf(std::ios::fixed, std::ios::floatfield);
    std::cout.precision(5);
    std::cout << "|...Elements function times (s)...." << std::endl;
    std::cout << "|   |...ctor.......................: ";
    std::cout << std::setw(12) << std::right << ctor.elapsed_time() << std::endl;
    std::cout << "|   |...compute_source_term........: ";
    std::cout << std::setw(12) << std::right << compute_source_term.elapsed_time() << std::endl;
    std::cout << "|   |...compute_charges............: ";
    std::cout << std::setw(12) << std::right << compute_charges.elapsed_time() << std::endl;
#ifdef OPENACC_ENABLED
    std::cout << "|   |...copyin_to_device...........: ";
    std::cout << std::setw(12) << std::right << copyin_to_device.elapsed_time() << std::endl;
    std::cout << "|   |...delete_from_device.........: ";
    std::cout << std::setw(12) << std::right << copyin_to_device.elapsed_time() << std::endl;
#endif
    std::cout << "|" << std::endl;
}


std::string Timers_Elements::get_durations() const
{
    std::string durations;
    durations.append(std::to_string(ctor                     .elapsed_time())).append(", ");
    durations.append(std::to_string(compute_source_term      .elapsed_time())).append(", ");
    durations.append(std::to_string(compute_charges          .elapsed_time())).append(", ");
    durations.append(std::to_string(copyin_to_device         .elapsed_time())).append(", ");
    durations.append(std::to_string(delete_from_device       .elapsed_time())).append(", ");
    
    return durations;
}


std::string Timers_Elements::get_headers() const
{
    std::string headers;
    headers.append("Elements ctor, ");
    headers.append("Elements compute_source_term, ");
    headers.append("Elements compute_charges, ");
    headers.append("Elements copyin_to_device, ");
    headers.append("Elements delete_from_device, ");
    
    return headers;
}
