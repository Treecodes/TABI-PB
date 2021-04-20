#include <algorithm>
#include <iostream>
#include <iomanip>
#include <fstream>

#include "tabipb_timers.h"
#include "coulombic_energy_compute.h"
#include "solvation_energy_compute.h"
#include "constants.h"
#include "output.h"


Output::Output(class Molecule& mol, class Elements& elem, const struct Params& params, struct Timers_Output& timers)
    : molecule_(mol), elements_(elem), params_(params), timers_(timers), potential_offset_(elements_.num())
{
    timers_.ctor.start();
    potential_.assign(2 * elements_.num(), 0.);
    timers_.ctor.stop();
}




void Output::compute_coulombic_energy()
{
    timers_.compute_coulombic_energy.start();

    double coulombic_energy = 0.;
    double epsp = params_.phys_eps_solute_;
    std::size_t num_atoms = molecule_.num();
    
    const double* __restrict mol_x_ptr = molecule_.x_ptr();
    const double* __restrict mol_y_ptr = molecule_.y_ptr();
    const double* __restrict mol_z_ptr = molecule_.z_ptr();
    const double* __restrict mol_q_ptr = molecule_.charge_ptr();

#ifdef OPENACC_ENABLED
    #pragma acc parallel loop gang present(mol_x_ptr, mol_y_ptr, mol_z_ptr, mol_q_ptr) \
                                   reduction(+:coulombic_energy)
#elif OPENMP_ENABLED
    #pragma omp parallel for reduction(+:coulombic_energy)
#endif
    for (std::size_t i = 0; i < num_atoms; ++i) {
        double xx = mol_x_ptr[i];
        double yy = mol_y_ptr[i];
        double zz = mol_z_ptr[i];
        double qq = mol_q_ptr[i];
        
#ifdef OPENACC_ENABLED
	#pragma acc loop vector reduction(+:coulombic_energy)
#endif
        for (std::size_t j = i+1; j < num_atoms; ++j) {
            double dx = xx - mol_x_ptr[j];
            double dy = yy - mol_y_ptr[j];
            double dz = zz - mol_z_ptr[j];
            coulombic_energy += qq * mol_q_ptr[j] / epsp / std::sqrt(dx*dx + dy*dy + dz*dz);
        }
    }

    coulombic_energy_ = coulombic_energy;

    timers_.compute_coulombic_energy.stop();
}




void Output::compute_coulombic_energy(const class InterpolationPoints& mol_interp_pts,  const class Tree& mol_tree,
                                      const class InteractionList& interaction_list)
{
    timers_.compute_coulombic_energy.start();

    class CoulombicEnergyCompute coulombic_energy(molecule_, mol_interp_pts, mol_tree,
                                                  interaction_list, params_.phys_eps_solute_);
                                                  
    coulombic_energy_ = coulombic_energy.compute();

    timers_.compute_coulombic_energy.stop();
}




void Output::compute_solvation_energy()
{
    timers_.compute_solvation_energy.start();

    double eps = params_.phys_eps_;
    double kappa = params_.phys_kappa_;
    double solvation_energy = 0.;
    std::size_t num_atoms = molecule_.num();
    std::size_t num_elems = elements_.num();
    
    const double* __restrict elem_x_ptr = elements_.x_ptr();
    const double* __restrict elem_y_ptr = elements_.y_ptr();
    const double* __restrict elem_z_ptr = elements_.z_ptr();
    
    const double* __restrict elem_nx_ptr = elements_.nx_ptr();
    const double* __restrict elem_ny_ptr = elements_.ny_ptr();
    const double* __restrict elem_nz_ptr = elements_.nz_ptr();
    
    const double* __restrict elem_area_ptr = elements_.area_ptr();
    
    const double* __restrict mol_x_ptr = molecule_.x_ptr();
    const double* __restrict mol_y_ptr = molecule_.y_ptr();
    const double* __restrict mol_z_ptr = molecule_.z_ptr();
    const double* __restrict mol_q_ptr = molecule_.charge_ptr();
    
    const double* __restrict potential_ptr = potential_.data();
    
#ifdef OPENACC_ENABLED
    #pragma acc enter data copyin(potential_ptr[0:potential_num])
    #pragma acc parallel loop gang present(mol_x_ptr, mol_y_ptr, mol_z_ptr, mol_q_ptr, \
                                      elem_x_ptr, elem_y_ptr, elem_z_ptr, \
                                      elem_nx_ptr, elem_ny_ptr, elem_nz_ptr, \
                                      elem_area_ptr) \
                                   reduction(+:solvation_energy)
#endif
    for (std::size_t i = 0; i < num_elems; ++i) {

#ifdef OPENACC_ENABLED
        #pragma acc loop vector reduction(+:solvation_energy)
#endif
        for (std::size_t j = 0; j < num_atoms; ++j) {
        
            double dx = elem_x_ptr[i] - mol_x_ptr[j];
            double dy = elem_y_ptr[i] - mol_y_ptr[j];
            double dz = elem_z_ptr[i] - mol_z_ptr[j];
            double r  = std::sqrt(dx*dx + dy*dy + dz*dz);

            double cos_theta   = (elem_nx_ptr[i] * dx
                                + elem_ny_ptr[i] * dy
                                + elem_nz_ptr[i] * dz) / r;

            double kappa_r     = kappa * r;
            double exp_kappa_r = std::exp(-kappa_r);

            double G0 = constants::ONE_OVER_4PI / r;
            double Gk = exp_kappa_r * G0;
            double G1 = cos_theta * G0 / r;
            double G2 = G1 * (1.0 + kappa_r) * exp_kappa_r;
        
            double L1 = G1 - eps * G2;
            double L2 = G0 - Gk;

            solvation_energy += mol_q_ptr[j] * elem_area_ptr[i]
                              * (L1 * potential_ptr[i] + L2 * potential_ptr[potential_offset_ + i]);
        }
    }
#ifdef OPENACC_ENABLED
    #pragma acc exit data delete(potential_ptr[0:potential_num])
#endif

    solvation_energy_ = solvation_energy;

    timers_.compute_solvation_energy.stop();
}




void Output::compute_solvation_energy(const class InterpolationPoints& elem_interp_pts, const class Tree& elem_tree,
                                      const class InterpolationPoints& mol_interp_pts,  const class Tree& mol_tree,
                                      const class InteractionList& interaction_list)
{
    timers_.compute_solvation_energy.start();
    const double* __restrict potential_ptr = potential_.data();
    
#ifdef OPENACC_ENABLED
    #pragma acc enter data copyin(potential_ptr[0:potential_num])
#endif
    
    class SolvationEnergyCompute solvation_energy(potential_,
                                                  elements_, elem_interp_pts, elem_tree,
                                                  molecule_, mol_interp_pts, mol_tree,
                                                  interaction_list, params_.phys_eps_, params_.phys_kappa_);
                                                  
    solvation_energy_ = solvation_energy.compute();
    
#ifdef OPENACC_ENABLED
    #pragma acc exit data delete(potential_ptr[0:potential_num])
#endif

    timers_.compute_solvation_energy.stop();
}




void Output::finalize()
{
    timers_.finalize.start();

    solvation_energy_ = constants::UNITS_PARA  * solvation_energy_;
    coulombic_energy_ = constants::UNITS_COEFF * coulombic_energy_;
    free_energy_      = solvation_energy_ + coulombic_energy_;
    
    constexpr double pot_scaling = constants::UNITS_COEFF * constants::PI * 4.;
    std::transform(std::begin(potential_), std::end(potential_),
                   std::begin(potential_), [=](double x){ return x * pot_scaling; });
                   
    auto pot_min_max = std::minmax_element(
        potential_.begin(), potential_.begin() + potential_.size() / 2);
                                           
    auto pot_normal_min_max = std::minmax_element(
        potential_.begin() + potential_.size() / 2, potential_.end());
        
    pot_min_ = *pot_min_max.first;
    pot_max_ = *pot_min_max.second;
    
    pot_normal_min_ = *pot_normal_min_max.first;
    pot_normal_max_ = *pot_normal_min_max.second;
    
    elements_.unorder(potential_);
    molecule_.unorder();
    
    timers_.finalize.stop();
}


void Output::files(const struct Timers& timers) const
{
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\n\n*** OUTPUT FOR TABI-PB RUN ***";
    std::cout << "\n\n    Solvation energy = " << solvation_energy_
                                               << " kJ/mol";
    std::cout << "\n         Free energy = "   << free_energy_
                                               << " kJ/mol";
    std::cout << "\n\nThe max and min potential and normal derivatives on vertices:";
    std::cout << "\n        Potential min: " << pot_min_ << ", "
                                     "max: " << pot_max_;
    std::cout << "\nNormal derivative min: " << pot_normal_min_ << ", "
                                     "max: " << pot_normal_max_ << "\n" << std::endl << std::endl;
    
    if (params_.output_vtk_) Output::output_VTK();
    if (params_.output_timers_) timers.print();

    if (params_.output_csv_headers_) {
        std::ofstream csv_headers("headers.csv");
        std::string headers;
        headers.append("num_atoms, ")            .append("mesh, ")
               .append("mesh_density, ")         .append("mesh_probe_radius, ")
               .append("tree_degree, ")          .append("tree_theta, ")
               .append("tree_max_per_leaf, ")    .append("precondition, ")
               .append("num_elements, ")         .append("surface_area, ")
               .append("num_iterations, ")       .append("residual, ")
               .append("solvation_energy, ")     .append("coulombic_energy, ")
               .append("free_energy, ")
               .append("potential_min, ")        .append("potential_max, ")
               .append("potential_normal_min, ") .append("potential_normal_max, ")
               .append(timers.get_headers());
        csv_headers << headers << std::endl;
        csv_headers.close();
    }

    if (params_.output_csv_) {
        std::ofstream csv_file("output.csv");
        csv_file << std::scientific << std::setprecision(12)
                 << molecule_.num()            << ", " << params_.mesh_              << ", "
                 << params_.mesh_density_      << ", " << params_.mesh_probe_radius_ << ", "
                 << params_.tree_degree_       << ", " << params_.tree_theta_        << ", "
                 << params_.tree_max_per_leaf_ << ", " << params_.precondition_      << ", "
                 << elements_.num()            << ", " << elements_.surface_area()   << ", "
                 << num_iter_                  << ", " << residual_                  << ", "
                 << solvation_energy_          << ", " << coulombic_energy_          << ", "
                 << free_energy_               << ", "
                 << pot_min_                   << ", " << pot_max_                   << ", "
                 << pot_normal_min_            << ", " << pot_normal_max_            << ", "
                 << timers.get_durations()    << std::endl;
        csv_file.close();
    }
}


void Output::output_VTK() const
{
    timers_.output_VTK.start();

    std::ofstream vtk_file("output.vtk");
    vtk_file << "# vtk DataFile Version 1.0\n";
    vtk_file << "vtk file output.vtk\n";
    vtk_file << "ASCII\n";
    vtk_file << "DATASET POLYDATA\n\n";
    
    vtk_file << "POINTS " << elements_.num() << " double\n";
    vtk_file << std::fixed << std::setprecision(6);
    for (std::size_t i = 0; i < elements_.num(); ++i)
        vtk_file << elements_.x_ptr()[i] << " " << elements_.y_ptr()[i] << " " << elements_.z_ptr()[i] << "\n";
        
    vtk_file << "POLYGONS " << elements_.num_faces() << " " << elements_.num_faces() * 4 << "\n";
    for (std::size_t i = 0; i < elements_.num_faces(); ++i)
        vtk_file << "3 " << elements_.face_x_ptr()[i]-1 << " " << elements_.face_y_ptr()[i]-1 << " "
                         << elements_.face_z_ptr()[i]-1 << "\n";

    // These are in KCAL. Multiplying by KCAL_TO_KJ would make them KJ.
    vtk_file << "\nPOINT_DATA " << elements_.num() << "\n";
    vtk_file << "SCALARS Potential double\n";
    vtk_file << "LOOKUP_TABLE default\n";
    std::copy(potential_.begin(), potential_.begin() + potential_offset_,
              std::ostream_iterator<double>(vtk_file, "\n"));
    
    // If we want induced surface charges, we can multiply NormalPotential by (1/eps + 1)
    vtk_file << "SCALARS NormalPotential double\n";
    vtk_file << "LOOKUP_TABLE default\n";
    std::copy(potential_.begin() + potential_offset_, potential_.end(),
              std::ostream_iterator<double>(vtk_file, "\n"));

    vtk_file << "\nNORMALS Normals double\n";
    for (std::size_t i = 0; i < elements_.num(); ++i)
        vtk_file << elements_.nx_ptr()[i] << " " << elements_.ny_ptr()[i] << " " << elements_.nz_ptr()[i] << "\n";
        
    vtk_file << std::endl;
    vtk_file.close();

    timers_.output_VTK.stop();
}


void Timers_Output::print() const
{
    std::cout.setf(std::ios::fixed, std::ios::floatfield);
    std::cout.precision(5);
    std::cout << "|...Output function times (s)......" << std::endl;
    std::cout << "|   |...ctor.......................: ";
    std::cout << std::setw(12) << std::right << ctor.elapsed_time() << std::endl;
    std::cout << "|   |...compute_coulombic_energy...: ";
    std::cout << std::setw(12) << std::right << compute_coulombic_energy.elapsed_time() << std::endl;
    std::cout << "|   |...compute_solvation_energy...: ";
    std::cout << std::setw(12) << std::right << compute_solvation_energy.elapsed_time() << std::endl;
    std::cout << "|   |...finalize...................: ";
    std::cout << std::setw(12) << std::right << finalize.elapsed_time() << std::endl;
    std::cout << "|       |...output_VTK.............: ";
    std::cout << std::setw(12) << std::right << output_VTK.elapsed_time() << std::endl;
    std::cout << "|" << std::endl;
}


std::string Timers_Output::get_durations() const
{
    std::string durations;
    durations.append(std::to_string(ctor                     .elapsed_time())).append(", ");
    durations.append(std::to_string(compute_coulombic_energy .elapsed_time())).append(", ");
    durations.append(std::to_string(compute_solvation_energy .elapsed_time())).append(", ");
    durations.append(std::to_string(finalize                 .elapsed_time())).append(", ");
    durations.append(std::to_string(output_VTK               .elapsed_time())).append(", ");
    return durations;
}


std::string Timers_Output::get_headers() const
{
    std::string headers;
    headers.append("Output ctor, ");
    headers.append("Output compute_coulombic_energy, ");
    headers.append("Output compute_solvation_energy, ");
    headers.append("Output finalize, ");
    headers.append("Output output_VTK, ");
    
    return headers;
}
