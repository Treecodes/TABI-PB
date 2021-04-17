#include <iostream>
#include <iomanip>
#include <fstream>
#include <array>

#include "boundary_element.h"
#include "tabipb_timers.h"

std::array<double, 3> Output(const BoundaryElement& bem, const Timers& timers)
{
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "\n\n*** OUTPUT FOR TABI-PB RUN ***";
    std::cout << "\n\n    Solvation energy = " << bem.solvation_energy_
                                               << " kJ/mol";
    std::cout << "\n         Free energy = "   << bem.free_energy_
                                               << " kJ/mol";
    std::cout << "\n\nThe max and min potential and normal derivatives on vertices:";
    std::cout << "\n        Potential min: " << bem.pot_min_ << ", "
                                     "max: " << bem.pot_max_;
    std::cout << "\nNormal derivative min: " << bem.pot_normal_min_ << ", "
                                     "max: " << bem.pot_normal_max_ << "\n" << std::endl << std::endl;

    if (bem.params_.output_timers_) timers.print();

    if (bem.params_.output_csv_headers_) {
        std::ofstream csv_headers("headers.csv");
        std::string headers;
        headers.append("num_atoms, ")            .append("mesh, ")
               .append("mesh_density, ")         .append("mesh_probe_radius, ")
               .append("tree_degree, ")          .append("tree_theta, ")
               .append("tree_max_per_leaf, ")    .append("precondition, ")
               .append("num_elements, ")        .append("surface_area, ")
               .append("num_iterations, ")       .append("residual, ")
               .append("solvation_energy, ")     .append("coulombic_energy, ")
               .append("free_energy, ")
               .append("potential_min, ")        .append("potential_max, ")
               .append("potential_normal_min, ") .append("potential_normal_max, ")
               .append(timers.get_headers());
        csv_headers << headers << std::endl;
        csv_headers.close();
    }

    if (bem.params_.output_csv_) {
        std::ofstream csv_file("output.csv");
        csv_file << std::scientific << std::setprecision(12)
                 << bem.molecule_.num()            << ", " << bem.params_.mesh_              << ", "
                 << bem.params_.mesh_density_      << ", " << bem.params_.mesh_probe_radius_ << ", "
                 << bem.params_.tree_degree_       << ", " << bem.params_.tree_theta_        << ", "
                 << bem.params_.tree_max_per_leaf_ << ", " << bem.params_.precondition_      << ", "
                 << bem.elements_.num()            << ", " << bem.elements_.surface_area()   << ", "
                 << bem.num_iter_                  << ", " << bem.residual_                  << ", "
                 << bem.solvation_energy_          << ", " << bem.coulombic_energy_          << ", "
                 << bem.free_energy_               << ", "
                 << bem.pot_min_                   << ", " << bem.pot_max_                   << ", "
                 << bem.pot_normal_min_            << ", " << bem.pot_normal_max_            << ", "
                 << timers.get_durations()         << std::endl;
        csv_file.close();
    }

    if (bem.params_.output_vtk_) bem.elements_.output_VTK(bem.potential_);

    return std::array<double, 3> {bem.solvation_energy_, bem.coulombic_energy_, bem.free_energy_};
}
