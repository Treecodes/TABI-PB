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
                                     
    if (bem.params_.output_csv_) {
        std::ofstream csv_file("output.csv");
        csv_file << bem.molecule_.num_atoms()      << ", " << bem.params_.mesh_ << ", "
                 << bem.params_.mesh_density_      << ", " << bem.params_.mesh_probe_radius_ << ", "
                 << bem.params_.tree_degree_       << ", " << bem.params_.tree_theta_ << ", "
                 << bem.params_.tree_max_per_leaf_ << ", " << bem.params_.precondition_ << ", "
                 << bem.particles_.num()           << ", " << bem.particles_.surface_area() << ", "
                 << bem.solvation_energy_          << ", " << bem.coulombic_energy_ << ", "
                 << bem.pot_min_                   << ", " << bem.pot_max_ << ", "
                 << bem.pot_normal_min_            << ", " << bem.pot_normal_max_ << ", "
                 << bem.num_iter_ << std::endl;
    }
    
    if (bem.params_.output_vtk_) bem.particles_.output_VTK(bem.potential_);

    return std::array<double, 3> {bem.solvation_energy_, bem.coulombic_energy_, bem.free_energy_};
}
