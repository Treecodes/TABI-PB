#ifndef H_TABIPB_ELEMENTS_STRUCT_H
#define H_TABIPB_ELEMENTS_STRUCT_H

#include <vector>
#include <cstdlib>

#include "timer.h"
#include "source_term_compute.h"
#include "molecule.h"
#include "particles.h"

struct Timers_Elements;

class Elements : public Particles
{
private:
    const class Molecule& molecule_;
    struct Timers_Elements& timers_;
    
    std::size_t num_faces_;
    double surface_area_;
    
    std::vector<std::size_t> face_x_;
    std::vector<std::size_t> face_y_;
    std::vector<std::size_t> face_z_;
    
    std::vector<double> nx_;
    std::vector<double> ny_;
    std::vector<double> nz_;
    
    std::vector<double> area_;
    std::vector<double> source_term_;
    
    std::vector<double> target_charge_;
    std::vector<double> target_charge_dx_;
    std::vector<double> target_charge_dy_;
    std::vector<double> target_charge_dz_;
    
    std::vector<double> source_charge_;
    std::vector<double> source_charge_dx_;
    std::vector<double> source_charge_dy_;
    std::vector<double> source_charge_dz_;
    
    void generate_elements(Params::Mesh, double, double);
    void update_source_term_on_host() const;
    
    
public:
    Elements(const class Molecule&, const struct Params&, struct Timers_Elements&);
    ~Elements() = default;
        
    std::size_t num_faces() const { return num_faces_; };
    double surface_area() const { return surface_area_; };
    
    const std::size_t* face_x_ptr() const { return face_x_.data(); };
    const std::size_t* face_y_ptr() const { return face_y_.data(); };
    const std::size_t* face_z_ptr() const { return face_z_.data(); };
    
    const double* nx_ptr() const { return nx_.data(); };
    const double* ny_ptr() const { return ny_.data(); };
    const double* nz_ptr() const { return nz_.data(); };
    
    const double* area_ptr() const { return area_.data(); };
    const double* source_term_ptr() const { return source_term_.data(); };
    
    const double* target_charge_ptr()    const { return target_charge_.data(); };
    const double* target_charge_dx_ptr() const { return target_charge_dx_.data(); };
    const double* target_charge_dy_ptr() const { return target_charge_dy_.data(); };
    const double* target_charge_dz_ptr() const { return target_charge_dz_.data(); };
    
    const double* source_charge_ptr()    const { return source_charge_.data(); };
    const double* source_charge_dx_ptr() const { return source_charge_dx_.data(); };
    const double* source_charge_dy_ptr() const { return source_charge_dy_.data(); };
    const double* source_charge_dz_ptr() const { return source_charge_dz_.data(); };
    
    void reorder() override;
    void unorder() override;
    void unorder(std::vector<double>& potential);
    
    void compute_source_term();
    void compute_source_term(const class InterpolationPoints& elem_interp_pts, const class Tree& elem_tree,
                             const class Molecule& molecule, const class InterpolationPoints& mol_interp_pts,
                             const class Tree& mol_tree, const class InteractionList& interaction_list);
    
    void compute_charges(const double* potential);
    
    void copyin_to_device() const override;
    void delete_from_device() const override;
};


struct Timers_Elements
{
    Timer ctor;
    Timer compute_source_term;
    Timer compute_charges;
    Timer copyin_to_device;
    Timer delete_from_device;
    Timer output_VTK;
    
    void print() const;
    std::string get_durations() const;
    std::string get_headers() const;

    Timers_Elements() = default;
    ~Timers_Elements() = default;
};

#endif /* H_TABIPB_ELEMENTS_STRUCT_H */
