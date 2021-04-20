#ifndef H_TABIPB_MOLECULE_STRUCT_H
#define H_TABIPB_MOLECULE_STRUCT_H

#include <vector>
#include <string>
#include <fstream>
#include <cstddef>

#include "timer.h"
#include "particles.h"

#ifdef TABIPB_APBS
    #include "generic/valist.h"
#endif

struct Timers_Molecule;

class Molecule : public Particles
{
private:
    struct Timers_Molecule& timers_;
    
    double coulombic_energy_;
    std::vector<double> charge_;
    std::vector<double> radius_;


public:
    Molecule(struct Params&, struct Timers_Molecule&);
    ~Molecule() = default;
    
#ifdef TABIPB_APBS
    Molecule(Valist*, struct Params&, struct Timers_Molecule&);
#endif
    
    void build_xyzr_file() const;
    
    const double* charge_ptr() const { return charge_.data(); };
    const double* radius_ptr() const { return radius_.data(); };
    
    void reorder() override;
    void unorder() override;
    
    void copyin_to_device() const override;
    void delete_from_device() const override;
};


struct Timers_Molecule
{
    Timer ctor;
    Timer build_xyzr_file;
    Timer copyin_to_device;
    Timer delete_from_device;
    
    void print() const;
    std::string get_durations() const;
    std::string get_headers() const;

    Timers_Molecule() = default;
    ~Timers_Molecule() = default;
};

#endif /* H_MOLECULE_STRUCT_H */
