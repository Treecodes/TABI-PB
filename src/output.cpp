#include <iostream>

#include "constants.h"
#include "treecode.h"

/********************************************************/
void Treecode::output()
{
    particles_.unorder(potential_);
    
    auto solvation_energy = constants::UNITS_PARA  * particles_.compute_solvation_energy(potential_);
    auto coulombic_energy = constants::UNITS_COEFF * molecule_.compute_coulombic_energy();

    constexpr double pot_scaling = constants::UNITS_COEFF * constants::PI * 4.;
    std::transform(std::begin(potential_), std::end(potential_),
                   std::begin(potential_), [](double x){ return x * pot_scaling; });
                   
    auto pot_min_max = std::minmax_element(
        potential_.begin(), potential_.begin() + potential_.size() / 2);
                                           
    auto pot_normal_min_max = std::minmax_element(
        potential_.begin() + potential_.size() / 2, potential_.end());
        
    std::cout << "\n\n*** OUTPUT FOR TABI-PB RUN ***";
    std::cout << "\n\n    Solvation energy = " << solvation_energy
                                               << " kJ/mol";
    std::cout << "\n         Free energy = "   << solvation_energy + coulombic_energy
                                               << " kJ/mol";
    std::cout << "\n\nThe max and min potential and normal derivatives on vertices:";
    std::cout << "\n        Potential min: " << *pot_min_max.first << ", "
                                     "max: " << *pot_min_max.second;
    std::cout << "\nNormal derivative min: " << *pot_normal_min_max.first << ", "
                                     "max: " << *pot_normal_min_max.second << "\n" << std::endl;
                                         
    //if (params_.output_dat_) Treecode::output_DAT()
    //if (params_.output_vtk_) Treecode::output_VTK()
    //if (params_.output_csv_) Treecode::output_CSV()
}

/*
int OutputDAT(char name[256], TABIPBvars *vars)
{
    char fname[256];
    int i;

    sprintf(fname, "%s.dat", name);

    FILE *fp = fopen(fname, "w");
    fprintf(fp, "%d %d\n", vars->nspt, vars->nface);

    for (i = 0; i < vars->nspt; i++)
        fprintf(fp, "%d %f %f %f %f %f %f %f %f\n", i,
                vars->vert[0][i], vars->vert[1][i], vars->vert[2][i],
                vars->snrm[0][i], vars->snrm[1][i], vars->snrm[2][i],
                vars->vert_ptl[i], vars->vert_ptl[i + vars->nspt]);

    for (i = 0; i < vars->nface; i++)
        fprintf(fp, "%d %d %d\n", vars->face[0][i], vars->face[1][i],
                                  vars->face[2][i]);
    fclose(fp);
    
    return 0;
}


int OutputVTK(char name[256], TABIPBvars *vars)
{
    char fname[256], nspt_str[20], nface_str[20], nface4_str[20];
    int i;
    
    sprintf(nspt_str, "%d", vars->nspt);
    sprintf(nface_str, "%d", vars->nface);
    sprintf(nface4_str, "%d", vars->nface * 4);

    sprintf(fname, "%s.vtk", name);

    FILE *fp = fopen(fname, "w");

    fprintf(fp, "# vtk DataFile Version 1.0\n");
    fprintf(fp, "vtk file %s\n", fname);
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET POLYDATA\n\n");

    fprintf(fp, "POINTS %s double\n", nspt_str);
    for (i = 0; i < vars->nspt; i++) {
        fprintf(fp, "%f %f %f\n", vars->vert[0][i], vars->vert[1][i],
                                    vars->vert[2][i]);
    }

    fprintf(fp, "POLYGONS %s %s\n", nface_str, nface4_str);
    for (i = 0; i < vars->nface; i++) {
        fprintf(fp, "3 %d %d %d\n", vars->face[0][i] - 1, vars->face[1][i] - 1,
                                    vars->face[2][i] - 1);
    }

    fprintf(fp, "\nPOINT_DATA %s\n", nspt_str);
    fprintf(fp, "SCALARS PotentialVert double\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (i = 0; i < vars->nspt; i++) {
        fprintf(fp, "%f\n", KCAL_TO_KJ * vars->vert_ptl[i]);
    }

    fprintf(fp, "SCALARS NormalPotentialVert double\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (i = 0; i < vars->nspt; i++) {
        fprintf(fp, "%f\n", KCAL_TO_KJ * vars->vert_ptl[vars->nspt + i]);
    }

    //if we want induced surface charges, we can multiply vertnorm by (1/eps + 1)
    fprintf(fp, "\nNORMALS VertNorms double\n");
    for (i = 0; i < vars->nspt; i++) {
        fprintf(fp, "%f %f %f\n", vars->snrm[0][i], vars->snrm[1][i],
                                    vars->snrm[2][i]);
    }

    fprintf(fp, "\nCELL_DATA %s\n", nface_str);
    fprintf(fp, "SCALARS PotentialFace double\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (i = 0; i < vars->nface; i++) {
        fprintf(fp, "%f\n", KCAL_TO_KJ * vars->xvct[i]);
    }

    //if we want induced surface charges, we can multiply vertnorm by (1/eps + 1)
    fprintf(fp, "SCALARS NormalPotentialFace double\n");
    fprintf(fp, "LOOKUP_TABLE default\n");
    for (i = 0; i < vars->nface; i++) {
        fprintf(fp, "%f\n", KCAL_TO_KJ * vars->xvct[vars->nface + i]);
    }
        
    fclose(fp);

    return 0;
}




int OutputCSV(TABIPBparm *parm, TABIPBvars *vars, double cpu_time)
{
    char timestr[64], meshtype[20], fname[256];
    time_t t = time(NULL);
    struct tm *tm = localtime(&t);
    
    strftime(timestr, sizeof(timestr), "%c", tm);

    sprintf(meshtype, "mesh flag: %d", parm->mesh_flag);
    sprintf(fname, "tabirunsdat.csv");

    FILE *fp = fopen(fname, "a");

    fprintf(fp, "%s, %s, %d, %s, %f, "
            "%f, %d, "
            "%e, %d, %d, %e, %e, "
            "%e, %e, %e, %e, "
            "%e, %d \n",
            timestr, parm->fname, vars->natm, meshtype, parm->density,
            parm->theta, parm->order,
            vars->soleng, vars->nspt, vars->nface, 1.0/vars->nface, vars->surface_area,
            vars->max_xvct, vars->min_xvct, vars->max_der_xvct, vars->min_der_xvct,
            cpu_time, vars->gmres_iter);

    fclose(fp);
    
    return 0;
}
*/
