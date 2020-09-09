#include <iostream>

#include "constants.h"
#include "treecode.h"

/********************************************************/


/*


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
