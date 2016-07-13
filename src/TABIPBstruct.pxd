cdef extern from "src/TABIPBstruct.h":

  ctypedef struct TABIPBparm:
    char fpath[256];
    char fname[5];
    double density;
    double probe_radius;
    double temp;
    double epsp;
    double epsw;
    double bulk_strength;
    int order;
    int maxparnode;
    double theta;
    int mesh_flag;
    int number_of_lines;

  ctypedef struct TABIPBvars:
    double soleng;
    int nspt, nface, natm;
    int **extr_v;
    int **extr_f;
    double *atmrad;
    double *atmchr;
    double *chrpos;
    double **vert;
    double **snrm;
    double *vert_ptl;
    double max_vert_ptl, min_vert_ptl;
    double max_der_vert_ptl, min_der_vert_ptl;
    double *xvct;
    double max_xvct, min_xvct;
    double max_der_xvct, min_der_xvct;
    int **face;

  int sphinx2tabipb(TABIPBparm *parm, TABIPBvars *vars);
