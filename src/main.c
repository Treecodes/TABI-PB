/**************************************************************************
* FILE NAME: main.c                                                       *
*                                                                         *
* PURPOSE: calls primary tabipb routine and prints output when running    *
*          as standalone                                                  *
*                                                                         *
* AUTHORS: Leighton Wilson, University of Michigan, Ann Arbor, MI         *
*          Jiahui Chen, Southern Methodist University, Dallas, TX         *
*                                                                         *
* BASED ON PACKAGE ORIGINALLY WRITTEN IN FORTRAN BY:                      *
*          Weihua Geng, Southern Methodist University, Dallas, TX         *
*          Robery Krasny, University of Michigan, Ann Arbor, MI           *
*                                                                         *
* DEVELOPMENT HISTORY:                                                    *
*                                                                         *
* Date        Author            Description Of Change                     *
* ----        ------            ---------------------                     *
* 02/10/2018  Leighton Wilson   Adding support for multiple PQRs          *
* 01/14/2018  Leighton Wilson   Fixing read in of PQRs                    *
* 07/14/2016  Jiahui Chen       Added Sphinx support                      *
* 06/30/2016  Jiahui Chen       Rebuilt wrapper architecture              *
* 06/23/2016  Leighton Wilson   Added NanoShaper support                  *
*                                                                         *
**************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#ifdef MPI_ENABLED
    #include <mpi.h>
#endif

#include "readin.h"
#include "tabipb.h"
#include "print_output.h"

#include "array.h"
#include "TABIPBstruct.h"

static int s_ReadInputFile(char **argv, TABIPBparm *main_parm,
                           int *num_mol, char list_mol[2][3][256]);
static int s_BroadcastParms(TABIPBparm *parm, int *num_mol,
                            char list_mol[2][3][256]);
static int s_SetupXYZRFile(TABIPBparm *main_parm);
static int s_SetupAtomicVars(TABIPBparm *parm, TABIPBvars *vars);

static int s_InitializeRun(int *argc, int *rank, int *num_procs,
                           TABIPBparm **parm, TABIPBvars ***vars);
static int s_FinalizeRun(TABIPBparm **parm, TABIPBvars ***vars, int num_mol);
static int s_FreeTABIPBvars(TABIPBvars **vars);


int main(int argc, char **argv)
{
  /* main reads the input file, writes xyzr file for msms and sets up position,
     radius and charges */

    char name[256], list_mol[2][3][256];
    int ierr, j, num_mol;
    double cpu_time;
    int rank, num_procs;
    
    TABIPBparm *main_parm;
    TABIPBvars **main_vars;
    
#ifdef MPI_ENABLED
    ierr = MPI_Init(&argc, &argv);
#endif

/********************************************************/

    //Initialize variables
    ierr = s_InitializeRun(&argc, &rank, &num_procs, &main_parm, &main_vars);
    
    //Read input file and broadcast to all processes
    ierr = s_ReadInputFile(argv, main_parm, &num_mol, list_mol);
    ierr = s_BroadcastParms(main_parm, &num_mol, list_mol);

/********************************************************/
    
    for (j = 0; j < num_mol; j++) {
    
        main_vars[j] = malloc(sizeof(TABIPBvars));
        
        strcpy(main_parm->fname, list_mol[j][0]);
        main_parm->mesh_flag = atoi(list_mol[j][1]);
        main_parm->density = atof(list_mol[j][2]);
    
        //Create XYZR file from PQR file and set number of lines in main_parm
        //Setup atomic information in main_vars from PQR file
        ierr = s_SetupXYZRFile(main_parm);
        ierr = s_SetupAtomicVars(main_parm, main_vars[j]);

        clock_t cpu_begin = clock();

        /* generate surface meshes from .xyzr and save in main_vars */
        /* run TABIPB */
        ierr = Readin(main_parm, main_vars[j]);
        ierr = TABIPB(main_parm, main_vars[j]);

        clock_t cpu_end = clock();
        cpu_time = (double)(cpu_end - cpu_begin) / CLOCKS_PER_SEC;

        if (rank == 0 && main_parm->output_datafile == 3) {
            ierr = OutputCSV(main_parm, main_vars[j], cpu_time);
        }
    }
    
    if (rank == 0) {
        for (j = 0; j < num_mol; j++) {
    
            sprintf(name, "tabipb_run_%d", j+1);
        
            ierr = OutputPrint(name, main_vars[j]);
    
            if (main_parm->output_datafile == 1) {
                ierr = OutputDAT(name, main_vars[j]);
            
            } else if (main_parm->output_datafile == 2) {
                ierr = OutputVTK(name, main_vars[j]);
            }
        }
    }

    ierr = s_FinalizeRun(&main_parm, &main_vars, num_mol);
    
#ifdef MPI_ENABLED
    ierr = MPI_Finalize();
#endif

    return 0;
}



static int s_ReadInputFile(char **argv, TABIPBparm *main_parm,
                           int *num_mol, char list_mol[2][3][256])
{
    FILE *fp;
    char c[256], c1[120], c2[120], c3[120];
    int ierr, rank = 0;

#ifdef MPI_ENABLED
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif
    
    *num_mol = 0;
    main_parm->output_datafile = 0;
    
    if (rank == 0) {
        fp = fopen(argv[1], "r");

        while (fgets(c, 256, fp) != NULL) {
            sscanf(c, "%s %s %s", c1, c2, c3);
        
            if (strncmp(c1, "mol", 3) == 0) {
                (*num_mol)++;
                if ((*num_mol) > 2) {
                    printf("Error: no more than two molecules supported\n");
                    exit(1);
                }
                strcpy(list_mol[(*num_mol)-1][0], c2);
             
                if (fgets(c, 256, fp) != NULL) {
                    sscanf(c, "%s %s %s", c1, c2, c3);
                
                    if (strcmp(c1, "mesh") == 0) {
                        if (strcmp(c2, "MSMS") == 0 || strcmp(c2, "msms") == 0) {
                            //main_parm->mesh_flag = 0;
                            strcpy(list_mol[(*num_mol)-1][1], "0");
                
                        } else if (strcmp(c2, "NanoShaper") == 0 ||
                                   strcmp(c2, "nanoshaper") == 0 ||
                                   strcmp(c2, "NANOSHAPER") == 0) {
                
                            if (strcmp(c3, "SES") == 0 ||
                                strcmp(c3, "ses") == 0 ||
                                strcmp(c3, "Ses") == 0) {
                                //main_parm->mesh_flag = 1;
                                strcpy(list_mol[(*num_mol)-1][1], "1");
                    
                            } else if (strcmp(c3, "SKIN") == 0 ||
                                       strcmp(c3, "skin") == 0 ||
                                       strcmp(c3, "Skin") == 0) {
                                //main_parm->mesh_flag = 2;
                                strcpy(list_mol[(*num_mol)-1][1], "2");
                            } else {
                                //main_parm->mesh_flag = 1;
                                strcpy(list_mol[(*num_mol)-1][1], "1");
                            }
                        }
                    } else {
                        printf("Error: must specify mesh type\n");
                        exit(1);
                    }
                } else {
                    printf("Error: incomplete input file\n");
                    exit(1);
                }
            
                if (fgets(c, 256, fp) != NULL) {
                    sscanf(c, "%s %s %s", c1, c2, c3);
                    
                    if (strcmp(c1, "sdens") == 0) {
                        //main_parm->density = atof(c2);
                        strcpy(list_mol[(*num_mol)-1][2], c2);
                    } else {
                        printf("Error: must specify mesh density\n");
                        exit(1);
                    }
                } else {
                    printf("Error: incomplete input file\n");
                    exit(1);
                }
            
            } else if (strcmp(c1, "pdie") == 0) {
                main_parm->epsp = atof(c2);
            
            } else if (strcmp(c1, "sdie") == 0) {
                main_parm->epsw = atof(c2);
            
            } else if (strcmp(c1, "bulk") == 0) {
                main_parm->bulk_strength = atof(c2);
            
            } else if (strcmp(c1, "temp") == 0) {
                main_parm->temp = atof(c2);
            
            } else if (strcmp(c1, "tree_order") == 0) {
                main_parm->order = atoi(c2);
            
            } else if (strcmp(c1, "tree_n0") == 0) {
                main_parm->maxparnode = atoi(c2);
            
            } else if (strcmp(c1, "tree_mac") == 0) {
                main_parm->theta = atof(c2);
            
            } else if (strcmp(c1, "mesh") == 0) {
                if (strcmp(c2, "MSMS") == 0 || strcmp(c2, "msms") == 0) {
                    main_parm->mesh_flag = 0;
                
                } else if (strcmp(c2, "NanoShaper") == 0 ||
                           strcmp(c2, "nanoshaper") == 0 ||
                           strcmp(c2, "NANOSHAPER") == 0) {
                
                    if (strcmp(c3, "SES") == 0 ||
                        strcmp(c3, "ses") == 0 ||
                        strcmp(c3, "Ses") == 0) {
                        main_parm->mesh_flag = 1;
                    
                    } else if (strcmp(c3, "SKIN") == 0 ||
                               strcmp(c3, "skin") == 0 ||
                               strcmp(c3, "Skin") == 0) {
                        main_parm->mesh_flag = 2;
                    
                    } else {
                        main_parm->mesh_flag = 1;
                    }
                }
                
            } else if (strcmp(c1, "sdens") == 0) {
                main_parm->density = atof(c2);
            
            } else if (strcmp(c1, "srad") == 0) {
                main_parm->probe_radius = atof(c2);
                
            } else if (strcmp(c1, "outdata") == 0) {
                if (strcmp(c2, "dat") == 0 || strcmp(c2, "DAT") == 0) {
                    main_parm->output_datafile = 1;
                
                } else if (strcmp(c2, "vtk") == 0 || strcmp(c2, "VTK") == 0) {
                    main_parm->output_datafile = 2;

                } else if (strcmp(c2, "csv") == 0 || strcmp(c2, "CSV") == 0) {
                    main_parm->output_datafile = 3;
                }
            }
        }

        fclose(fp);
    }

    return 0;
}



static int s_FreeTABIPBvars(TABIPBvars **vars)
{
    //allocated in SetupAtomicVars
    free_vector((*vars)->atmchr);
    free_vector((*vars)->chrpos);
    free_vector((*vars)->atmrad);
    
    //allocated in Readin
    free_vector((*vars)->vert_ptl);
    free_vector((*vars)->xvct);
    free_matrix((*vars)->vert);
    free_matrix((*vars)->snrm);
    free_matrix((*vars)->face);
    
    free(*vars);

    return 0;
}



static int s_BroadcastParms(TABIPBparm *main_parm, int *num_mol,
                            char list_mol[2][3][256])
{
#ifdef MPI_ENABLED
    int ierr;
    int nitems = 5;
    int blocklengths[5] = {512, 9, 2, 1, 3};
    MPI_Datatype types[5] = {MPI_CHAR, MPI_DOUBLE, MPI_INT,
                             MPI_DOUBLE, MPI_INT};
    MPI_Datatype mpi_tabipbparm_type;
    MPI_Aint offsets[5];

    offsets[0] = offsetof(TABIPBparm, fpath);
    offsets[1] = offsetof(TABIPBparm, density);
    offsets[2] = offsetof(TABIPBparm, order);
    offsets[3] = offsetof(TABIPBparm, theta);
    offsets[4] = offsetof(TABIPBparm, mesh_flag);

    MPI_Type_create_struct(nitems, blocklengths, offsets, types,
                           &mpi_tabipbparm_type);
    MPI_Type_commit(&mpi_tabipbparm_type);

    MPI_Bcast(main_parm, 1, mpi_tabipbparm_type, 0, MPI_COMM_WORLD);
    MPI_Bcast(list_mol, 10*256, MPI_CHAR, 0, MPI_COMM_WORLD);
    MPI_Bcast(num_mol, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    ierr = MPI_Type_free(&mpi_tabipbparm_type);
#endif
    
    return 0;
}



static int s_SetupXYZRFile(TABIPBparm *parm)
{
    FILE *fp, *wfp;
    char c[256], c1[120], c2[120], c3[120], c4[10], c5[10];
    double a1, a2, a3, b1, b2;
    int ierr, rank = 0;

#ifdef MPI_ENABLED
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    if (rank == 0) {

        printf("\n\n*** BEGINNING RUN FOR: pqr file %s ***\n", parm->fname);
    
        fp = fopen(parm->fname, "r");
        wfp = fopen("molecule.xyzr", "w");

        parm->number_of_lines = 0;
        
        while (fgets(c, sizeof(c), fp)) {
            sscanf(c, "%s %s %s %s %s %lf %lf %lf %lf %lf",
                   c1, c2, c3, c4, c5, &a1, &a2, &a3, &b1, &b2);

            if (strncmp(c1, "ATOM", 4) == 0) {
                fprintf(wfp, "%f %f %f %f\n", a1, a2, a3, b2);
                parm->number_of_lines++;
            }
        }
    
        fclose(fp);
        fclose(wfp);
    }
    
#ifdef MPI_ENABLED
    MPI_Bcast(&parm->number_of_lines, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif
    
    return 0;
}



static int s_InitializeRun(int *argc, int *rank, int *num_procs,
                           TABIPBparm **parm, TABIPBvars ***vars)
{
    int ierr;
    
    *parm = malloc(sizeof **parm);
    *vars = malloc(3 * sizeof **vars);

#ifdef MPI_ENABLED
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, rank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, num_procs);
#endif

    if (rank == 0) {
        if (*argc < 2) {
            printf("No input file specified. Exiting.\n");
            exit(1);
        }
    }

    return 0;
}



static int s_FinalizeRun(TABIPBparm **parm, TABIPBvars ***vars, int num_mol)
{
    int j, ierr;
    
    for (j = 0; j < num_mol; j++) {
        ierr = s_FreeTABIPBvars(&((*vars)[j]));
    }
    
    free(*parm);
    free(*vars);
    
    return 0;
}



static int s_SetupAtomicVars(TABIPBparm *parm, TABIPBvars *vars)
{

    FILE *fp;
    char c[256], c1[120], c2[120], c3[120], c4[10], c5[10];
    double a1, a2, a3, b1, b2;
    int i, ierr, rank = 0;

#ifdef MPI_ENABLED
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    vars->natm = parm->number_of_lines;
    make_vector(vars->chrpos, 3 * vars->natm);
    make_vector(vars->atmchr, vars->natm);
    make_vector(vars->atmrad, vars->natm);

    if (rank == 0) {
        fp = fopen(parm->fname, "r");
        i = 0;
    
        while (fgets(c, sizeof(c), fp)) {
            sscanf(c, "%s %s %s %s %s %lf %lf %lf %lf %lf",
                    c1, c2, c3, c4, c5, &a1, &a2, &a3, &b1, &b2);

            if (strncmp(c1, "ATOM", 4) == 0) {
                vars->chrpos[3*i] = a1;
                vars->chrpos[3*i + 1] = a2;
                vars->chrpos[3*i + 2] = a3;
                vars->atmchr[i] = b1;
                vars->atmrad[i] = b2;
                i++;
            }
        }
            
        fclose(fp);
    }

#ifdef MPI_ENABLED
    MPI_Bcast(vars->chrpos, 3 * vars->natm, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(vars->atmchr, vars->natm, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(vars->atmrad, vars->natm, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

    return 0;
}
