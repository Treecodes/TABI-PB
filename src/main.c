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

#ifdef MPI_ENABLED
    #include <mpi.h>
#endif

#include "tabipb.h"
#include "print_output.h"

#include "array.h"
#include "TABIPBstruct.h"

int main(int argc, char **argv)
{
  /* main reads the input file, writes xyzr file for msms and sets up position,
     radius and charges */

    FILE *fp, *wfp;
    char c[256], fname_tp[256], list_mol[10][256];
    char c1[120], c2[120], c3[120], c4[10], c5[10];
    double a1, a2, a3, b1, b2;
    int ierr, i, j, num_mol = 0;
    
    int rank = 0, num_procs = 1;
    
#ifdef MPI_ENABLED
    MPI_Status status;
#endif

  /* timing functions for *nix systems */
#ifndef _WIN32                                                                     
    extern void timer_start();
    extern void timer_end();
#endif

#ifndef _WIN32                                                                     
    timer_start("TOTAL_TIME");
#endif

#ifdef MPI_ENABLED
    ierr = MPI_Init(&argc, &argv);
    
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
#endif

    if (rank == 0) {
        if (argc < 2) {
            printf("No input file specified. Exiting.\n");
            return 1;
        }
    }
    
    TABIPBparm *main_parm = malloc(sizeof *main_parm);
    TABIPBvars *main_vars;

/********************************************************/
    if (rank == 0) {
    
        fp = fopen(argv[1], "r");
    
        main_parm->output_datafile = 0;
    
        while (fgets(c, 256, fp) != NULL) {
            sscanf(c, "%s %s %s", c1, c2, c3);
        
            if (strncmp(c1, "mol", 3) == 0) {
                num_mol++;
                strcpy(list_mol[num_mol-1], c2);
            
            } else if (strcmp(c1, "sdens") == 0) {
                main_parm->density = atof(c2);
            
            } else if (strcmp(c1, "srad") == 0) {
                main_parm->probe_radius = atof(c2);
            
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
                
            } else if (strcmp(c1, "outdata") == 0) {
                if (strcmp(c2, "dat") == 0 || strcmp(c2, "DAT") == 0) {
                    main_parm->output_datafile = 1;
                
                } else if (strcmp(c2, "vtk") == 0 || strcmp(c2, "VTK") == 0) {
                    main_parm->output_datafile = 2;
                }
            }
        }
    
        fclose(fp);
    }
    
    //COPY PARM TO OTHER PROCESSORS
    
#ifdef MPI_ENABLED
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
    MPI_Bcast(&num_mol, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

/********************************************************/

    
    for (j = 0; j < num_mol; j++) {
    
        main_vars = malloc(sizeof(TABIPBvars));
        strcpy(main_parm->fname, list_mol[j]);
    
        if (rank == 0) {
            printf("\n\n*** BEGINNING RUN %d: molecule file %s ***\n",
                   j+1, list_mol[j]);
        
            fp = fopen(main_parm->fname, "r");
        
            sprintf(fname_tp, "molecule.xyzr");
            wfp = fopen(fname_tp, "w");

            main_parm->number_of_lines = 0;
            while (fscanf(fp, "%s %s %s %s %s %lf %lf %lf %lf %lf",
                   c1, c2, c3, c4, c5, &a1, &a2, &a3, &b1, &b2) != EOF) {
                if (strncmp(c1, "ATOM", 4) == 0) {
                    fprintf(wfp, "%f %f %f %f\n", a1, a2, a3, b2);
                    main_parm->number_of_lines++;
                }
            }

            fclose(wfp);
        }
        
#ifdef MPI_ENABLED
        MPI_Bcast(&main_parm->number_of_lines, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

        make_vector(main_vars->chrpos, 3 * main_parm->number_of_lines);
        make_vector(main_vars->atmchr, main_parm->number_of_lines);
        make_vector(main_vars->atmrad, main_parm->number_of_lines);

        if (rank == 0) {
            rewind(fp);
            i = 0;
    
            while (fscanf(fp, "%s %s %s %s %s %lf %lf %lf %lf %lf",
                   c1, c2, c3, c4, c5, &a1, &a2, &a3, &b1, &b2) != EOF) {
                if (strncmp(c1, "ATOM", 4) == 0) {
                    main_vars->chrpos[3*i] = a1;
                    main_vars->chrpos[3*i + 1] = a2;
                    main_vars->chrpos[3*i + 2] = a3;
                    main_vars->atmchr[i] = b1;
                    main_vars->atmrad[i] = b2;
                    i++;
                }
            }

            fclose(fp);
        }

#ifdef MPI_ENABLED
        MPI_Bcast(main_vars->chrpos, 3 * main_parm->number_of_lines,
                  MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(main_vars->atmchr, main_parm->number_of_lines,
                  MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Bcast(main_vars->atmrad, main_parm->number_of_lines,
                  MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

        ierr = TABIPB(main_parm, main_vars);
        
        if (rank == 0) {
            ierr = OutputPrint(main_vars);
            
            if (main_parm->output_datafile == 1) {
                ierr = OutputDAT(main_parm, main_vars);
            } else if (main_parm->output_datafile == 2) {
                ierr = OutputVTK(main_parm, main_vars);
            }
        }

        free_vector(main_vars->atmchr);
        free_vector(main_vars->chrpos);
        free_vector(main_vars->atmrad);
        free_vector(main_vars->vert_ptl); // allocate in output_potential()
        free_vector(main_vars->xvct);
        free_matrix(main_vars->vert);
        free_matrix(main_vars->snrm);
        free_matrix(main_vars->face);
        free(main_vars);
    }
    
    free(main_parm);

#ifndef _WIN32
    timer_end();
#endif

#ifdef MPI_ENABLED
    ierr = MPI_Type_free(&mpi_tabipbparm_type);
    ierr = MPI_Finalize();
#endif
    
    return 0;
}
