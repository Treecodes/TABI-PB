/**************************************************************************
* FILE NAME: readin.c                                                     *
*                                                                         *
* PURPOSE: meshes surface with msms or NanoShaper and reads in surface    *
*          data                                                           *
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
* 01/12/2018  Leighton Wilson   Created                                   *
* 08/05/2016  Leighton Wilson   Adding Windows support for msms and       *
*                               NanoShaper                                *
*                                                                         *
**************************************************************************/

#include <time.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#ifdef MPI_ENABLED
    #include <mpi.h>
#endif

#include "readin.h"
#include "utilities.h"

#include "global_params.h"
#include "array.h"

#include "TABIPBstruct.h"

/* function to read in molecule information */
int Readin(TABIPBparm *parm, TABIPBvars *vars)
{
    FILE *fp, *nsfp;
    char c;
    char fname_tp[256];

    int i, j, i1, i2, i3, j1, j2, j3, ii, jj, kk;
    int nfacenew, ichanged, ierr;
    int natm_msms;

    double den, prob_rds, a1, a2, a3, b1, b2, b3;
    int jface[3], iface[3];
    double dist_local, area_local;
    double r[3][3], xx[3], yy[3];

    int nface = -1;
    int *face_copy[3], *face[3];
    
    int rank = 0, num_procs = 1;
    
#ifdef MPI_ENABLED
    ierr = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    ierr = MPI_Comm_size(MPI_COMM_WORLD, &num_procs);
#endif

  /* Run msms */
    if (rank == 0) {
        if (parm->mesh_flag == 0) {

    #ifdef _WIN32
            sprintf(fname_tp, "msms.exe -if molecule.xyzr -prob %f -dens %f -of molecule",
                    parm->probe_radius,parm->density);
    #else
            sprintf(fname_tp, "msms -if molecule.xyzr -prob %f -dens %f -of molecule",
                    parm->probe_radius,parm->density);
    #endif

            printf("%s\n", fname_tp);
            printf("Running MSMS...\n");
            ierr = system(fname_tp);

  /* Run NanoShaper */
        } else if (parm->mesh_flag == 1 ||
                   parm->mesh_flag == 2 ||
                   parm->mesh_flag == 3) {
            nsfp = fopen("surfaceConfiguration.prm", "w");
            fprintf(nsfp, "Grid_scale = %f\n", parm->density);
            fprintf(nsfp, "Grid_perfil = %f\n", 90.0);
            fprintf(nsfp, "XYZR_FileName = molecule.xyzr\n");
            fprintf(nsfp, "Build_epsilon_maps = false\n");
            fprintf(nsfp, "Build_status_map = false\n");
            fprintf(nsfp, "Save_Mesh_MSMS_Format = true\n");
            fprintf(nsfp, "Compute_Vertex_Normals = true\n");

            if (parm->mesh_flag == 1)
                fprintf(nsfp, "Surface = ses\n");
            else if (parm->mesh_flag == 2)
                fprintf(nsfp, "Surface = skin\n");
            else if (parm->mesh_flag == 3)
                fprintf(nsfp, "Surface = blobby\n");

            fprintf(nsfp, "Smooth_Mesh = true\n");
            fprintf(nsfp, "Skin_Surface_Parameter = %f\n", 0.45);

            fprintf(nsfp, "Cavity_Detection_Filling = false\n");
            fprintf(nsfp, "Conditional_Volume_Filling_Value = %f\n", 11.4);
            fprintf(nsfp, "Keep_Water_Shaped_Cavities = false\n");
            fprintf(nsfp, "Probe_Radius = %f\n", parm->probe_radius);
            fprintf(nsfp, "Accurate_Triangulation = true\n");
            fprintf(nsfp, "Triangulation = true\n");
            fprintf(nsfp, "Check_duplicated_vertices = true\n");
            fprintf(nsfp, "Save_Status_map = false\n");
            fprintf(nsfp, "Save_PovRay = false\n");
            fprintf(nsfp, "Max_ses_patches_per_auxiliary_grid_2d_cell = %d\n", 800);

            fclose(nsfp);

            printf("Running NanoShaper...\n");

    #ifdef _WIN32
            ierr = system("NanoShaper.exe");
    #else
            ierr = system("NanoShaper");
    #endif

            rename("triangulatedSurf.face", "molecule.face");
            rename("triangulatedSurf.vert", "molecule.vert");

            remove("stderror.txt");
            remove("surfaceConfiguration.prm");
            remove("triangleAreas.txt");
            remove("exposed.xyz");
            remove("exposedIndices.txt");
        }

  /* read in vert */
        sprintf(fname_tp, "molecule.vert");

  /* open the file and read through the first two rows */
        fp = fopen(fname_tp, "r");
        for (i = 1; i <= 2; i++) {
            while ((c = getc(fp)) != '\n') {
            }
        }


        if (parm->mesh_flag == 0) {
            ierr = fscanf(fp, "%d %d %lf %lf ", &vars->nspt, &natm_msms,
                          &den, &prob_rds);
        } else if (parm->mesh_flag == 1 ||
                   parm->mesh_flag == 2 ||
                   parm->mesh_flag == 3) {
            ierr = fscanf(fp, "%d ", &vars->nspt);
        }
    }
    
#ifdef MPI_ENABLED
    MPI_Bcast(&vars->nspt, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

   /*allocate variables for vertices file*/
    vars->natm = parm->number_of_lines;    //natm is number of lines in .xyzr

    //make_matrix(vars->vert, 3, vars->nspt);
    vars->vert[0] = (double *)malloc(vars->nspt * sizeof(double));
    vars->vert[1] = (double *)malloc(vars->nspt * sizeof(double));
    vars->vert[2] = (double *)malloc(vars->nspt * sizeof(double));

    //make_matrix(vars->snrm, 3, vars->nspt);
    vars->snrm[0] = (double *)malloc(vars->nspt * sizeof(double));
    vars->snrm[1] = (double *)malloc(vars->nspt * sizeof(double));
    vars->snrm[2] = (double *)malloc(vars->nspt * sizeof(double));

    //make_vector(vars->vert_ptl, 2 * vars->nspt);
    vars->vert_ptl = (double *)malloc(2 * vars->nspt * sizeof(double));

    if (rank == 0) {
        for (i = 0; i <= vars->nspt-1; i++) {
            ierr = fscanf(fp, "%lf %lf %lf %lf %lf %lf %d %d %d",
                      &a1, &a2, &a3, &b1, &b2, &b3, &i1, &i2, &i3);

            vars->vert[0][i] = a1;
            vars->vert[1][i] = a2;
            vars->vert[2][i] = a3;
            vars->snrm[0][i] = b1;
            vars->snrm[1][i] = b2;
            vars->snrm[2][i] = b3;
        }

        fclose(fp);
        printf("Finished reading .vert file...\n");

  /* read in faces */

        sprintf(fname_tp, "molecule.face");
        fp = fopen(fname_tp, "r");
        for (i = 1; i < 3; i++) {
            while ((c = getc(fp)) != '\n') {
            }
        }

        if (parm->mesh_flag == 0) {
            ierr = fscanf(fp, "%d %d %lf %lf ", &nface, &natm_msms,
                          &den, &prob_rds);
        } else if (parm->mesh_flag == 1 ||
                   parm->mesh_flag == 2 ||
                   parm->mesh_flag == 3) {
            ierr = fscanf(fp, "%d ", &nface);
        }

        //make_matrix(face, 3, nface);
        face[0] = (int *)malloc(nface * sizeof(int));
        face[1] = (int *)malloc(nface * sizeof(int));
        face[2] = (int *)malloc(nface * sizeof(int));

        for (i = 0; i < nface; i++) {
            ierr = fscanf(fp,"%d %d %d %d %d",&j1,&j2,&j3,&i1,&i2);
            face[0][i] = j1;
            face[1][i] = j2;
            face[2][i] = j3;
        }

        fclose(fp);
        printf("Finished reading .face file...\n");

  /*
   * Delete triangles that either:
   *    + have extremely small areas, or
   *    + are too close to each other
   */
        nfacenew = nface;

        //make_matrix(face_copy, 3, nface);
        face_copy[0] = (int *)malloc(nface * sizeof(int));
        face_copy[1] = (int *)malloc(nface * sizeof(int));
        face_copy[2] = (int *)malloc(nface * sizeof(int));

        for (i = 0; i < 3; i++)
            memcpy(face_copy[i], face[i], nface * sizeof(int));

        for (i = 0; i < nface; i++) {
            for (ii = 0; ii < 3; ii++) {
                iface[ii] = face[ii][i];
                xx[ii] = 0.0;
            }

            for (ii = 0; ii < 3; ii++) {
                for (kk = 0; kk < 3; kk++) {
                    r[kk][ii] = vars->vert[kk][iface[ii]-1];
                    xx[kk] = xx[kk] + 1.0/3.0 * r[kk][ii];
                }
            }

            area_local = TriangleArea(r);

            for (j = i-10; (j >= 0 && j < i); j++) {
                for (jj = 0; jj < 3; jj++) {
                    jface[jj] = face[jj][j];
                    yy[jj] = 0.0;
                }

                for (jj = 0; jj < 3; jj++) {
                    for (kk = 0; kk < 3; kk++) {
                        r[kk][jj] = vars->vert[kk][jface[jj]-1];
                        yy[kk] = yy[kk] + 1.0/3.0 * r[kk][jj];
                    }
                }

                dist_local = 0.0;

                for (jj = 0; jj < 3; jj++) {
                    dist_local += (xx[jj]-yy[jj]) * (xx[jj]-yy[jj]);
                }

                dist_local = sqrt(dist_local);

                if (dist_local < 1e-5) {
                    goto exit;
                }
            }

            if (area_local < 1e-5) {
                goto exit;
            }

            continue;
            exit: ichanged = nface - nfacenew;

            for (ii = i - ichanged; ii < nface; ii++) {
                for (jj = 0; jj < 3; jj++) {
                    face_copy[jj][ii] = face_copy[jj][ii+1];
                }
            }

            nfacenew = nfacenew - 1;
        }
        printf("\n%d faces have been deleted.\n", nface - nfacenew);
        nface = nfacenew;
        vars->nface = nface;
    }
    
#ifdef MPI_ENABLED
    MPI_Bcast(&vars->nface, 1, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    //make_matrix(vars->face, 3, vars->nface);
    vars->face[0] = (int *)malloc(vars->nface * sizeof(int));
    vars->face[1] = (int *)malloc(vars->nface * sizeof(int));
    vars->face[2] = (int *)malloc(vars->nface * sizeof(int));

    //make_vector(vars->xvct, 2 * vars->nface);
    vars->xvct = (double *)malloc(vars->nface * 2 * sizeof(double));

    if (rank == 0) {
        for (i = 0; i < vars->nface; i++) {
            for (j = 0; j < 3; j++) {
                vars->face[j][i] = face_copy[j][i];
            }
        }

        free(face[0]);
        free(face[1]);
        free(face[2]);
        free(face_copy[0]);
        free(face_copy[1]);
        free(face_copy[2]);
    
        remove("molecule.xyzr");
        remove("molecule.vert");
        remove("molecule.face");
    }
  
#ifdef MPI_ENABLED
    MPI_Bcast(vars->vert[0], vars->nspt, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(vars->vert[1], vars->nspt, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(vars->vert[2], vars->nspt, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(vars->snrm[0], vars->nspt, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(vars->snrm[1], vars->nspt, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(vars->snrm[2], vars->nspt, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(vars->face[0], vars->nface, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(vars->face[1], vars->nface, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(vars->face[2], vars->nface, MPI_INT, 0, MPI_COMM_WORLD);
#endif

    return 0;
}
