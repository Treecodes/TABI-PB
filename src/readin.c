/*
 * C routine to mesh surfaces and read in surface data for tabipb
 *
 * C version authored by:
 * Jiahui Chen, Southern Methodist University, Dallas, TX
 *
 * Additional modifications and updates by:
 * Leighton Wilson, University of Michigan, Ann Arbor, MI
 *
 * Based on package originally written in FORTRAN by:
 * Weihua Geng, Southern Methodist University, Dallas, TX
 * Robery Krasny, University of Michigan, Ann Arbor, MI
 *
 * Last changed at 8/05/2016: Adding Windows support for msms and NanoShaper
 */

#include <time.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "global_params.h"
#include "array.h"

#include "TABIPBstruct.h"
#include "particle_struct.h"

/* function computing the area of a triangle given vertices coodinates */
double TriangleArea(double v[3][3])
{
    int i;
    double a[3], b[3], c[3], aa, bb, cc, ss, area;

    for (i = 0; i <= 2; i++) {
        a[i] = v[i][0] - v[i][1];
        b[i] = v[i][0] - v[i][2];
        c[i] = v[i][1] - v[i][2];
    }

    aa = sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
    bb = sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
    cc = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);

    ss = 0.5 * (aa + bb + cc);
    area = sqrt(ss * (ss-aa) * (ss-bb) * (ss-cc));

    return(area);
}

/* function read in molecule information */
int Readin(TABIPBparm *parm, TABIPBvars *vars, ParticleStruct *particles)
{
    FILE *fp, *wfp, *nsfp;
    char c, c1[10], c2[10], c3[10], c4[10], c5[10];
    char fname_tp[256];

    int i, j, k, i1, i2, i3, j1, j2, j3, ii, jj, kk;
    int nfacenew, nface, ichanged, ierr;
    int natm_msms;

    int namelength = 4;

    double den, prob_rds, a1, a2, a3, b1, b2, b3, v0_norm;
    double r0[3], v0[3], v[3][3], r[3][3];
    int idx[3], jface[3], iface[3];
    double rs, rd[3], pot = 0.0, sum = 0.0, pot_temp = 0.0;
    double dist_local, area_local;
    double cos_theta, G0, tp1, G1, r_s[3];
    double xx[3], yy[3];

    int **face_copy, **face;

    vars->natm = parm->number_of_lines;    //natm is number of lines in .xyzr

  /* Run msms */
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

        //sprintf(fname_tp,"mv triangulatedSurf.face molecule.face\n");
        //ierr = system(fname_tp);
        //sprintf(fname_tp,"mv triangulatedSurf.vert molecule.vert\n");
        //ierr = system(fname_tp);
        //ierr = system("rm -f stderror.txt");
        //ierr = system("rm -f surfaceConfiguration.prm");
        //ierr = system("rm -f triangleAreas.txt");
        //ierr = system("rm -f exposed.xyz");
        //ierr = system("rm -f exposedIndices.txt");
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
        ierr = fscanf(fp,"%d %d %lf %lf ",&nspt,&natm_msms,&den,&prob_rds);
    } else if (parm->mesh_flag == 1 ||
               parm->mesh_flag == 2 ||
               parm->mesh_flag == 3) {
        ierr = fscanf(fp,"%d ",&nspt);
    }


  /*allocate variables for vertices file*/

    make_matrix(vars->vert, 3, nspt);
    make_matrix(vars->snrm, 3, nspt);

    for (i = 0; i <= nspt-1; i++) {
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
        ierr=fscanf(fp,"%d %d %lf %lf ",&nface,&natm_msms,&den,&prob_rds);
        //printf("nface=%d, natm=%d, den=%lf, prob=%lf\n", nface,natm,den,prob_rds);
    } else if (parm->mesh_flag == 1 ||
               parm->mesh_flag == 2 ||
               parm->mesh_flag == 3) {
        ierr=fscanf(fp,"%d ",&nface);
        //printf("nface=%d, natm=%d, den=%lf, prob=%lf\n", nface,natm,den,prob_rds);
    }

    make_matrix(face, 3, nface);

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

    make_matrix(face_copy, 3, nface);

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
//                     printf("particles %d and %d are too close: %e\n",
//                            i,j,dist_local);
                goto exit;
            }
        }

        if (area_local < 1e-5) {
//                     printf("Triangle %d has small area: %e\n",
//                            i, area_local);
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
    printf("\n%d faces have been deleted...\n", nface - nfacenew);
    nface = nfacenew;

    free_matrix(face);
    make_matrix(vars->face, 3, nface);

    for (i = 0; i < nface; i++) {
        for (j = 0; j < 3; j++) {
            vars->face[j][i] = face_copy[j][i];
        }
    }

    free_matrix(face_copy);
    
    vars->nface = nface;
    vars->nspt = nspt;
    
    remove("molecule.xyzr");
    remove("molecule.vert");
    remove("molecule.face");

    //sprintf(fname_tp, "rm -f molecule.xyzr");
    //ierr = system(fname_tp);
    //sprintf(fname_tp, "rm -f molecule.vert");
    //ierr = system(fname_tp);
    //sprintf(fname_tp, "rm -f molecule.face");
    //ierr = system(fname_tp);






/*  tr_xyz: The position of the particles on surface */
/*    tr_q: The normal direction at the particle location */
/* tr_area: The triangular area of each element */

    make_matrix(particles->particle_position, 3, vars->nface);
    make_matrix(particles->particle_normal, 3, vars->nface);
    make_vector(particles->particle_area, vars->nface);
    make_vector(particles->source_term, 2 * vars->nface);
    
    //tr_xyz = (double *) calloc(3*nface, sizeof(double));
    //tr_q = (double *) calloc(3*nface, sizeof(double));
    //tr_area = (double *) calloc(nface, sizeof(double));

    for (i = 0; i < vars->nface; i++) {
        for (j = 0; j < 3; j++) {
            idx[j] = vars->face[j][i];
        }

        for (j = 0; j < 3; j++) {
            r0[j] = 0;
            v0[j] = 0;

            for (k = 0; k < 3; k++) {
                r0[j] = r0[j] + vars->vert[j][idx[k]-1] / 3.0;
                v0[j] = v0[j] + vars->snrm[j][idx[k]-1] / 3.0;
                r[j][k] = vars->vert[j][idx[k]-1];
                v[j][k] = vars->snrm[j][idx[k]-1];
            }
        }

        v0_norm = sqrt(v0[0]*v0[0] + v0[1]*v0[1] + v0[2]*v0[2]);

        for (k = 0; k<3; k++) {
            v0[k] = v0[k] / v0_norm;
        }

        for (j = 0; j < 3; j++) {
            particles->particle_position[j][i] = r0[j];
            particles->particle_normal[j][i] = v0[j];
            //tr_xyz[3*i + j] = r0[j];
            //tr_q[3*i + j] = v0[j];
        }

        particles->particle_area[i] = TriangleArea(r);
        //tr_area[i] = TriangleArea(r);
        sum += particles->particle_area[i];
    }

    printf("Total suface area = %.17f\n",sum);

    return 0;
}
