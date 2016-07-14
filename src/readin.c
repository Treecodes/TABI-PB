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
 * Last changed at 6/29/2016
 */

#include <time.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "gl_variables.h"
#include "array.h"
#include "TABIPBstruct.h"

/* function computing the area of a triangle given vertices coodinates */
double triangle_area(double v[3][3])
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
int readin(TABIPBparm *parm, TABIPBvars *vars)
{
        FILE *fp, *wfp, *nsfp;
        char c, c1[10], c2[10], c3[10], c4[10], c5[10];
        char fname_tp[256];

        int i, j, k, i1, i2, i3, j1, j2, j3, ii, jj, kk;
        int nfacenew, ichanged, ierr;

        int namelength = 4;

        double den, prob_rds, a1, a2, a3, b1, b2, b3, a_norm, r0_norm, v0_norm;
        double r0[3], v0[3], v[3][3], r[3][3];
        int idx[3], jface[3], iface[3];
        double rs, rd[3], pot = 0.0, sum = 0.0, pot_temp = 0.0;
        double temp_x, temp_q, tchg, tpos[3], dist_local, area_local;
        double cos_theta, G0, tp1, G1, r_s[3];
        double xx[3], yy[3];

        int **face_copy;

        natm = parm->number_of_lines;    //natm is number of lines in .xyzr

  /* Run msms */
        if (parm->mesh_flag == 0) {
                //sprintf(fname_tp, "msms -if %s%s.xyzr -prob %s -dens %s -of %s%s ",
                //        parm->fpath, parm->fname, parm->probe_radius, parm->density, parm->fpath, parm->fname);
                sprintf(fname_tp, "msms -if molecule.xyzr -prob %f -dens %f -of molecule",
                        parm->probe_radius,parm->density);
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
                //fprintf(nsfp, "XYZR_FileName = %s%s.xyzr\n", parm->fpath, parm->fname);
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
                ierr = system("NanoShaper");
                //sprintf(fname_tp,"mv triangulatedSurf.face %s%s.face\n", parm->fpath, parm->fname);
                sprintf(fname_tp,"mv triangulatedSurf.face molecule.face\n");
                ierr = system(fname_tp);
                //sprintf(fname_tp,"mv triangulatedSurf.vert %s%s.vert\n", parm->fpath, parm->fname);
                sprintf(fname_tp,"mv triangulatedSurf.vert molecule.vert\n");
                ierr = system(fname_tp);
                ierr = system("rm -f stderror.txt");
                ierr = system("rm -f surfaceConfiguration.prm");
                ierr = system("rm -f triangleAreas.txt");
                ierr = system("rm -f exposed.xyz");
                ierr = system("rm -f exposedIndices.txt");
        }

  /* read in vert */
        //sprintf(fname_tp, "%s%s.vert", parm->fpath, parm->fname);
        sprintf(fname_tp, "molecule.vert");

  /* open the file and read through the first two rows */
        fp = fopen(fname_tp, "r");
        for (i = 1; i <= 2; i++) {
                while ((c = getc(fp)) != '\n') {
                }
        }


        if (parm->mesh_flag == 0) {
                ierr = fscanf(fp,"%d %d %lf %lf ",&nspt,&natm,&den,&prob_rds);
                //printf("nspt=%d, natm=%d, den=%lf, prob=%lf\n", nspt,natm,den,prob_rds);

        } else if (parm->mesh_flag == 1 ||
                   parm->mesh_flag == 2 || 
                   parm->mesh_flag == 3) {
                ierr = fscanf(fp,"%d ",&nspt);
                //printf("nspt=%d, natm=%d, den=%lf, prob=%lf\n", nspt,natm,den,prob_rds);
        }


  /*allocate variables for vertices file*/

        make_matrix(vars->extr_v, 3, nspt);
        make_matrix(vert, 3, nspt);
        make_matrix(snrm, 3, nspt);

        for (i = 0; i <= nspt-1; i++) {
                ierr = fscanf(fp, "%lf %lf %lf %lf %lf %lf %d %d %d",
                              &a1, &a2, &a3, &b1, &b2, &b3, &i1, &i2, &i3);

  /*radial projection to improve msms accuracy, ONLY FOR SPHERE!!!!!!!!
  a_norm=sqrt(a1*a1+a2*a2+a3*a3);
  b1=a1/a_norm;
  a1=b1*rds;
  b2=a2/a_norm;
  a2=b2*rds;
  b3=a3/a_norm;
  a3=b3*rds;*/

                vert[0][i] = a1;
                vert[1][i] = a2;
                vert[2][i] = a3;
                snrm[0][i] = b1;
                snrm[1][i] = b2;
                snrm[2][i] = b3;
                vars->extr_v[0][i] = i1;
                vars->extr_v[1][i] = i2;
                vars->extr_v[2][i] = i3;
        }

        fclose(fp);
        printf("Finished reading .vert file...\n");

  /* read in faces */

        //sprintf(fname_tp, "%s%s.face", parm->fpath, parm->fname);
        sprintf(fname_tp, "molecule.face");
        fp = fopen(fname_tp, "r");
        for (i = 1; i < 3; i++) {
                while ((c = getc(fp)) != '\n') {
                }
        }


        if (parm->mesh_flag == 0) {
                ierr=fscanf(fp,"%d %d %lf %lf ",&nface,&natm,&den,&prob_rds);
                //printf("nface=%d, natm=%d, den=%lf, prob=%lf\n", nface,natm,den,prob_rds);

        } else if (parm->mesh_flag == 1 || 
                   parm->mesh_flag == 2 || 
                   parm->mesh_flag == 3) {
                ierr=fscanf(fp,"%d ",&nface);
                //printf("nface=%d, natm=%d, den=%lf, prob=%lf\n", nface,natm,den,prob_rds);
        }


        make_matrix(vars->extr_f, 2, nface);
        make_matrix(face, 3, nface);


        for (i = 0; i < nface; i++) {
                ierr = fscanf(fp,"%d %d %d %d %d",&j1,&j2,&j3,&i1,&i2);
                face[0][i] = j1;
                face[1][i] = j2;
                face[2][i] = j3;
                vars->extr_f[0][i] = i1;
                vars->extr_f[1][i] = i2;
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
                                r[kk][ii] = vert[kk][iface[ii]-1];
                                xx[kk] = xx[kk] + 1.0/3.0 * r[kk][ii];
                        }
                }

                area_local = triangle_area(r);

                for (j = i-10; (j >= 0 && j < i); j++) {
                        for (jj = 0; jj < 3; jj++) {
                                jface[jj] = face[jj][j];
                                yy[jj] = 0.0;
                        }

                        for (jj = 0; jj < 3; jj++) {
                                for (kk = 0; kk < 3; kk++) {
                                        r[kk][jj] = vert[kk][jface[jj]-1];
                                        yy[kk] = yy[kk] + 1.0/3.0 * r[kk][jj];
                                }
                        }

                        dist_local = 0.0;

                        for (jj = 0; jj < 3; jj++)
                                dist_local += (xx[jj]-yy[jj]) * (xx[jj]-yy[jj]);

                        dist_local = sqrt(dist_local);

                        if (dist_local < 1e-5) {
                                printf("particles %d and %d are too close: %e\n",
                                                i,j,dist_local);
                                goto exit;
                        }
                }

                if (area_local < 1e-5) {
                        printf("Triangle %d has small area: %e\n",
                               i, area_local);
                        goto exit;
                }

                continue;
                exit: ichanged = nface - nfacenew;

                for (ii = i - ichanged; ii < nface; ii++) {
                        for (jj = 0; jj < 3; jj++)
                                face_copy[jj][ii] = face_copy[jj][ii+1];
                }

                nfacenew = nfacenew - 1;
        }
        printf("\n%d faces have been deleted...\n", nface - nfacenew);
        nface = nfacenew;


        free_matrix(face);

        make_matrix(face, 3, nface);

        for (i = 0; i < nface; i++) {
                for (j = 0; j < 3; j++)
                        face[j][i] = face_copy[j][i];
        }

        free_matrix(face_copy);

/*  tr_xyz: The position of the particles on surface */
/*    tr_q: The normal direction at the particle location */
/* tr_area: The triangular area of each element */

        tr_xyz = (double *) calloc(3*nface, sizeof(double));
        tr_q = (double *) calloc(3*nface, sizeof(double));
        tr_area = (double *) calloc(nface, sizeof(double));
        bvct = (double *) calloc(2*nface, sizeof(double));

        for (i = 0; i < nface; i++) {
                for (j = 0; j < 3; j++)
                        idx[j] = face[j][i];

                for (j = 0; j < 3; j++) {
                        r0[j] = 0;
                        v0[j] = 0;

                        for (k = 0; k < 3; k++) {
                                r0[j] = r0[j] + vert[j][idx[k]-1] / 3.0;
                                v0[j] = v0[j] + snrm[j][idx[k]-1] / 3.0;
                                r[j][k] = vert[j][idx[k]-1];
                                v[j][k] = snrm[j][idx[k]-1];
                        }
                }

                v0_norm = sqrt(v0[0]*v0[0] + v0[1]*v0[1] + v0[2]*v0[2]);

                for (k = 0; k<3; k++)
                        v0[k] = v0[k] / v0_norm;


  /* radial projection for sphere only!!!*/
  /*  r0_norm=sqrt(r0[0]*r0[0]+r0[1]*r0[1]+r0[2]*r0[2]);
    for (k=0;k<=2;k++){
      r0[k]=r0[k]/r0_norm*rds;
    }*/

                for (j = 0; j < 3; j++) {
                        tr_xyz[3*i + j] = r0[j];
                        tr_q[3*i + j] = v0[j];
                }

                tr_area[i] = triangle_area(r);
                sum += tr_area[i];
        }

        printf("Total suface area = %.17f\n",sum);

        //sprintf(fname_tp, "rm -f %s%s.xyzr", parm->fpath, parm->fname);
        sprintf(fname_tp, "rm -f molecule.xyzr");
        ierr = system(fname_tp);
        //sprintf(fname_tp, "rm -f %s%s.vert", parm->fpath, parm->fname);
        sprintf(fname_tp, "rm -f molecule.vert");
        ierr = system(fname_tp);
        //sprintf(fname_tp, "rm -f %s%s.face", parm->fpath, parm->fname);
        sprintf(fname_tp, "rm -f molecule.face");
        ierr = system(fname_tp);

        return 0;

}
