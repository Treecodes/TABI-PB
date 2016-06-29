/*
 * C routine to interface tabipb with apbs
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
 * Rebuild the architecture of wrapper
 *
 */

int main(int argc, char *argv[]){

  int tabipb(char** apbs_pqr_filename,int* nion,double* ionc,
             double* ionq, double* ionr, double* pdie,
             double* sdie, double* sdens, double* temp, double* srad,
             int* tree_order, int* tree_n0, double* mac, int* mesh);






  return 0;
}
