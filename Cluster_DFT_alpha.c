/**********************************************************************
  Cluster_DFT_alpha.c:

     Cluster_DFT_alpha.c is a subroutine to calulate the alpha
     based on cluster calculations

  Log of Cluster_DFT_Dosout.c:

     Feb/2020  Released by J.Parq

***********************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "openmx_common.h"
#include "mpi.h"

#define measure_time 0

#define UNDER1 0.7 // regareded as one if larger than this
#define OVER1 1.4  // regareded as one if less than this
#define B_LIM nm*10          // max. size of bphase

void Set_MOD_Grid(double **moc, double **qp, double *MO_Grid_D, double *MOSD_Grid_D, int nI, double sdr);
void Set_MO_VHart_Grid(double *MO_VHart_Grid_B, double *ReRhok, double *ImRhok);

void Calc_EXC_EH(double *MOD_Grid_D, double *Null_Grid_D, double *MOSD_Grid_D, double *MO_VHart_Grid_B, int nI, double **qp, double sdr, double Usic[2]);
void Calc_EX_EH(double *MOD_Grid_D, double *Null_Grid_D, double *MOSD_Grid_D, double *MO_VHart_Grid_B, int nI, double **qp, double sdr, double Usic[2]);

double angle_xc(double denR, double r, double dRdr, double *dFdGr, int i, int pi, int ni, double dr, int xc_switch);

double V_XC_GGA(double den0, double dndx, double dndy, double dndz, double F[2]);
double V_XC_PW_LSDA(double den0);

double RefAOUh(double **qp, double sdr);
double RefAOUxc(double **qp, double sdr);


void Cluster_DFT_alpha( int SpinP_switch,  double *****nh, double ****CntOLP, double *ReDenk, double *ImDenk)
{
  static int firsttime=1;
  int n,i,j,k,l,m,n2,k2;
  int *MP;
  int n1min,iemin,iemax,spin,i1,j1,k1,iemin0,iemax0,n1,ie,ie1,deg_i,deg_l;
  int GA_AN, MA_AN, wanA, tnoA, Anum, LB_AN, GB_AN, MB_AN, wanB, tnoB, Bnum;
  int LC_AN, GC_AN, wanC, tnoC, Cnum, Rnh;
  int noA,noB,m1,Acount,Bcount,Ccount;
  int nI, tl, MaxL;
  int i_vec[10];  
  int file_ptr_size;

  double EV_cut0, *ev_gap;
  double **ko; int N_ko, i_ko[10];
  double ***H; int N_H, i_H[10];
  double ***C; int N_C, i_C[10];
  double **Ctmp; int N_Ctmp, i_Ctmp[10];
  double **SD; int N_SD=2, i_SD[10];
  double TStime,TEtime;

  double sum,dum;
  float *fSD; 

  int DN,BN;
  double *LOD_Grid_D, *LOSD_Grid_D, *Null_Grid_D, *LO_VHart_Grid_B;

  int **bphase,*cphase,**tbphase;
  int *loindex;
  int *SizeOrder,*multi_m;
  int nm,nlo,nonzeros,remainder,match,unitsign,po,loop_count;
  int action,count,pcount,ncount,impact,record;
  int tmpnlo,tmpi1,tmpi2;
  double p,tp,top,q,tq,tosp,p2,tp2;     // t... = total ...
  double criterion,aupperb; // aupperb = approximate upper bound
  double Ccutoff,upperb,lowerb,middle,Cmax;   // upperb = upper bound, lowerb = lower bound
  double diff, prevdiff, prevdiff2;  // diff = upperb - lowerb, prev = previous
  double drop;
  double **b, **Q;
  double *U_t;    //  target row of the unitary matrix
  double *tmp_UtSize;  // temporary array for |Ut|
  double *cop;   // population projected on concerned orbitals
  double pdsum,psdsum,dsum,sdsum,sdr;
  double *lop,*tC,*lop0;
  double lopmin;
  double FermiF,x,max_x = 60.0;
  double loph,lopxc,copplo;  // copplo = concerned orbital population projected on localized orbital
  double oph, opxc, Ulo[2];
  double opx;
  double sie,denom; // sie = self-interaction energy
  double alpha_h,alpha_xc,average;
  double coe0,coe1,U_eff;
  double **a_target_Uh, **a_target_Uxc;   // average asymmetric

  char buf1[fp_bsize];          /* setvbuf */
  char buf2[fp_bsize];          /* setvbuf */
  char buf3[fp_bsize];          /* setvbuf */
  /*  char file_eig[YOUSO10],file_ev[YOUSO10];
      FILE *fp_eig, *fp_ev; */
  int numprocs,myid,ID,tag;
  MPI_Status stat;
  MPI_Request request;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  if (myid==Host_ID && 0<level_stdout) {
    printf("\n Cluster_DFT_alpha: start\n"); fflush(stdout);
  }

  MPI_Barrier(mpi_comm_level1);
  dtime(&TStime);

  /****************************************************
             calculation of the array size
  ****************************************************/

  n = 0;
  for (i=1; i<=atomnum; i++){
    wanA  = WhatSpecies[i];
    n  = n + Spe_Total_CNO[wanA];
  }
  n2 = n + 2;

  GA_AN = target_atom;
  wanA = WhatSpecies[GA_AN];  
  if (target_l[1] < 0) nI = 1;
  else nI = 2;

  nm = 2*target_l[0] + 1;

  /****************************************************************************
   Allocation

   int     MP[List_YOUSO[1]]
   double  ko[List_YOUSO[23]][n2]
   double  H[List_YOUSO[23]][n2][n2]
   double  C[List_YOUSO[23]][n2][n2]
   double  Ctmp[n2][n2]
   double  SD[List_YOUSO[1]][List_YOUSO[7]]
   float   fSD[List_YOUSO[7]]


  ****************************************************************************/

  MP = (int*)malloc(sizeof(int)*List_YOUSO[1]);

  N_ko=2; i_ko[0]=List_YOUSO[23]; i_ko[1]=n2;
  ko=(double**)malloc_multidimarray("double",N_ko,i_ko);

  N_H=3; i_H[0]=List_YOUSO[23]; i_H[1]=n2; i_H[2]=n2;
  H=(double***)malloc_multidimarray("double",N_H, i_H);

  N_C=3; i_C[0]=List_YOUSO[23]; i_C[1]=n2; i_C[2]=n2;
  C=(double***)malloc_multidimarray("double",N_C, i_C);

  N_Ctmp=2; i_Ctmp[0]=n2; i_Ctmp[1]=n2;
  Ctmp=(double**)malloc_multidimarray("double",N_Ctmp, i_Ctmp);

  N_SD=2; i_SD[0]=List_YOUSO[1]; i_SD[1]=List_YOUSO[7];
  SD = (double**)malloc_multidimarray("double",N_SD, i_SD);

  fSD=(float*)malloc(sizeof(float)*List_YOUSO[7]);
  /* ev_gap=(double*)malloc(sizeof(double)*n2); */

  /* grids for LOs */

  LOD_Grid_D = (double*)malloc(sizeof(double)*My_NumGridD); 
  LOSD_Grid_D = (double*)malloc(sizeof(double)*My_NumGridD); // symmetric
  Null_Grid_D = (double*)malloc(sizeof(double)*My_NumGridD); 
  for (DN=0; DN<My_NumGridD; DN++) Null_Grid_D[DN] = 0.0; 

  LOD_Grid_B = (double*)malloc(sizeof(double)*My_NumGridB_AB);
  LOSD_Grid_B = (double*)malloc(sizeof(double)*My_NumGridB_AB); // symmetric

  LO_VHart_Grid_B = (double*)malloc(sizeof(double)*My_Max_NumGridB); 
  for (i=0; i<My_Max_NumGridB; i++) LO_VHart_Grid_B[i] = 0.0;

  /* formal orbital phases for target l */

  bphase = (int**)malloc(sizeof(int*)*B_LIM);
  for (i=0; i<B_LIM; i++) {
    bphase[i] = (int*)malloc(sizeof(int)*(nm+1));
  }

  tbphase = (int**)malloc(sizeof(int*)*B_LIM);  // temporary array
  for (i=0; i<B_LIM; i++) {
    tbphase[i] = (int*)malloc(sizeof(int)*(nm+1));
  }

  cphase = (int*)malloc(sizeof(int)*(nm+1));

  /* auxiliary arrays to determine formal phases of target l */
  multi_m = (int*)malloc(sizeof(int)*(nm+1));
  SizeOrder = (int*)malloc(sizeof(int)*n);  // additional usage: E diff. order

  /* index assigning a KS orbital to a LO */
  loindex = (int*)malloc(sizeof(int)*n2);

  /* LO populations */
  lop = (double*)malloc(sizeof(double)*B_LIM);  // only occupied
  lop0 = (double*)malloc(sizeof(double)*nm*3);  // occupied + unoccupied

  tC = (double*)malloc(sizeof(double)*(nm+1)); // temporary array

  /* LO's orbital coefficients for density matrix */
  b = (double**)malloc(sizeof(double*)*atomnum);
  for (k=0; k<atomnum; k++){
    GB_AN = k+1;
    wanB = WhatSpecies[GB_AN];
    tnoB = Spe_Total_CNO[wanB];

    b[k] = (double*)malloc(sizeof(double)*tnoB);
  }

  /* density coefficeints for symmetric density */
  Q = (double**)malloc(sizeof(double*)*(Spe_MaxL_Basis[wanA]+1));
  for (l=0; l<=Spe_MaxL_Basis[wanA]; l++) {
    Q[l] = (double*)malloc(sizeof(double)*Spe_Num_Basis[wanA][l]);
  }

  U_t = (double*)malloc(sizeof(double)*n2);
  cop = (double*)malloc(sizeof(double)*n2);
  tmp_UtSize = (double*)malloc(sizeof(double)*n);

  a_target_Uh = (double**)malloc(sizeof(double*)*nI);
  for (i=0; i<nI; i++){
    a_target_Uh[i] = (double*)malloc(sizeof(double)*(nm+2));
  }

  a_target_Uxc = (double**)malloc(sizeof(double*)*nI);
  for (i=0; i<nI; i++){
    a_target_Uxc[i] = (double*)malloc(sizeof(double)*(nm+2));
  }

  if (firsttime){

    PrintMemory("alpha: LOD_Grid_B",   sizeof(double)*My_NumGridB_AB, NULL);
    PrintMemory("alpha: LOSD_Grid_B",   sizeof(double)*My_NumGridB_AB, NULL);
    PrintMemory("alpha: LO_VHart_Grid_B",    sizeof(double)*My_Max_NumGridB,    NULL);
    PrintMemory("alpha: LOD_Grid_D",   sizeof(double)*My_NumGridD,    NULL);
    PrintMemory("alpha: LOSD_Grid_D",   sizeof(double)*My_NumGridD,    NULL);
    PrintMemory("alpha: Null_Grid_D",       sizeof(double)*My_NumGridD,    NULL);

    firsttime = 0;
  }

  /****************************************************
              setting a_target_Uh & a_target_Uxc
  ****************************************************/

  for (i=0; i<nI; i++){
    tl = target_l[i];

    if (i==1 && target_l[0]==target_l[1]) i1 = 1;
    else i1 = 0;

    oph = 0.0;
    opxc = 0.0;
    upperb = 0.0;
      
    for (m=0; m<(2*tl+1); m++) {

      /* calculation of b */

      for (k=0; k<atomnum; k++){
	GB_AN = k+1;
	wanB = WhatSpecies[GB_AN];

	noB = 0;
	for (j=0; j<=Spe_MaxL_Basis[wanB]; j++){
	  for (j1=0; j1<Spe_Num_Basis[wanB][j]; j1++){

	    for (m1=0; m1<(2*j+1); m1++){

	      if (GA_AN==GB_AN && tl==j && i1==j1 && m==m1) b[k][noB] = 1.0;
	      else b[k][noB] = 0.0;
	      noB++;

	    }  /* m1 */
	    
	  } /* j1 */
	} /* j */
      } /* k */

	/* calculation of Q */

      k = GA_AN-1;
      noA = 0;

      for (l=0; l<=Spe_MaxL_Basis[wanA]; l++) {
	for (j1=0; j1<Spe_Num_Basis[wanA][l]; j1++) {

	  q = 0.0;
	  for (m1=0; m1<(2*l+1); m1++) {
	    q += b[k][noA]*b[k][noA];
	    noA++;
	  }

	  Q[l][j1] = q;

	} /* j1 */
      } /* l */

	/* generation of LOD_Grid, Null_Grid, & LO_VHart_Grid */

      Set_MOD_Grid(b, Q, LOD_Grid_D, LOSD_Grid_D, nI, 1.0);

	/* pdsum = 0.0;
	psdsum =0.0;
	for (BN=0; BN<My_NumGridB_AB; BN++){
	  pdsum += LOD_Grid_B[BN];
	  psdsum += LOSD_Grid_B[BN];
	}
	printf("[%d] (%lf %lf) ",myid,pdsum*GridVol,psdsum*GridVol); fflush(stdout); */

      for (DN=0; DN<My_NumGridD; DN++) Null_Grid_D[DN] = 0.0; 

      Set_MO_VHart_Grid(LO_VHart_Grid_B,ReDenk,ImDenk);

	/* calculation of Uh & Uxc */

      Calc_EXC_EH(LOD_Grid_D, Null_Grid_D, LOSD_Grid_D, LO_VHart_Grid_B, nI, Q,1.0,Ulo);

      oph += Ulo[0];
      opxc += Ulo[1];
      sie = Ulo[0] + Ulo[1];
      if (myid==Host_ID) printf("%d: %lf %lf = %lf\n",m,Ulo[0],Ulo[1],sie);

      a_target_Uh[i][m] = Ulo[0];
      a_target_Uxc[i][m] = Ulo[1];
      if (sie > upperb) {
	upperb = sie;
	a_target_Uh[i][nm] = Ulo[0];
	a_target_Uxc[i][nm] = Ulo[1];
      }

      Calc_EX_EH(LOD_Grid_D, Null_Grid_D, LOSD_Grid_D, LO_VHart_Grid_B, nI, Q,1.0,Ulo);
      opx += Ulo[1];
      if (myid==Host_ID) printf("%d: %lf \n",m,Ulo[1]);

    } /* m */

    U_eff = eV2Hartree*(oph+opxc)/(double)nm;
    if (myid==Host_ID) printf("%d %d %d %lf %lf SIE_av = %lf eV\n",i,i1,m,oph,opxc,U_eff);

  } /* i */

  /****************************************************
                  finding S matrix 
  ****************************************************/

  Overlap_Cluster(CntOLP,S,MP);

  if (myid==Host_ID){

    n = S[0][0];

    Eigen_lapack(S,ko[0],n,n);

    /* S[i][j] contains jth eigenvector, not ith ! */

    /****************************************************
              searching of negative eigenvalues
    ****************************************************/

    /* minus eigenvalues to 1.0e-14 */
    
    for (l=1; l<=n; l++){
      if (ko[0][l]<0.0) ko[0][l] = 1.0e-14;
      EV_S[l] = ko[0][l];
    }

    /* print to the standard output */

    if (2<=level_stdout && myid==Host_ID){
      for (l=1; l<=n; l++){
	printf(" <Cluster_DFT_alpha>  Eigenvalues of OLP  %2d  %18.15f\n",l,ko[0][l]);
      }
    }

    /* calculate S*1/sqrt(ko) */

    for (l=1; l<=n; l++){
      IEV_S[l] = 1.0/sqrt(ko[0][l]);
    }

    /*********************************************************************
     A = U^+ S U    : A diagonal
     A[n] delta_nm =  U[j][n] S[j][i] U[i][m] =U^+[n][j] S[j][i] U[i][m]
     1 = U A^-1 U^+ S
     S^-1 =U A^-1 U^+
     S^-1[i][j]= U[i][n] A^-1[n] U^+[n][j]
     S^-1[i][j]= U[i][n] A^-1[n] U[j][n]
    **********************************************************************/
 
    /*
    printf("Error check S S^{-1} =\n");
    for (i=1;i<=n;i++) {
      for (j=1;j<=n;j++) {
        sum=0.0;
        for (k=1;k<=n;k++) {
           sum+= Ovlp[i][k]*  Sinv[k][j];
        }
        printf("%lf ",sum);
      }
      printf("\n");
    }
   */

  }  /*  if (myid==Host_ID */

#if 0
      for (i=1;i<=n;i++) {
        printf("%d: ",myid);
        for (j=1;j<=n;j++) {
          sum=S[i][j];
          printf("%lf ",sum);
        }
        printf("\n");
      }

  MPI_Finalize();
  exit(0);
#endif

  /****************************************************
    Calculations of eigenstaes for up and down spins

     Note:
         MP indicates the starting position of
              atom i in arraies H and S
  ****************************************************/

  n1min=n;
  iemin=1; iemax=n;
  EV_cut0 = Threshold_OLP_Eigen;

  for (spin=0; spin<=SpinP_switch; spin++){

    Hamiltonian_Cluster(nh[spin],H[spin],MP);

    if (myid==Host_ID){

      for (i1=1; i1<=n; i1++){
        for (j1=1; j1<=n; j1++){
          sum = 0.0;
          for (l=1; l<=n; l++){
            sum = sum + H[spin][i1][l]*S[l][j1]*IEV_S[j1]; 
          }
          C[spin][i1][j1] = sum;
        }
      }

      for (i1=1; i1<=n; i1++){
        for (j1=1; j1<=n; j1++){
          sum = 0.0;
          for (l=1; l<=n; l++){
            sum = sum + IEV_S[i1]*S[l][i1]*C[spin][l][j1];
          }
          H[spin][i1][j1] = sum;
        }
      }

      /*****   H -> B_nl in the note  *****/

      for (i1=1; i1<=n; i1++){
        for (j1=1; j1<=n; j1++){
          C[spin][i1][j1] = H[spin][i1][j1];
        }
      }

      /* penalty for ill-conditioning states */

      for (i1=1; i1<=n; i1++){

        if (EV_S[i1]<EV_cut0){
          C[spin][i1][i1] += pow((EV_S[i1]/EV_cut0),-2.0) - 1.0;
        }

        /* cutoff the interaction between the ill-conditioned state */
 
        if (1.0e+3<C[spin][i1][i1]){
          for (j1=1; j1<=n; j1++){
            C[spin][i1][j1] = 0.0;
            C[spin][j1][i1] = 0.0;
	  }
          C[spin][i1][i1] = 1.0e+4;
        }
      }

      /* diagonalize the matrix */

      n1 = n;
      Eigen_lapack(C[spin],ko[spin],n1,n1);

      for (i1=1; i1<=n; i1++){
        for (j1=1; j1<=n1; j1++){
	  H[spin][i1][j1] = C[spin][i1][j1];
        }
      }

      /* print to the standard output */

      if (2<=level_stdout && myid==Host_ID){
        for (l=1; l<=n; l++){
	  printf(" <Cluster_DFT_alpha>  Eigenvalues of H spin=%2d  %2d  %18.15f\n",
                    spin,l,ko[spin][l]);
        }
      }

      /****************************************************
          Transformation to the original eigenvectors.
                        AIST NOTE 244P
      ****************************************************/

      for (i1=1; i1<=n; i1++){
        for (j1=1; j1<=n; j1++){
          C[spin][i1][j1] = 0.0;
        }
      }

      for (i1=1; i1<=n; i1++){
        for (j1=1; j1<=n1; j1++){
          sum = 0.0;
          for (l=1; l<=n; l++){
            sum = sum + S[i1][l]*IEV_S[l]*H[spin][l][j1];
          }
          C[spin][i1][j1] = sum;
        }
      }

    } /* if (myid==Host_ID) */
  }   /* spin */

  /****************************************************
    MPI:

    C
  ****************************************************/

  for (spin=0; spin<=SpinP_switch; spin++){

    for (i1=0; i1<=n; i1++){
      for (j1=0; j1<=n; j1++){
        Ctmp[i1][j1] = C[spin][i1][j1];
      }
    }

    for (i1=0; i1<=n; i1++){
      MPI_Bcast(&Ctmp[i1][0], n+1, MPI_DOUBLE, Host_ID, mpi_comm_level1);
    }

    for (i1=0; i1<=n; i1++){
      for (j1=0; j1<=n; j1++){
        C[spin][i1][j1] = Ctmp[i1][j1];
      }
    }

    MPI_Bcast(&ko[spin][0], n+1, MPI_DOUBLE, Host_ID, mpi_comm_level1);

  }

#if  0
  printf("%d: Bcast C end %d %d\n",myid,iemin,iemax);
  MPI_Barrier(mpi_comm_level1);
#endif


  /*************************************************************
     adjustment of C matrix to rotate degenerate eigenvectors

     Only double degeneracy is handled in the current version.
  *************************************************************/

  for (spin=0; spin<=SpinP_switch; spin++){

    /* GA_AN = target_atom;
       wanA = WhatSpecies[GA_AN];  
       MA_AN = F_G2M[GA_AN]; */
    Anum = MP[GA_AN]; 

    tl = target_l[0];
    noA = 0;
    for (l=0; l<tl; l++){
      for (i1=0; i1<Spe_Num_Basis[wanA][l]; i1++) {
	for (m1=0; m1<(2*l+1); m1++) noA++;
      }
    }
    Acount = Anum + noA;

    po = 0;
    for (ie=iemin; ie<=iemax; ie++){

      if (ie<n) dum = fabs(ko[spin][ie] - ko[spin][ie+1]);
      else dum = 1e4;

      /* just for safety */
      if (dum >= EV_cut0 && po == 1) po = 2;
      MPI_Reduce(&po, &pcount, 1, MPI_INT, MPI_SUM, Host_ID, mpi_comm_level1);
      if (myid==Host_ID) {
	if (pcount%numprocs != 0) {
	  printf("Inconsidtent eigenvalues!\n");
	  MPI_Finalize();
	  exit(0);
	}
      }

      if (dum < EV_cut0) {
	if (po) count++;
	else {
	  po = 1;
	  count = 2;
	  deg_i = ie;
	}
      }
      else if (po) {
	deg_l = deg_i + count;

	if (count==2) {
	  for (i1=0; i1<=n; i1++){
	    Ctmp[i1][1] = C[spin][i1][deg_i];
	    Ctmp[i1][2] = C[spin][i1][deg_i+1];
	  }

	  Cmax = 0.0;
	  for (m1=0; m1<nm; m1++) {
	    dum = fabs(Ctmp[Acount+m1][1]);
	    if (dum>Cmax) {
	      Cmax = dum;
	      action = m1;
	    }
	  }
	  tC[1] = Cmax;

	  Cmax = 0.0;
	  for (m1=0; m1<nm; m1++) {
	    dum = fabs(Ctmp[Acount+m1][2]);
	    if (dum>Cmax) {
	      Cmax = dum;
	      impact = m1;
	    }
	  }
	  tC[2] = Cmax;

	  if (tC[1]>EV_cut0 && tC[2]>EV_cut0 && action!=impact) {

	    coe0 = Ctmp[Acount+action][1];
	    if (coe0*Ctmp[Acount+impact][2]>0.0) coe1 = Ctmp[Acount+impact][1];
	    else coe1 = -Ctmp[Acount+impact][1];

	    for (i1=0; i1<=n; i1++){
	      C[spin][i1][deg_i] = coe0*Ctmp[i1][1] - coe1*Ctmp[i1][2];
	      C[spin][i1][deg_i+1] = coe1*Ctmp[i1][1] + coe0*Ctmp[i1][2];
	    }

	    for (ie1=deg_i; ie1<deg_l; ie1++) {

	      p = 0.0;
	      for (MB_AN=1; MB_AN<=Matomnum; MB_AN++){
		GB_AN = M2G[MB_AN];
		wanB = WhatSpecies[GB_AN];
		Bnum = MP[GB_AN];
		tnoB = Spe_Total_CNO[wanB];

		for (LC_AN=0; LC_AN<=FNAN[GB_AN]; LC_AN++){
		  GC_AN = natn[GB_AN][LC_AN];
		  wanC = WhatSpecies[GC_AN];
		  Cnum = MP[GC_AN];
		  tnoC = Spe_Total_CNO[wanC];

		  for(i=0; i<tnoB; i++) {
		    for (j=0; j<tnoC; j++) {
		      p += C[spin][Bnum+i][ie1]*C[spin][Cnum+j][ie1]*CntOLP[MB_AN][LC_AN][i][j];
		    }
		  }

		} /* LC_AN */
	      } /* MB_AN */

	      MPI_Allreduce(&p,
			    &tp, 
			    1, MPI_DOUBLE, MPI_SUM, 
			    mpi_comm_level1);

	      MPI_Barrier(mpi_comm_level1);

	      dum = 1.0/sqrt(tp);

	      for (GB_AN=1; GB_AN<=atomnum; GB_AN++){
		wanB = WhatSpecies[GB_AN];
		Bnum = MP[GB_AN];
		tnoB = Spe_Total_CNO[wanB];

		for(i=0; i<tnoB; i++) {
		  C[spin][Bnum+i][ie1] *= dum;
		}
	      }

	    } /* ie1 */
	    if (myid==Host_ID) printf(".");

	  } /* if */

	} /* if (count==2) */
	else {
	  printf("Warning! %d folds at %d\n",count,myid);
	}

	po = 0;
      } /* else if (po) */

    } /* ie */

    /* temporary 

    for (ie=iemin; ie<=iemax; ie++){

      p = 0.0;

      for (MA_AN=1; MA_AN<=Matomnum; MA_AN++){
	GA_AN = M2G[MA_AN];
	wanA = WhatSpecies[GA_AN];
	Anum = MP[GA_AN];
	tnoA = Spe_Total_CNO[wanA];

	for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
	  GB_AN = natn[GA_AN][LB_AN];
	  wanB = WhatSpecies[GB_AN];
	  Bnum = MP[GB_AN];
	  tnoB = Spe_Total_CNO[wanB];

	  for(i=0; i<tnoA; i++) {
	    for (j=0; j<tnoB; j++) {
	      Bcount = Bnum + j;
	      p += C[spin][Anum+i][ie]*C[spin][Bcount][ie]*CntOLP[MA_AN][LB_AN][i][j];
	    }
	  }

	} /* LB_AN 
      } /* MA_AN 

      MPI_Allreduce(&p,
		    &tp, 
		    1, MPI_DOUBLE, MPI_SUM, 
		    mpi_comm_level1);

      MPI_Barrier(mpi_comm_level1);

      if (myid==Host_ID) {


	printf("%lf     ",tp);
      }

    } /* ie */

    if (myid==Host_ID) printf("\n");

  } /* spin */


  /****************************************************
            calculation of alpha and U_eff
                      beginning
  ****************************************************/

  /* GA_AN = target_atom;
     wanA = WhatSpecies[GA_AN];  */
  MA_AN = F_G2M[GA_AN];
  Anum = MP[GA_AN];
  top = 0.0;
  oph = 0.0;
  opxc = 0.0;
  denom = 0.0;

  for (spin=0; spin<=SpinP_switch; spin++){

    tosp = 0.0;

    tl = target_l[0];
    noA = 0;
    for (l=0; l<tl; l++){
      for (i1=0; i1<Spe_Num_Basis[wanA][l]; i1++) {
	for (m1=0; m1<(2*l+1); m1++) noA++;
      }
    }
    Acount = Anum + noA;

    /* calculation of U_t & cop */

    for (ie=iemin; ie<=iemax; ie++){

      /* ev_gap

      if (ie>1) ev_gap[ie] = fabs(ko[spin][ie] - ko[spin][ie-1]);
      else ev_gap[ie] = 1.0e+4;

      if (ie<n) {
	dum = fabs(ko[spin][ie] - ko[spin][ie+1]);
	if (ev_gap[ie]>dum) ev_gap[ie] = dum;
      }

      if (ev_gap_min>ev_gap[ie]) ev_gap_min = ev_gap[ie];

      /* U_t & cop */

      p = 0.0;
      p2 = 0.0;

      if (MA_AN<=Matomnum && MA_AN>0) {
	tl = target_l[0];

	noA = 0;
	for (l=0; l<=Spe_MaxL_Basis[wanA]; l++){
	  for (i1=0; i1<Spe_Num_Basis[wanA][l]; i1++) {

	    for (m1=0; m1<(2*l+1); m1++) {

	      for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
		GB_AN = natn[GA_AN][LB_AN];
		wanB = WhatSpecies[GB_AN];
		Bnum = MP[GB_AN];
		tnoB = Spe_Total_CNO[wanB];

		for (j=0; j<tnoB; j++) {
		  Bcount = Bnum + j;
		  if (l==tl) {
		    if (i1==0) p += C[spin][Anum+noA][ie]*C[spin][Bcount][ie]*CntOLP[MA_AN][LB_AN][noA][j];
		    else p2 += C[spin][Anum+noA][ie]*C[spin][Bcount][ie]*CntOLP[MA_AN][LB_AN][noA][j];
		  }
		}
	      }
	      noA++;

	    } /* m1 */
	  } /* i1 */
	} /* l */
      } /* if */

      MPI_Allreduce(&p,
		    &tp, 
		    1, MPI_DOUBLE, MPI_SUM, 
		    mpi_comm_level1);

      MPI_Barrier(mpi_comm_level1);

      /* MPI_Allreduce(&p2,
		    &tp2, 
		    1, MPI_DOUBLE, MPI_SUM, 
		    mpi_comm_level1);
      */

      if (tp>1e-14) {

	cop[ie] = tp;
	U_t[ie] = tp;
	tosp += tp;
      }
      else {
	cop[ie] = 0.0;
	U_t[ie] = 0.0;
      }

      loindex[ie] = -1;

    } /* ie */

    criterion = largest(top,tosp)*0.01;
    if (criterion < 1e-14) criterion = 1e-14;

    /* calculation of bphase & nlo */

    Ccutoff = EV_cut0;
    upperb = 1.5;
    lowerb = 1e-14;
    diff = upperb - lowerb;
    prevdiff = 0.0;
    prevdiff2 = lowerb;
    drop = 1e-14*lonif;

    nlo = 0;
    for (m1=0; m1<nm; m1++) multi_m[m1] = 0;

    for (ie=iemin; ie<=iemax; ie++){

      if (cop[ie]>1e-14) {

	for (m1=0; m1<nm; m1++) {

	  if (C[spin][Acount+m1][ie]>Ccutoff) cphase[m1] = 1;
	  else if (C[spin][Acount+m1][ie]<-Ccutoff) cphase[m1] = -1;
	  else cphase[m1] = 0;
	}

	nonzeros = 0;
	for (m1=0; m1<nm; m1++) {
	  if (cphase[m1]!=0) nonzeros++;
	}
	cphase[nm] = nonzeros;

	j = 0;
	po = 0;
	while (j<nlo && po==0) {
	  match = 0;
	  for (m1=0; m1<nm; m1++) match += bphase[j][m1]*cphase[m1];
	  if (match<0) {
	    match = -match;
	  }

	  if (match==bphase[j][nm] && match==cphase[nm]) {
	    po++;
	    lop[j] += cop[ie];
	  }

	  j++;
	} /* while */

	if (po==0) {
	  if (cphase[nm]!=0) {
	    for (m1=0; m1<=nm; m1++) bphase[nlo][m1] = cphase[m1];
	    lop[nlo] = cop[ie];
	    multi_m[cphase[nm]]++;
	    nlo++;
	  }
	}

      }

    } /* ie */

    while (diff > 2e-14)  {

      prevdiff2 = prevdiff;
      prevdiff = diff;

      if (nlo>nm) {
	lowerb = Ccutoff;
	Ccutoff *= (double)nlo/(double)nm;
	middle = (upperb+lowerb)*0.5;
	if (Ccutoff>upperb || Ccutoff<middle) Ccutoff = middle;
      }
      else if (nlo<nm) {
	upperb = Ccutoff;
	Ccutoff *= (double)nlo/(double)nm;
	middle = (upperb+lowerb)*0.5;
	if (Ccutoff<lowerb || Ccutoff>middle) Ccutoff = middle;	
      }
      else {
	middle = (upperb+lowerb)*0.5;
	upperb = Ccutoff+5e-15;
	if (middle<upperb) Ccutoff = middle;
      }

      nlo = 0;
      for (m1=0; m1<nm; m1++) multi_m[m1] = 0;

      for (ie=iemin; ie<=iemax; ie++){

	if (cop[ie]>1e-14) {

	  for (m1=0; m1<nm; m1++) {

	    if (C[spin][Acount+m1][ie]>Ccutoff) cphase[m1] = 1;
	    else if (C[spin][Acount+m1][ie]<-Ccutoff) cphase[m1] = -1;
	    else cphase[m1] = 0;
	  }

	  nonzeros = 0;
	  for (m1=0; m1<nm; m1++) {
	    if (cphase[m1]!=0) nonzeros++;
	  }
	  cphase[nm] = nonzeros;

	  j = 0;
	  po = 0;
	  while (j<nlo && po==0) {
	    match = 0;
	    for (m1=0; m1<nm; m1++) match += bphase[j][m1]*cphase[m1];
	    if (match<0) {
	      match = -match;
	    }

	    if (match==bphase[j][nm] && match==nonzeros) {
	      po++;
	      lop[j] += cop[ie];
	    }

	    j++;
	  } /* while */

	  if (po==0) {
	    if (nonzeros!=0) {
	      for (m1=0; m1<=nm; m1++) bphase[nlo][m1] = cphase[m1];
	      lop[nlo] = cop[ie];
	      multi_m[nonzeros]++;
	      nlo++;
	    }
	  }

	}

      } /* ie */

      tmpnlo = nlo;

      count = multi_m[1];
      remainder = 0;
      m = 2;
      for (m1=2; m1<4; m1++) {
	count += multi_m[m1]/m + multi_m[m1]%m;
	remainder += multi_m[m1]%2;
	m *= 2;
      }
      for (m1=4; m1<nm; m1++) remainder += multi_m[m1]%2;
      count %= nm;
      if (count != 0) {
	  nlo -= remainder;
	  Ccutoff -= drop/2.0; 
	  lowerb -= drop;
      }           // encouraging to increase nlo

      diff = upperb - lowerb;
      if (fabs(prevdiff-diff)<1e-16 || fabs(prevdiff2-diff)<1e-16) drop -= 1e-14*lonif_decay;

    } /* while (diff > 2e-14) */

      /* reduction */

    for (i=nm-1; i>1; i--) {	
      remainder = multi_m[i]%2;

      if (remainder>0 && tmpnlo-remainder>=nm) {

	k = -1;
	lopmin = 2.0;
	for (j=0; j<nlo; j++) {
	  if (bphase[j][nm] == i && lop[j]<lopmin) {
	    lopmin = lop[j];
	    k = j;
	  }
	}
	
	if (lop[k] < 0.5) {
	  for (m1=0; m1<=nm; m1++) bphase[k][m1] = -1; // deletion
	  multi_m[i]--;
	  tmpnlo--;
	}

      } /*  if  */
    }  /* i */

    if (myid==Host_ID) {
      printf(" %d %d  %e %e %e\n",nlo,nm,Ccutoff,upperb,lowerb); 
      for (j=0; j<nlo; j++) printf(" %lf %d ",lop[j],bphase[j][nm]);
      printf("\n");
    }

    /* reorder of bphase */

    for (i=1; i<=nlo; i++) {
      nonzeros = 0;
      for (j=0; j<nlo; j++) {
	if (bphase[j][nm]>nonzeros) {
	  action = j;
	  nonzeros = bphase[j][nm];
	}
	else if (bphase[j][nm]==nonzeros) {
	  if (lop[j]<lop[action]) {
	    action = j;
	    nonzeros = bphase[j][nm];
	  }
	}
      }

      for (m1=0; m1<=nm; m1++) {
	tbphase[i][m1] = bphase[action][m1];
      }
      bphase[action][nm] = 0;
    }

    for (i=1; i<=nlo; i++) {
      lop[i] = 0.0;
      lop0[i] = 0.0;
      for (m1=0; m1<=nm; m1++) {
	bphase[i][m1] = tbphase[i][m1];
      }
    }

    nlo++;
    lop[0] = 0.0;
    lop0[0] = 0.0;
    for (m1=0; m1<=nm; m1++) {
      bphase[0][m1] = 0;
    }

    /* calculation of lop */
    
    for (ie=iemin; ie<=iemax; ie++){

      if (cop[ie]>1e-14) {

	x = (ko[spin][ie] - ChemP)*Beta;
	if (x<=-max_x) x = -max_x;
	if (max_x<=x)  x = max_x;
	FermiF = 1.0/(1.0 + exp(x));

	for (m1=0; m1<nm; m1++)  tC[m1] = C[spin][Acount+m1][ie];
	for (j=0; j<nm; j++) {
	  Cmax = 0.0;
	  action = -1;
	  for (m1=0; m1<nm; m1++) {
	    dum = fabs(tC[m1]);
	    if (dum>Cmax && dum>1e-14) {
	      Cmax = dum;
	      action = m1;
	    }
	  }
	  SizeOrder[j] = action;
	  if (action>-1) tC[action] = 0.0;
	}

	for (m1=0; m1<nm; m1++)  tC[m1] = C[spin][Acount+m1][ie];

	impact = 0;
	i = 1;
	po = 0;
	while (i<nlo && po==0) {

	  nonzeros = bphase[i][nm];
	  count = 0;

	  for (j=0; j<nonzeros; j++) {
	    action = SizeOrder[j];
	    if (action==-1 || bphase[i][action]==0) j += nonzeros;
	    else if (j==0) {
	      count++;
	      if (tC[action]>0.0) unitsign = bphase[i][action];
	      else unitsign = -bphase[i][action];
	    }
	    else {
	      match = bphase[i][action]*unitsign;
	      dum = tC[action]*(double)match;
	      if (dum>1e-14) {
		count++;
	      }
	      else j += nonzeros;
	    } 
	  } /* j */

	  if (count==nonzeros) {
	    po++;
	    loindex[ie] = i;
	    lop0[i] += cop[ie];
	    lop[i] += cop[ie]*FermiF;
	    U_t[ie] = (double)unitsign*sqrt(U_t[ie]);
	  }
	  else impact += count;

	  i++;

	} /* while */

	if (po==0 && impact>0) {

	  action = SizeOrder[0];
	  i = 1;
	  while (i<nlo && po==0 && action>-1) {

	    nonzeros = bphase[i][nm];
	    if (tC[action]>0.0) unitsign = bphase[i][action];
	    else unitsign = -bphase[i][action];

	    count = 0;
	    for (m1=0; m1<nm; m1++) {
	      if (tC[m1]*(double)(bphase[i][m1]*unitsign) > 1e-14) count++;
	    }

	    if (count==nonzeros) {
	      po++;
	      loindex[ie] = i;
	      lop0[i] += cop[ie];
	      lop[i] += cop[ie]*FermiF;
	      U_t[ie] = (double)unitsign*sqrt(U_t[ie]);
	    }

	    i++;

	  }
	    
	} /* if (po==0 && impact>0) */

	if (po==0) {
	  q = 0.0;
	  for (ie1=iemin; ie1<ie; ie1++){
	    if (loindex[ie1]==0) {
	      for (m1=0; m1<nm; m1++)
		q += U_t[ie1]*C[spin][Acount+m1][ie1]*C[spin][Acount+m1][ie];
	    }
	  }
	  if (q<0.0) unitsign = -1;
	  else unitsign = 1;

	  loindex[ie] = 0;
	  lop0[0] += cop[ie];
	  lop[0] += cop[ie]*FermiF;
	  U_t[ie] = (double)unitsign*sqrt(U_t[ie]);
	}

      } /* if (cop[ie]>1e-14) */

    } /* ie */

    if (myid==Host_ID) {
      for (j=0; j<nlo; j++) {
	printf(" %lf %d ",lop0[j],bphase[j][nm]);
      }
      printf(".\n");
    }

    /****************/
    /* merging lops */
    /****************/

    tosp = 0.0;
    for (i=1; i<nlo; i++) tosp += lop0[i];

    pcount = 0;
    for (i=1; i<nlo; i++) {
      if (lop0[i]>OVER1) pcount++;
    }

    ncount = 0;
    for (i=1; i<nlo; i++) {
      if (lop0[i]<UNDER1 && lop0[i]>1e-14) ncount++;
    }

    loop_count = 0;

    while (pcount>loop_count || ncount>0) {

      /* eleminationg too large los */
    
      for (j=0; j<pcount; j++) {
	action = -1;
	upperb = OVER1;
	for (i=1; i<nlo; i++) {
	  if (lop0[i]>upperb) {
	    action = i;
	    upperb = lop0[i];
	  }
	}

	/* reset */
	for (ie=iemin; ie<=iemax; ie++){
	  if (loindex[ie]==action) loindex[ie] = -1;
	}

	if (bphase[action][nm]>1) {

	  /* move to the last order */

	  for (m1=0; m1<=nm; m1++) tbphase[action][m1] = bphase[action][m1];

	  for (k=action+1; k<nlo; k++) {
	    lop0[k-1] = lop0[k];
	    lop[k-1] = lop[k];
	    for (m1=0; m1<=nm; m1++) bphase[k-1][m1] = bphase[k][m1];
	    for (ie=iemin; ie<=iemax; ie++){
	      if (loindex[ie]==k) loindex[ie] = k-1;
	    }
	  }

	  for (m1=0; m1<=nm; m1++) bphase[nlo-1][m1] = tbphase[action][m1];

	  /* reset */
	  lop0[nlo-1] = 0.0;
	  lop[nlo-1] = 0.0;
	}
	else {

	  /* reset */
	  lop0[action] = 0.0;
	  lop[action] = 0.0;
	}
      }

      /* eleminating too small los */

      for (j=0; j<ncount; j++) {
	action = -1;
	lowerb = UNDER1;
	for (i=1; i<nlo; i++) {
	  if (lop0[i]<lowerb && lop0[i]>1e-14) {
	    action = i;
	    lowerb = lop0[i];
	  }
	}
	k = bphase[action][nm];
	
	if (lowerb<0.5) {
	  count = 0;
	  for (i=1; i<nlo; i++) {
	    if (bphase[i][nm]==k) count++;
	  }

	  if (count>0 && multi_m[k]==0) bphase[action][nm] = -1;  // deletion
	}

	/* reset */

	lop0[action] = 0.0;
	lop[action] = 0.0;

	for (ie=iemin; ie<=iemax; ie++){
	  if (loindex[ie]==action) loindex[ie] = -1;
	}
      }

      /* recalculating lops */

      for (ie=iemin; ie<=iemax; ie++){

	if (loindex[ie]==-1) {

	  x = (ko[spin][ie] - ChemP)*Beta;
	  if (x<=-max_x) x = -max_x;
	  if (max_x<=x)  x = max_x;
	  FermiF = 1.0/(1.0 + exp(x));

	  for (m1=0; m1<nm; m1++)  tC[m1] = C[spin][Acount+m1][ie];

	  Cmax = 0.0;
	  action = -1;
	  for (m1=0; m1<nm; m1++) {
	    dum = fabs(tC[m1]);
	    if (dum>Cmax && dum>1e-14) {
	      Cmax = dum;
	      action = m1;
	    }
	  }

	  po = 0;
	  i = 1;
	  while (i<nlo && po==0 && action>-1) {

	    nonzeros = bphase[i][nm];
	    if (tC[action]>0.0) unitsign = bphase[i][action];
	    else unitsign = -bphase[i][action];

	    count = 0;
	    for (m1=0; m1<nm; m1++) {
	      if (tC[m1]*(double)(bphase[i][m1]*unitsign) > 1e-14) count++;
	    }

	    if (count==nonzeros) {
	      po++;
	      loindex[ie] = i;
	      lop0[i] += cop[ie];
	      lop[i] += cop[ie]*FermiF;
	      U_t[ie] = (double)unitsign*fabs(U_t[ie]);
	    }

	    i++;

	  }

	  if (po==0) {
	    q = 0.0;
	    for (ie1=iemin; ie1<ie; ie1++){
	      if (loindex[ie1]==0) {
		for (m1=0; m1<nm; m1++)
		  q += U_t[ie1]*C[spin][Acount+m1][ie1]*C[spin][Acount+m1][ie];
	      }
	    }
	    if (q<0.0) unitsign = -1;
	    else unitsign = 1;

	    loindex[ie] = 0;
	    lop0[0] += cop[ie];
	    lop[0] += cop[ie]*FermiF;
	    U_t[ie] = (double)unitsign*fabs(U_t[ie]);
	  }

	} /* if ((loindex[ie]==-1 || loindex[ie]==-2) */
      } /* ie */

      /* reordering */

      if (lop0[0]>1e-14) count = 1;
      else count = 0;

      for (i=1; i<nlo; i++) {

	for (j=i+1; j<nlo; j++) {
	  if (lop0[i]>lop0[j] && bphase[i][nm]==bphase[j][nm]) {
	    dum = lop0[i];
	    lop0[i] = lop0[j];
	    lop0[j] = dum;

	    dum = lop[i];
	    lop[i] = lop[j];
	    lop[j] = dum;

	    for (m1=0; m1<nm; m1++) {
	      Bnum = bphase[i][m1];
	      bphase[i][m1] = bphase[j][m1];
	      bphase[j][m1] = Bnum;
	    }

	    for (ie=iemin; ie<=iemax; ie++){
	      if (loindex[ie]==i) loindex[ie]=j;
	      else if (loindex[ie]==j) loindex[ie]=i;
	    }

	    j += nlo;
	  }
	} /* j */

	if (lop0[i]>1e-14) count++;
      } /* i */

      if (myid==Host_ID) {
	for (j=0; j<nlo; j++) {
	  printf(" %lf %d ",lop0[j],bphase[j][nm]);
	}
	printf("..\n");
      }

      pcount = 0;
      for (i=1; i<nlo; i++) {
	if (lop0[i]>OVER1) pcount++;
      }

      if (ncount==0) {
	loop_count++;
	printf("----------Warning! Incomplete characterization.\n");
      }
      else ncount = 0;

    } /*   while (pcount>loop_count || ncount>0)   */

    /* merging again */

    pcount = count;
    po = 0;

    while (count>nm) {
      for (i=1; i<nlo; i++) {

	if (lop0[i]>1e-14) {

	  action = -1;
	  for (j=i+1; j<nlo; j++) {
	    if (bphase[j][nm]!=bphase[i][nm] && lop0[j]>1e-14) {

	      match = 0;
	      for (m1=0; m1<nm; m1++) match += bphase[j][m1]*bphase[i][m1];

	      if (match<0) {
		unitsign = -1;
		match = -match;
	      }
	      else unitsign = 1;

	      if (match==bphase[j][nm] || match==bphase[i][nm]) {
		dum = fabs(lop0[i]+lop0[j]-1.0);
		if (fabs((lop0[j]-1.0)>dum && fabs(lop0[i]-1.0)>dum) || (lop0[i]+lop0[j])<1.5) {
		  action = j;
		  j += nlo;
		}
	      }

	    }
	  } /* j */

	  if (action!=-1) {
	    for (ie=iemin; ie<=iemax; ie++){
	      if (loindex[ie]==i) {
		loindex[ie] = action;
		U_t[ie] *= (double)unitsign;
	      }
	    }
	    lop0[action] += lop0[i];
	    lop[action] += lop[i];
	    lop0[i] = 0.0;
	    lop[i] = 0.0;
	    count--;
	  }
	  else if (lop0[0]>1e-14) {
	    dum = fabs(lop0[i]+lop0[0]-1.0);
	    if (fabs(lop0[0]-1.0)>dum  && fabs(lop0[i]-1.0)>dum) {

	      for (ie=iemin; ie<=iemax; ie++){
		if (loindex[ie]==i || loindex[ie]==0) {
		  q = 0.0;
		  for (ie1=iemin; ie1<ie; ie1++){
		    if (loindex[ie1]==0) {
		      for (m1=0; m1<nm; m1++)
			q += U_t[ie1]*C[spin][Acount+m1][ie1]*C[spin][Acount+m1][ie];
		    }
		  }
		  if (q<0.0) unitsign = -1;
		  else unitsign = 1;

		  loindex[ie] = 0;
		  U_t[ie] = (double)unitsign*fabs(U_t[ie]);
		}
	      }

	      lop0[0] += lop0[i];
	      lop[0] += lop[i];
	      lop0[i] = 0.0;
	      lop[i] = 0.0;
	      count--;

	    }
	  }
	  else if (lop0[i]<0.5 && bphase[i][nm]==1 && po) {

	    for (ie=iemin; ie<=iemax; ie++){
	      if (loindex[ie]==i) loindex[ie] = 0;
	    }

	    lop0[0] += lop0[i];
	    lop[0] += lop[i];
	    lop0[i] = 0.0;
	    lop[i] = 0.0;
	    po = 0;

	  } /* if & else if & else if */

	}

      } /* i */

      if (pcount==count) {
	if (po) {
	  count--;
	  po = 0;
	}
	else po = 1;
      }
      pcount = count;

    } /* while */

    if (myid==Host_ID) {
      for (j=0; j<nlo; j++) printf(" %lf %d ",lop0[j],bphase[j][nm]);
      printf("...\n");
    }

    /* just for safety */
    MPI_Bcast(&nlo, 1, MPI_INT, Host_ID, mpi_comm_level1); 
    MPI_Bcast(lop, nlo, MPI_DOUBLE, Host_ID, mpi_comm_level1); 
    for (j=0; j<nlo; j++) MPI_Bcast(bphase[j], nm+1, MPI_INT, Host_ID, mpi_comm_level1); 

    /* if (myid==Host_ID) {
      for (ie=iemin; ie<=iemax; ie++){
	for (m1=0; m1<nm; m1++) printf("%lf ",C[spin][Acount+m1][ie]);
	printf("<%d> %lf\n",loindex[ie],U_t[ie]);
      }
      }*/

    /* calculation of oph & opxc */

    tosp = 0.0;

    MPI_Barrier(mpi_comm_level1);

    for (i=0; i<nlo; i++) {

      if (lop[i] > criterion) {

	/* reordering by |Ut| */

	count = 0;
	for (ie=iemin; ie<=iemax; ie++){
	  if (loindex[ie]==i) {
	    x = (ko[spin][ie] - ChemP)*Beta;
	    if (x<=-max_x) x = -max_x;
	    if (max_x<=x)  x = max_x;
	    FermiF = 1.0/(1.0 + exp(x));
	    dum = cop[ie]*FermiF;
	    cop[ie] = dum;

	    if (dum>1e-14) {

	      po = 1;
	      m = 0;
	      while (po && m<count) {
		if (dum > tmp_UtSize[m]) {
		  for (m1=count; m1>m; m1--){
		    SizeOrder[m1] = SizeOrder[m1-1];
		    tmp_UtSize[m1] = tmp_UtSize[m1-1];
		  }
		  SizeOrder[m] = ie;
		  tmp_UtSize[m] = dum;
		  po = 0;
		}
		m++;
	      }

	      if (po) {
		SizeOrder[m] = ie;
		tmp_UtSize[m] = dum;
	      }

	      count++;
	    } /* if (dum>1e-14) */
	  }
	} /* ie */

	/* just for safety */
	MPI_Bcast(&count, 1, MPI_INT, Host_ID, mpi_comm_level1); 

	/* normalization constant 

	sum = 0.0;
	for (ie=iemin; ie<=iemax; ie++){
	  if (loindex[ie]==i) {
	    if (alpha_OCC) {
	      sum += cop[ie]; 
	    }
	    else {
	      sum += U_t[ie]*U_t[ie]; 
	    }
	  }
	}
	*/
	sum = lop[i];
	coe0 = 1.0/sqrt(sum);

	upperb = 0.0;
	copplo = 0.0;
	po = 1;
	loop_count = 0;
	aupperb = (a_target_Uh[0][nm] + a_target_Uxc[0][nm])*lop[i];

	MPI_Barrier(mpi_comm_level1);

	while (po && count>0) {

	  count--;

	  /* calculation of LO coefficents b */

	  for (k=0; k<atomnum; k++){
	    GB_AN = k+1;
	    wanB = WhatSpecies[GB_AN];
	    MB_AN = F_G2M[GB_AN];
	    Bnum = MP[GB_AN];

	    noB = 0;
	    for (j=0; j<=Spe_MaxL_Basis[wanB]; j++){
	      for (j1=0; j1<Spe_Num_Basis[wanB][j]; j1++){

		for (m1=0; m1<(2*j+1); m1++){

		  Bcount = Bnum + noB;
		  q = 0.0;
		  for (ie=iemin; ie<=iemax; ie++){
		    if (loindex[ie]==i) {
		      if (alpha_OCC) {
			x = (ko[spin][ie] - ChemP)*Beta;
			if (x<=-max_x) x = -max_x;
			if (max_x<=x)  x = max_x;
			FermiF = 1.0/(1.0 + exp(x));
			q += U_t[ie]*C[spin][Bcount][ie]*sqrt(FermiF); 
		      }
		      else {
			q += U_t[ie]*C[spin][Bcount][ie]; 
		      }
		    }
		  }
		  b[k][noB] = coe0*q;
		  noB++;

		}  /* m1 */
	    
	      } /* j1 */
	    } /* j */
	  } /* k */

	  /* calculation of target AO population on LO */

	  p = 0.0;
	  k = GA_AN-1;
	  MA_AN = F_G2M[GA_AN];

	  noA = 0;
	  for (l=0; l<=Spe_MaxL_Basis[wanA]; l++){
	    for (i1=0; i1<Spe_Num_Basis[wanA][l]; i1++){

	      for (m1=0; m1<(2*l+1); m1++){

		Acount = Anum + noA;
		q = b[k][noA];

		if (l==target_l[0]) {
		  for (LB_AN=0; LB_AN<=FNAN[GA_AN]; LB_AN++){
		    GB_AN = natn[GA_AN][LB_AN];
		    k2 = GB_AN-1;
		    wanB = WhatSpecies[GB_AN];
		    tnoB = Spe_Total_CNO[wanB];

		    if (i1==0) {
		      for (m=0; m<tnoB; m++) {
			p += q*b[k2][m]*CntOLP[MA_AN][LB_AN][noA][m];
		      }
		    }
		  }
		}

		noA++;

	      }  /* m1 */
	    } /* i1 */
	  } /* l */

	  MPI_Allreduce(&p,
			&tp, 
			1, MPI_DOUBLE, MPI_SUM, 
			mpi_comm_level1);

	  /* calculation of Q */

	  noA = 0;

	  for (l=0; l<=Spe_MaxL_Basis[wanA]; l++) {
	    for (i1=0; i1<Spe_Num_Basis[wanA][l]; i1++) {

	      q = 0.0;
	      for (m1=0; m1<(2*l+1); m1++) {
		q += b[k][noA]*b[k][noA];
		noA++;
	      }

	      Q[l][i1] = q;

	    } /* i1 */
	  } /* l */

	  Set_MOD_Grid(b, Q, LOD_Grid_D, LOSD_Grid_D, nI, 1.0);

	  pdsum = 0.0;
	  psdsum =0.0;
	  for (BN=0; BN<My_NumGridB_AB; BN++){
	    pdsum += LOD_Grid_B[BN];
	    psdsum += LOSD_Grid_B[BN];
	  }
	  /* printf("[%d] (%lf %lf) ",myid,pdsum,psdsum); fflush(stdout);*/

	  MPI_Allreduce(&pdsum,
			&dsum, 
			1, MPI_DOUBLE, MPI_SUM, 
			mpi_comm_level1);

	  MPI_Allreduce(&psdsum,
			&sdsum, 
			1, MPI_DOUBLE, MPI_SUM, 
			mpi_comm_level1);

	  MPI_Barrier(mpi_comm_level1);
	  sdr = dsum/sdsum;
	  /* printf("sdr=%lf [%d]\n",sdr,myid); fflush(stdout); */

	  /* generation of LOD_Grid, Null_Grid, & LO_VHart_Grid */

	  Set_MOD_Grid(b, Q, LOD_Grid_D, LOSD_Grid_D, nI, sdr);

	  /* psdsum =0.0;
	     for (BN=0; BN<My_NumGridB_AB; BN++){
	     psdsum += LOSD_Grid_B[BN];
	     }

	     MPI_Allreduce(&psdsum,
		      &sdsum, 
		      1, MPI_DOUBLE, MPI_SUM, 
		      mpi_comm_level1);

		      MPI_Barrier(mpi_comm_level1);

	if (myid==Host_ID) printf("  TotSymDen = %lf\n",sdsum*GridVol); */

	  Set_MO_VHart_Grid(LO_VHart_Grid_B,ReDenk,ImDenk);

	  /* calculation of Uh & Uxc */

	  Calc_EXC_EH(LOD_Grid_D, Null_Grid_D, LOSD_Grid_D, LO_VHart_Grid_B, nI, Q, sdr,Ulo);

	  sie = Ulo[0]+Ulo[1];

	  if (2<=level_stdout && myid==Host_ID){
	    printf("%lf x %lf = %lf ;  lop = %lf, aupperb = %lf\n",tp,sie,tp*sie,lop[i],aupperb);
	  }

	  /* just for safety */
	  MPI_Bcast(&sie, 1, MPI_DOUBLE, Host_ID, mpi_comm_level1); 
	  MPI_Bcast(&aupperb, 1, MPI_DOUBLE, Host_ID, mpi_comm_level1); 

	  MPI_Barrier(mpi_comm_level1);

	  if ((upperb > 1e-14 && tp < 0.5*lop[i]) || sie < 1e-14) po = 0;    // to reduce computational time and errors
	  else if ((tp*sie > aupperb && count>0) || sie <= upperb) {

	    if (loop_count>0) {
	      U_t[SizeOrder[count+1]] = -U_t[SizeOrder[count+1]];
	      U_t[SizeOrder[count]] = -U_t[SizeOrder[count]];
	    }
	    else if (count>0) U_t[SizeOrder[count]] = -U_t[SizeOrder[count]];

	  }
	  else {
	    if (myid==Host_ID) {
	      /* printing d orbital components of b */
	      k = GA_AN-1;
	      noA = 0;
	      for (l=0; l<tl; l++){
		for (i1=0; i1<Spe_Num_Basis[wanA][l]; i1++) {
		  for (m1=0; m1<(2*l+1); m1++) noA++;
		}
	      }
	      for (m1=0; m1<nm; m1++) printf("%d %lf ",bphase[i][m1],b[k][noA+m1]);
	      printf("%d lop0[%d]=%lf (%lf) SIE=%lf\n",bphase[i][nm],i,lop0[i],tp,sie); 
	    }

	    upperb = sie;
	    loph = Ulo[0];
	    lopxc = Ulo[1];
	    copplo = tp;
	    if (count>0) U_t[SizeOrder[count]] = -U_t[SizeOrder[count]];
	  }

	  loop_count++;

	}  /* while (po && count>0) */

	k = GA_AN-1;
	noA = 0;
	for (l=0; l<tl; l++){
	  for (i1=0; i1<Spe_Num_Basis[wanA][l]; i1++) {
	    for (m1=0; m1<(2*l+1); m1++) noA++;
	  }
	}

	if (upperb > 1e-14) {
	  if (alpha_OCC) {
	    oph += loph*copplo;
	    opxc += lopxc*copplo;

	    top += lop[i];
	    tosp += lop[i];

	    aupperb = 0.0;
	    sum = 0.0;
	    for (m1=0; m1<nm; m1++) {
	      dum = b[k][noA+m1]*b[k][noA+m1];
	      sum += dum;
	      aupperb += dum*(a_target_Uh[0][m1] + a_target_Uxc[0][m1]);
	    }
	    aupperb *= lop[i]/sum;
	    denom += aupperb;
	  }
	  else {
	    oph += lop0[i]*loph;
	    opxc += lop0[i]*lopxc;

	    top += lop0[i];
	    tosp += lop0[i];
	  }
	}    /*	if (upperb > 1e-14) */
	else {
	  top += lop[i];
	  tosp += lop[i];

	  aupperb = 0.0;
	  sum = 0.0;
	  for (m1=0; m1<nm; m1++) {
	    dum = b[k][noA+m1]*b[k][noA+m1];
	    sum += dum;
	    aupperb += dum*(a_target_Uh[0][m1] + a_target_Uxc[0][m1]);
	  }
	  aupperb *= lop[i]/sum;
	  denom += aupperb;
	}

	if (myid==Host_ID) printf(" [%lf %lf wp=%lf] %lf\n",loph,lopxc,lop[i],sdr);

      } /* if (lop[i]>1e-14) */

    } /* i */

    if (myid==Host_ID) printf("  spin = %lf\n",tosp);    

  } /* spin */

  /* calcultion of alpha & U_eff */

 Finishing:    if (myid==Host_ID) printf("\n");

  for (i=0; i<nI; i++){
    tl = target_l[i];

    if (i==1 && target_l[0]==target_l[1]) i1 = 1;
    else i1 = 0;

    if (top>1e-14) coe0 = 1.0/top;
    else coe0 = 0.0;

    alpha_h = coe0*oph/a_target_Uh[i][nm];
    alpha_xc = coe0*opxc/a_target_Uxc[i][nm];
    U_eff = coe0*(oph + opxc);
    if (U_eff<0.0) U_eff = 0.0;
    average = (oph + opxc)/denom;
    U_eff *= eV2Hartree;

    if (myid==Host_ID){
      printf("%lf   ",top);
      printf("%lf   %lf    ",a_target_Uh[i][nm],a_target_Uxc[i][nm]);
      printf("%lf   ",denom);
      printf("%lf   %lf\n",oph,opxc);

      printf("[l=%d, mul=%d] alpha_H: %lf   alpha_XC: %lf\n",tl,i1,alpha_h,alpha_xc);
    }
    printf(" Average alpha: %lf    U_eff = %lf eV\n",average,U_eff);

    if (myid==Host_ID) printf("\n");
  } /* i */

  /****************************************************
                          Free
  ****************************************************/

#if 0
  printf("%d: free start\n",myid); 
#endif

  free(MP);

  for (spin=0; spin<i_ko[0]; spin++){
    free(ko[spin]);
  }
  free(ko);

  /* free(ev_gap); */

  for (spin=0; spin<i_H[0]; spin++){
    for (i=0; i<i_H[1]; i++){
      free(H[spin][i]);
    }
    free(H[spin]);
  } 
  free(H);

  for (spin=0; spin<i_C[0]; spin++){
    for (i=0; i<i_C[1]; i++){
      free(C[spin][i]);
    }
    free(C[spin]);
  } 
  free(C);

  for (i=0; i<i_Ctmp[0]; i++){
    free(Ctmp[i]);
  }
  free(Ctmp);

  for (i=0; i<i_SD[0]; i++){
    free(SD[i]);
  }
  free(SD);

  free(fSD);

  free(LOD_Grid_D);
  free(LOSD_Grid_D);
  free(Null_Grid_D);

  free(LOD_Grid_B);
  free(LOSD_Grid_B);
  free(LO_VHart_Grid_B);

  free(tC);

  free(lop0);
  free(lop);
  free(loindex);

  for (j=0; j<B_LIM; j++) free(bphase[j]);
  free(bphase);

  for (j=0; j<B_LIM; j++) free(tbphase[j]);
  free(tbphase);

  free(cphase);

  free(SizeOrder);
  free(multi_m);

  for (k=0; k<atomnum; k++){
    free(b[k]);
  }
  free(b);

  for (l=0; l<=Spe_MaxL_Basis[wanA]; l++) free(Q[l]);
  free(Q);

  free(U_t);
  free(cop);
  free(tmp_UtSize);

  for (i=0; i<nI; i++){
    free(a_target_Uh[i]);
  }
  free(a_target_Uh);

  for (i=0; i<nI; i++){
    free(a_target_Uxc[i]);
  }
  free(a_target_Uxc);

#if 0
   printf("%d: alpha  Barrier start\n",myid);
   MPI_Barrier(mpi_comm_level1);
   printf("%d: alpha Barrier end\n",myid);

#endif

  /* for elapsed time */

  dtime(&TEtime);
  /*
  printf("myid=%2d Elapsed Time (s) = %15.12f\n",myid,TEtime-TStime);
  MPI_Finalize();
  exit(0);
  */
}


void Set_MOD_Grid(double **moc, double **qp, double *MOD_Grid_D, double *MOSD_Grid_D, int nI, double sdr)
{
  int al,L0,Mul0,M0,p,size1,size2;
  int Gc_AN,Mc_AN,Mh_AN,LN,AN,BN,CN;
  int DN,gp,GN_AB,GNc;
  int GA_AN,MA_AN,wanA,l,tl;
  int n1,n2,n3,k1,k2,k3,N3[4];
  int Cwan,NO0,NO1,Rn,N,Hwan,i,j,k,n;
  int NN_S,NN_R;
  unsigned long long int N2D,n2D,GN; 
  int Max_Size,My_Max;
  int size_Tmp_Den_Grid;
  int size_Den_Snd_Grid_A2B;
  int size_Den_Rcv_Grid_A2B;
  int h_AN,Gh_AN,Rnh,Nc,GRc,Nh,Nog;
  int Nc_0,Nc_1,Nc_2,Nc_3,Nh_0,Nh_1,Nh_2,Nh_3;

  double threshold;
  double tmp0,tmp1,sk1,sk2,sk3,tot_den,sum;
  double tmp0_0,tmp0_1,tmp0_2,tmp0_3;
  double sum_0,sum_1,sum_2,sum_3;
  double d1,d2,d3,cop,sip,sit,cot;
  double x,y,z,Cxyz[4];
  double dx,dy,dz,r;
  double TStime,TEtime;
  double **Tmp_Den_Grid;
  double **Den_Snd_Grid_A2B;
  double **Den_Rcv_Grid_A2B;
  double *tmp_array;
  double *tmp_array2;
  double *orbs0,*orbs1;
  double *orbs0_0,*orbs0_1,*orbs0_2,*orbs0_3;
  double *orbs1_0,*orbs1_1,*orbs1_2,*orbs1_3;
  double **tmp_CDM;
  double *Work_Array_Snd_Grid_B2D;
  double *Work_Array_Rcv_Grid_B2D;

  int *Snd_Size,*Rcv_Size;
  int numprocs,myid,tag=999,ID,IDS,IDR;
  double Stime_atom, Etime_atom;
  double time0,time1,time2;

  MPI_Status stat;
  MPI_Request request;
  MPI_Status *stat_send;
  MPI_Status *stat_recv;
  MPI_Request *request_send;
  MPI_Request *request_recv;

  /* for OpenMP */
  int OMPID,Nthrds;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);
  
  dtime(&TStime);

  /* allocation of arrays */

  size_Tmp_Den_Grid = 0;
  Tmp_Den_Grid = (double**)malloc(sizeof(double*)*(Matomnum+1)); 
  Tmp_Den_Grid[0] = (double*)malloc(sizeof(double)*1); 
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = F_M2G[Mc_AN];
    Tmp_Den_Grid[Mc_AN] = (double*)malloc(sizeof(double)*GridN_Atom[Gc_AN]);
	  
    for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
      Tmp_Den_Grid[Mc_AN][Nc] = 0.0;
    }
      
    size_Tmp_Den_Grid += GridN_Atom[Gc_AN];
  }

  size_Den_Snd_Grid_A2B = 0; 
  Den_Snd_Grid_A2B = (double**)malloc(sizeof(double*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Den_Snd_Grid_A2B[ID] = (double*)malloc(sizeof(double)*Num_Snd_Grid_A2B[ID]);
    size_Den_Snd_Grid_A2B += Num_Snd_Grid_A2B[ID];
  }  

  size_Den_Rcv_Grid_A2B = 0;   
  Den_Rcv_Grid_A2B = (double**)malloc(sizeof(double*)*numprocs);
  for (ID=0; ID<numprocs; ID++){
    Den_Rcv_Grid_A2B[ID] = (double*)malloc(sizeof(double)*Num_Rcv_Grid_A2B[ID]);
    size_Den_Rcv_Grid_A2B += Num_Rcv_Grid_A2B[ID];   
  }

  /**********************************************
              calculate LOSD_Grid_B
  ***********************************************/

  for (BN=0; BN<My_NumGridB_AB; BN++){
      LOSD_Grid_B[BN] = 0.0;
  }

  tmp0 = sdr*0.25/PI;

  GA_AN = target_atom;
  wanA = WhatSpecies[GA_AN];
  MA_AN = F_G2M[GA_AN];

  if (MA_AN <= Matomnum && MA_AN>0) {

#pragma omp parallel shared(GridN_Atom,GridListAtom,CellListAtom,atv,Gxyz,MA_AN,GA_AN,qp,tmp0,wanA,Tmp_Den_Grid,Spe_MaxL_Basis,Spe_Num_Basis) private(OMPID,Nthrds,Nc_0,GNc,GRc,Cxyz,dx,dy,dz,l,p,r,sum)
    {

      /* get info. on OpenMP */ 

      OMPID = omp_get_thread_num();
      Nthrds = omp_get_num_threads();

      for (Nc_0=OMPID*GridN_Atom[GA_AN]/Nthrds; Nc_0<(OMPID+1)*GridN_Atom[GA_AN]/Nthrds; Nc_0++){

	GNc = GridListAtom[MA_AN][Nc_0];
	GRc = CellListAtom[MA_AN][Nc_0];

	Get_Grid_XYZ(GNc,Cxyz);
	dx = Cxyz[1] + atv[GRc][1] - Gxyz[GA_AN][1];
	dy = Cxyz[2] + atv[GRc][2] - Gxyz[GA_AN][2];
	dz = Cxyz[3] + atv[GRc][3] - Gxyz[GA_AN][3];

	r = sqrt(dx*dx + dy*dy + dz*dz);
	sum = 0.0;
	for (l=0; l<=Spe_MaxL_Basis[wanA]; l++) {
	  for (p=0; p<Spe_Num_Basis[wanA][l]; p++) {
	    sum += qp[l][p]*DenRF(wanA,l,p,r);
	  }
	}
	Tmp_Den_Grid[MA_AN][Nc_0] = tmp0*sum;

      }

#pragma omp flush(Tmp_Den_Grid)

    } /* #pragma omp parallel */

  }
  
  /* copy Tmp_Den_Grid to Den_Snd_Grid_A2B */

  for (ID=0; ID<numprocs; ID++) Num_Snd_Grid_A2B[ID] = 0;
  
  N2D = Ngrid1*Ngrid2;

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    Gc_AN = M2G[Mc_AN];

    for (AN=0; AN<GridN_Atom[Gc_AN]; AN++){

      GN = GridListAtom[Mc_AN][AN];
      GN2N(GN,N3);
      n2D = N3[1]*Ngrid2 + N3[2];
      ID = (int)(n2D*(unsigned long long int)numprocs/N2D);

      Den_Snd_Grid_A2B[ID][Num_Snd_Grid_A2B[ID]] = Tmp_Den_Grid[Mc_AN][AN];

      Num_Snd_Grid_A2B[ID]++;
    }
  }    

  /* MPI: A to B */  

  request_send = malloc(sizeof(MPI_Request)*NN_A2B_S);
  request_recv = malloc(sizeof(MPI_Request)*NN_A2B_R);
  stat_send = malloc(sizeof(MPI_Status)*NN_A2B_S);
  stat_recv = malloc(sizeof(MPI_Status)*NN_A2B_R);

  NN_S = 0;
  NN_R = 0;

  tag = 999;
  for (ID=1; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (Num_Snd_Grid_A2B[IDS]!=0){
      MPI_Isend( &Den_Snd_Grid_A2B[IDS][0], Num_Snd_Grid_A2B[IDS], 
	         MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request_send[NN_S]);
      NN_S++;
    }

    if (Num_Rcv_Grid_A2B[IDR]!=0){
      MPI_Irecv( &Den_Rcv_Grid_A2B[IDR][0], Num_Rcv_Grid_A2B[IDR], 
  	         MPI_DOUBLE, IDR, tag, mpi_comm_level1, &request_recv[NN_R]);
      NN_R++;
    }
  }

  if (NN_S!=0) MPI_Waitall(NN_S,request_send,stat_send);
  if (NN_R!=0) MPI_Waitall(NN_R,request_recv,stat_recv);

  free(request_send);
  free(request_recv);
  free(stat_send);
  free(stat_recv);

  /* for myid */
  for (i=0; i<Num_Rcv_Grid_A2B[myid]; i++){
    Den_Rcv_Grid_A2B[myid][i] = Den_Snd_Grid_A2B[myid][i];
  }

  for (ID=0; ID<numprocs; ID++){
    for (LN=0; LN<Num_Rcv_Grid_A2B[ID]; LN++){

      BN = Index_Rcv_Grid_A2B[ID][3*LN+0];      
      LOSD_Grid_B[BN] += Den_Rcv_Grid_A2B[ID][LN];

    } /* LN */ 
  } /* ID */  

  /**********************************************
              calculate Tmp_Den_Grid
  ***********************************************/
    
  dtime(&time1);
  
  /* AITUNE */
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    Gc_AN = F_M2G[Mc_AN];
	  
    for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
      Tmp_Den_Grid[Mc_AN][Nc] = 0.0;
    }
  }
  
  /* AITUNE ========================== */ 
  int OneD_Nloop = 0;
  int ai_MaxNc = 0;
  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){
    int Gc_AN = M2G[Mc_AN];    
    for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){
      OneD_Nloop++;
      if(ai_MaxNc < GridN_Atom[Gc_AN]) {ai_MaxNc = GridN_Atom[Gc_AN];}
    }
  }  
  /* ai_MaxNc is maximum of GridN_Atom[] */

  int gNthrds;
#pragma omp parallel
  {
    gNthrds = omp_get_num_threads();
  }

  double** ai_tmpDG_all = (double**)malloc(sizeof(double*)*gNthrds);
	
  /* ========================== AITUNE */ 

#pragma omp parallel shared(myid,G2ID,Orbs_Grid_FNAN,List_YOUSO,time_per_atom,Tmp_Den_Grid,Orbs_Grid,COrbs_Grid,GListTAtoms2,GListTAtoms1,NumOLG,moc,WhatSpecies,ncn,F_G2M,natn,Spe_Total_CNO,M2G,MOD_Grid_D) private(OMPID,Nthrds,Mc_AN,h_AN,Stime_atom,Etime_atom,Gc_AN,Cwan,NO0,Gh_AN,Mh_AN,Rnh,Hwan,NO1,i,j,tmp_CDM,Nog,Nc_0,Nc_1,Nc_2,Nc_3,Nh_0,Nh_1,Nh_2,Nh_3,orbs0_0,orbs0_1,orbs0_2,orbs0_3,orbs1_0,orbs1_1,orbs1_2,orbs1_3,sum_0,sum_1,sum_2,sum_3,tmp0_0,tmp0_1,tmp0_2,tmp0_3,Nc,Nh,orbs0,orbs1,sum,tmp0)
  {

    orbs0 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    orbs1 = (double*)malloc(sizeof(double)*List_YOUSO[7]);

    orbs0_0 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    orbs0_1 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    orbs0_2 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    orbs0_3 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    orbs1_0 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    orbs1_1 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    orbs1_2 = (double*)malloc(sizeof(double)*List_YOUSO[7]);
    orbs1_3 = (double*)malloc(sizeof(double)*List_YOUSO[7]);

    tmp_CDM = (double**)malloc(sizeof(double*)*List_YOUSO[7]); 
    for (j=0; j<List_YOUSO[7]; j++){
      tmp_CDM[j] = (double*)malloc(sizeof(double)*List_YOUSO[7]); 
    }

    /* get info. on OpenMP */ 

    OMPID = omp_get_thread_num();
    Nthrds = omp_get_num_threads();

	
    /* AITUNE ========================== */  


    double *ai_tmpDGs;
    {
      ai_tmpDGs = (double*)malloc(sizeof(double)* ai_MaxNc);
    }
    ai_tmpDG_all[OMPID] = ai_tmpDGs;
    /* ==================================== AITUNE */


    /* for (Mc_AN=(OMPID+1); Mc_AN<=Matomnum; Mc_AN+=Nthrds){ AITUNE */
    for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

      dtime(&Stime_atom);

      /* set data on Mc_AN */

      Gc_AN = M2G[Mc_AN];
      Cwan = WhatSpecies[Gc_AN];
      NO0 = Spe_Total_CNO[Cwan]; 
	  
      int Nc;
      for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
	ai_tmpDGs[Nc] = 0.0;
      }

      for (h_AN=0; h_AN<=FNAN[Gc_AN]; h_AN++){

	/* set data on h_AN */
    
	Gh_AN = natn[Gc_AN][h_AN];
	Mh_AN = F_G2M[Gh_AN];
	Rnh = ncn[Gc_AN][h_AN];
	Hwan = WhatSpecies[Gh_AN];
	NO1 = Spe_Total_CNO[Hwan];

	/* store CDM into tmp_CDM */

	for (i=0; i<NO0; i++){
	  for (j=0; j<NO1; j++){
	    tmp_CDM[i][j] = moc[Gc_AN-1][i]*moc[Gh_AN-1][j];
	  }
	}

	/* summation of non-zero elements */
	/* for (Nog=0; Nog<NumOLG[Mc_AN][h_AN]; Nog++){ */
#pragma omp for
	for (Nog=0; Nog<NumOLG[Mc_AN][h_AN]-3; Nog+=4){

	  Nc_0 = GListTAtoms1[Mc_AN][h_AN][Nog];
	  Nc_1 = GListTAtoms1[Mc_AN][h_AN][Nog+1];
	  Nc_2 = GListTAtoms1[Mc_AN][h_AN][Nog+2];
	  Nc_3 = GListTAtoms1[Mc_AN][h_AN][Nog+3];
	  
	  Nh_0 = GListTAtoms2[Mc_AN][h_AN][Nog];
	  Nh_1 = GListTAtoms2[Mc_AN][h_AN][Nog+1];
	  Nh_2 = GListTAtoms2[Mc_AN][h_AN][Nog+2];
	  Nh_3 = GListTAtoms2[Mc_AN][h_AN][Nog+3];
	  
	  for (i=0; i<NO0; i++){
	    orbs0_0[i] = Orbs_Grid[Mc_AN][Nc_0][i];
	    orbs0_1[i] = Orbs_Grid[Mc_AN][Nc_1][i];
	    orbs0_2[i] = Orbs_Grid[Mc_AN][Nc_2][i];
	    orbs0_3[i] = Orbs_Grid[Mc_AN][Nc_3][i]; 
	  }

	  if (G2ID[Gh_AN]==myid){
	    for (j=0; j<NO1; j++){
	      orbs1_0[j] = Orbs_Grid[Mh_AN][Nh_0][j];
	      orbs1_1[j] = Orbs_Grid[Mh_AN][Nh_1][j];
	      orbs1_2[j] = Orbs_Grid[Mh_AN][Nh_2][j];
	      orbs1_3[j] = Orbs_Grid[Mh_AN][Nh_3][j]; 
	    }
	  }
	  else{
	    for (j=0; j<NO1; j++){
	      orbs1_0[j] = Orbs_Grid_FNAN[Mc_AN][h_AN][Nog  ][j];
	      orbs1_1[j] = Orbs_Grid_FNAN[Mc_AN][h_AN][Nog+1][j];
	      orbs1_2[j] = Orbs_Grid_FNAN[Mc_AN][h_AN][Nog+2][j];
	      orbs1_3[j] = Orbs_Grid_FNAN[Mc_AN][h_AN][Nog+3][j]; 
	    }
	  }
	  
	  /* Tmp_Den_Grid */

	  sum_0 = 0.0;
	  sum_1 = 0.0;
	  sum_2 = 0.0;
	  sum_3 = 0.0;

	  for (i=0; i<NO0; i++){

	    tmp0_0 = 0.0;
	    tmp0_1 = 0.0;
	    tmp0_2 = 0.0;
	    tmp0_3 = 0.0;

	    for (j=0; j<NO1; j++){
	      tmp0_0 += orbs1_0[j]*tmp_CDM[i][j];
	      tmp0_1 += orbs1_1[j]*tmp_CDM[i][j];
	      tmp0_2 += orbs1_2[j]*tmp_CDM[i][j];
	      tmp0_3 += orbs1_3[j]*tmp_CDM[i][j];
	    }

	    sum_0 += orbs0_0[i]*tmp0_0;
	    sum_1 += orbs0_1[i]*tmp0_1;
	    sum_2 += orbs0_2[i]*tmp0_2;
	    sum_3 += orbs0_3[i]*tmp0_3;
	  }
		
	  ai_tmpDGs[Nc_0] += sum_0;
	  ai_tmpDGs[Nc_1] += sum_1;
	  ai_tmpDGs[Nc_2] += sum_2;
	  ai_tmpDGs[Nc_3] += sum_3;

	} /* Nog */

#pragma omp for
	for (Nog = NumOLG[Mc_AN][h_AN] - (NumOLG[Mc_AN][h_AN] % 4); Nog<NumOLG[Mc_AN][h_AN]; Nog++){
	
	  Nc = GListTAtoms1[Mc_AN][h_AN][Nog];
	  Nh = GListTAtoms2[Mc_AN][h_AN][Nog]; 
 
	  for (i=0; i<NO0; i++){
	    orbs0[i] = Orbs_Grid[Mc_AN][Nc][i];
	  }

	  if (G2ID[Gh_AN]==myid){
	    for (j=0; j<NO1; j++){
	      orbs1[j] = Orbs_Grid[Mh_AN][Nh][j];
	    }
	  }
	  else{
	    for (j=0; j<NO1; j++){
	      orbs1[j] = Orbs_Grid_FNAN[Mc_AN][h_AN][Nog][j];
	    }
	  }
 
	  sum = 0.0;
	  for (i=0; i<NO0; i++){
	    tmp0 = 0.0;
	    for (j=0; j<NO1; j++){
	      tmp0 += orbs1[j]*tmp_CDM[i][j];
	    }
	    sum += orbs0[i]*tmp0;
	  }
 
	  ai_tmpDGs[Nc] += sum;

	} /* Nog */
	
      } /* h_AN */

      /* AITUNE   merge temporary buffer for all omp threads */	
#pragma omp for
      for (Nc=0; Nc<GridN_Atom[Gc_AN]; Nc++){
	double sum = 0.0;
	int th;
	for(th = 0; th < Nthrds; th++){
	  sum += ai_tmpDG_all[th][Nc];
	}
	Tmp_Den_Grid[Mc_AN][Nc] += sum;
      }
		

      dtime(&Etime_atom);
      time_per_atom[Gc_AN] += Etime_atom - Stime_atom;

    } /* Mc_AN */

    /* freeing of arrays */ 

    free(orbs0);
    free(orbs1);

    free(orbs0_0);
    free(orbs0_1);
    free(orbs0_2);
    free(orbs0_3);
    free(orbs1_0);
    free(orbs1_1);
    free(orbs1_2);
    free(orbs1_3);

    for (j=0; j<List_YOUSO[7]; j++){
      free(tmp_CDM[j]);
    }
    free(tmp_CDM);
	
    free(ai_tmpDGs); /* AITUNE */

#pragma omp flush(Tmp_Den_Grid)

  } /* #pragma omp parallel */
  
  free(ai_tmpDG_all);

  dtime(&time2);
  if(myid==0 && measure_time){
    printf("Time for Part1=%18.5f\n",(time2-time1));fflush(stdout);
  }

  /******************************************************
      MPI communication from the partitions A to B 
  ******************************************************/
  
  /* copy Tmp_Den_Grid to Den_Snd_Grid_A2B */

  for (ID=0; ID<numprocs; ID++) Num_Snd_Grid_A2B[ID] = 0;
  
  N2D = Ngrid1*Ngrid2;

  for (Mc_AN=1; Mc_AN<=Matomnum; Mc_AN++){

    Gc_AN = M2G[Mc_AN];

    for (AN=0; AN<GridN_Atom[Gc_AN]; AN++){

      GN = GridListAtom[Mc_AN][AN];
      GN2N(GN,N3);
      n2D = N3[1]*Ngrid2 + N3[2];
      ID = (int)(n2D*(unsigned long long int)numprocs/N2D);

      Den_Snd_Grid_A2B[ID][Num_Snd_Grid_A2B[ID]] = Tmp_Den_Grid[Mc_AN][AN];

      Num_Snd_Grid_A2B[ID]++;
    }
  }    

  /* MPI: A to B */  

  request_send = malloc(sizeof(MPI_Request)*NN_A2B_S);
  request_recv = malloc(sizeof(MPI_Request)*NN_A2B_R);
  stat_send = malloc(sizeof(MPI_Status)*NN_A2B_S);
  stat_recv = malloc(sizeof(MPI_Status)*NN_A2B_R);

  NN_S = 0;
  NN_R = 0;

  tag = 999;
  for (ID=1; ID<numprocs; ID++){

    IDS = (myid + ID) % numprocs;
    IDR = (myid - ID + numprocs) % numprocs;

    if (Num_Snd_Grid_A2B[IDS]!=0){
      MPI_Isend( &Den_Snd_Grid_A2B[IDS][0], Num_Snd_Grid_A2B[IDS], 
	         MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request_send[NN_S]);
      NN_S++;
    }

    if (Num_Rcv_Grid_A2B[IDR]!=0){
      MPI_Irecv( &Den_Rcv_Grid_A2B[IDR][0], Num_Rcv_Grid_A2B[IDR], 
  	         MPI_DOUBLE, IDR, tag, mpi_comm_level1, &request_recv[NN_R]);
      NN_R++;
    }
  }

  if (NN_S!=0) MPI_Waitall(NN_S,request_send,stat_send);
  if (NN_R!=0) MPI_Waitall(NN_R,request_recv,stat_recv);

  free(request_send);
  free(request_recv);
  free(stat_send);
  free(stat_recv);

  /* for myid */
  for (i=0; i<Num_Rcv_Grid_A2B[myid]; i++){
    Den_Rcv_Grid_A2B[myid][i] = Den_Snd_Grid_A2B[myid][i];
  }

  /******************************************************
   superposition of rho_i to calculate charge density 
   in the partition B.
  ******************************************************/

  /* initialize arrays */

  for (BN=0; BN<My_NumGridB_AB; BN++){
      LOD_Grid_B[BN] = 0.0;
  }

  /* superposition of densities rho_i */

  for (ID=0; ID<numprocs; ID++){

    for (LN=0; LN<Num_Rcv_Grid_A2B[ID]; LN++){

      BN    = Index_Rcv_Grid_A2B[ID][3*LN+0];      
      Gc_AN = Index_Rcv_Grid_A2B[ID][3*LN+1];        
      GRc   = Index_Rcv_Grid_A2B[ID][3*LN+2]; 

      if (Solver!=4 || (Solver==4 && atv_ijk[GRc][1]==0 )){

	/* spin collinear non-polarization */
	LOD_Grid_B[BN] += Den_Rcv_Grid_A2B[ID][LN];

      } /* if (Solve!=4.....) */           

    } /* AN */ 
  } /* ID */  

  /* freeing of arrays */

  for (Mc_AN=0; Mc_AN<=Matomnum; Mc_AN++){
    free(Tmp_Den_Grid[Mc_AN]);
  }
  free(Tmp_Den_Grid);

  for (ID=0; ID<numprocs; ID++){
    free(Den_Snd_Grid_A2B[ID]);
  }  
  free(Den_Snd_Grid_A2B);

  for (ID=0; ID<numprocs; ID++){
    free(Den_Rcv_Grid_A2B[ID]);
  }
  free(Den_Rcv_Grid_A2B);

  /******************************************************
             MPI: from the partitions B to D
  ******************************************************/

  /* allocation of arrays */
  
  Work_Array_Snd_Grid_B2D = (double*)malloc(sizeof(double)*GP_B2D_S[NN_B2D_S]); 
  Work_Array_Rcv_Grid_B2D = (double*)malloc(sizeof(double)*GP_B2D_R[NN_B2D_R]); 

  /******************************************************
                     MPI: LOD_Grid
  ******************************************************/

  request_send = malloc(sizeof(MPI_Request)*NN_B2D_S);
  request_recv = malloc(sizeof(MPI_Request)*NN_B2D_R);
  stat_send = malloc(sizeof(MPI_Status)*NN_B2D_S);
  stat_recv = malloc(sizeof(MPI_Status)*NN_B2D_R);

  NN_S = 0;
  NN_R = 0;

  /* MPI_Irecv */

  for (ID=0; ID<NN_B2D_R; ID++){

    IDR = ID_NN_B2D_R[ID];
    gp = GP_B2D_R[ID];

    if (IDR!=myid){ 
      MPI_Irecv( &Work_Array_Rcv_Grid_B2D[gp], Num_Rcv_Grid_B2D[IDR],
                 MPI_DOUBLE, IDR, tag, mpi_comm_level1, &request_recv[NN_R]);
      NN_R++;
    }
  }

  /* MPI_Isend */

  for (ID=0; ID<NN_B2D_S; ID++){

    IDS = ID_NN_B2D_S[ID];
    gp = GP_B2D_S[ID];

    /* copy LOD_Grid_B to Work_Array_Snd_Grid_B2D */

    for (LN=0; LN<Num_Snd_Grid_B2D[IDS]; LN++){

      BN = Index_Snd_Grid_B2D[IDS][LN];

      Work_Array_Snd_Grid_B2D[gp+LN]       = LOD_Grid_B[BN];

    } /* LN */        

    if (IDS!=myid){
      MPI_Isend( &Work_Array_Snd_Grid_B2D[gp], Num_Snd_Grid_B2D[IDS], 
		 MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request_send[NN_S]);
      NN_S++;
    }
  }

  /* MPI_Waitall */

  if (NN_S!=0) MPI_Waitall(NN_S,request_send,stat_send);
  if (NN_R!=0) MPI_Waitall(NN_R,request_recv,stat_recv);

  free(request_send);
  free(request_recv);
  free(stat_send);
  free(stat_recv);

  /* copy Work_Array_Rcv_Grid_B2D to LOD_Grid_D */

  for (ID=0; ID<NN_B2D_R; ID++){

    IDR = ID_NN_B2D_R[ID];

    if (IDR==myid){

      gp = GP_B2D_S[ID];

      for (LN=0; LN<Num_Rcv_Grid_B2D[IDR]; LN++){

	DN = Index_Rcv_Grid_B2D[IDR][LN];

	MOD_Grid_D[DN] = Work_Array_Snd_Grid_B2D[gp+LN];

      } /* LN */   

    }

    else{

      gp = GP_B2D_R[ID];

      for (LN=0; LN<Num_Rcv_Grid_B2D[IDR]; LN++){

	DN = Index_Rcv_Grid_B2D[IDR][LN];

	MOD_Grid_D[DN] = Work_Array_Rcv_Grid_B2D[gp+LN];
      }

    }
  }
  /******************************************************
                     MPI: LOSD_Grid
  ******************************************************/

  request_send = malloc(sizeof(MPI_Request)*NN_B2D_S);
  request_recv = malloc(sizeof(MPI_Request)*NN_B2D_R);
  stat_send = malloc(sizeof(MPI_Status)*NN_B2D_S);
  stat_recv = malloc(sizeof(MPI_Status)*NN_B2D_R);

  NN_S = 0;
  NN_R = 0;

  /* MPI_Irecv */

  for (ID=0; ID<NN_B2D_R; ID++){

    IDR = ID_NN_B2D_R[ID];
    gp = GP_B2D_R[ID];

    if (IDR!=myid){ 
      MPI_Irecv( &Work_Array_Rcv_Grid_B2D[gp], Num_Rcv_Grid_B2D[IDR],
                 MPI_DOUBLE, IDR, tag, mpi_comm_level1, &request_recv[NN_R]);
      NN_R++;
    }
  }

  /* MPI_Isend */

  for (ID=0; ID<NN_B2D_S; ID++){

    IDS = ID_NN_B2D_S[ID];
    gp = GP_B2D_S[ID];

    /* copy LOSD_Grid_B to Work_Array_Snd_Grid_B2D */

    for (LN=0; LN<Num_Snd_Grid_B2D[IDS]; LN++){

      BN = Index_Snd_Grid_B2D[IDS][LN];

      Work_Array_Snd_Grid_B2D[gp+LN]       = LOSD_Grid_B[BN];

    } /* LN */        

    if (IDS!=myid){
      MPI_Isend( &Work_Array_Snd_Grid_B2D[gp], Num_Snd_Grid_B2D[IDS], 
		 MPI_DOUBLE, IDS, tag, mpi_comm_level1, &request_send[NN_S]);
      NN_S++;
    }
  }

  /* MPI_Waitall */

  if (NN_S!=0) MPI_Waitall(NN_S,request_send,stat_send);
  if (NN_R!=0) MPI_Waitall(NN_R,request_recv,stat_recv);

  free(request_send);
  free(request_recv);
  free(stat_send);
  free(stat_recv);

  /* copy Work_Array_Rcv_Grid_B2D to LOSD_Grid_D */

  for (ID=0; ID<NN_B2D_R; ID++){

    IDR = ID_NN_B2D_R[ID];

    if (IDR==myid){

      gp = GP_B2D_S[ID];

      for (LN=0; LN<Num_Rcv_Grid_B2D[IDR]; LN++){

	DN = Index_Rcv_Grid_B2D[IDR][LN];

	MOSD_Grid_D[DN] = Work_Array_Snd_Grid_B2D[gp+LN];

      } /* LN */   

    }

    else{

      gp = GP_B2D_R[ID];

      for (LN=0; LN<Num_Rcv_Grid_B2D[IDR]; LN++){

	DN = Index_Rcv_Grid_B2D[IDR][LN];

	MOSD_Grid_D[DN] = Work_Array_Rcv_Grid_B2D[gp+LN];
      }

    }
  }

  /* freeing of arrays */

  free(Work_Array_Snd_Grid_B2D);
  free(Work_Array_Rcv_Grid_B2D);

}


void Set_MO_VHart_Grid(double *MO_VHart_Grid_B, double *ReRhok, double *ImRhok)
{
  int k1,k2,k3;
  int N2D,GNs,GN,BN_CB;
  int N3[4];
  double tmp0,sk1,sk2,sk3;
  double Gx,Gy,Gz,fac_invG2;
  double TStime,TEtime,etime;
  int numprocs,myid,tag=999,ID,IDS,IDR;
  MPI_Status stat;
  MPI_Request request;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  MPI_Barrier(mpi_comm_level1);
  dtime(&TStime);

  /****************************************************
            FFT of difference charge density 
  ****************************************************/

  etime = FFT_Density(5,ReRhok,ImRhok);

  /****************************************************
                       4*PI/G2/N^3
  ****************************************************/

  tmp0 = 4.0*PI/(double)(Ngrid1*Ngrid2*Ngrid3);

  N2D = Ngrid3*Ngrid2;
  GNs = ((myid*N2D+numprocs-1)/numprocs)*Ngrid1;

  for (BN_CB=0; BN_CB<My_NumGridB_CB; BN_CB++){

    GN = BN_CB + GNs;     
    k3 = GN/(Ngrid2*Ngrid1);    
    k2 = (GN - k3*Ngrid2*Ngrid1)/Ngrid1;
    k1 = GN - k3*Ngrid2*Ngrid1 - k2*Ngrid1; 

    if (k1<Ngrid1/2) sk1 = (double)k1;
    else             sk1 = (double)(k1 - Ngrid1);

    if (k2<Ngrid2/2) sk2 = (double)k2;
    else             sk2 = (double)(k2 - Ngrid2);

    if (k3<Ngrid3/2) sk3 = (double)k3;
    else             sk3 = (double)(k3 - Ngrid3);

    Gx = sk1*rtv[1][1] + sk2*rtv[2][1] + sk3*rtv[3][1];
    Gy = sk1*rtv[1][2] + sk2*rtv[2][2] + sk3*rtv[3][2]; 
    Gz = sk1*rtv[1][3] + sk2*rtv[2][3] + sk3*rtv[3][3];
    fac_invG2 = tmp0/(Gx*Gx + Gy*Gy + Gz*Gz);

    if (k1==0 && k2==0 && k3==0){
      ReRhok[BN_CB] = 0.0;
      ImRhok[BN_CB] = 0.0;
    }
    else{
      ReRhok[BN_CB] *= fac_invG2;
      ImRhok[BN_CB] *= fac_invG2;
    }
  }  

  /****************************************************
        find the Hartree potential in real space
  ****************************************************/
  
  Get_Value_inReal(0,MO_VHart_Grid_B,MO_VHart_Grid_B,ReRhok,ImRhok);

  /* for time */
  MPI_Barrier(mpi_comm_level1);
  dtime(&TEtime);

}


void Calc_EXC_EH(double *MOD_Grid_D, double *Null_Grid_D, double *MOSD_Grid_D, double *MO_VHart_Grid_B, int nI, double **qp, double sdr, double Usic[2])
{
  /****************************************************
     EXC = \sum_{\sigma}
           n_{\sigma}(\epsilon_{xc}-\mu_{xc,\sigma})

     EH = 1/2\int n(r) V_H dr
  ****************************************************/

  static int firsttime=1;
  int i,XC_P_switch,xcs;
  int n,n1,n2,n3,Ng1,Ng2,Ng3,j,k;
  int GN,GNs,BN,DN,LN,N2D,n2D,N3[4];
  double EXC,EH;
  double My_EXC,My_EH;
  double sum,tmp1,tmp2;
  int numprocs,myid,tag=999,ID,IDS,IDR;
  double sden,aden;

  MPI_Status stat;
  MPI_Request request;
  /* for OpenMP */
  int OMPID,Nthrds,Nthrds0,Nprocs,Nloop;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /****************************************************
   set Vxc_Grid
  ****************************************************/

  XC_P_switch = 4;
  xcs = XC_switch;
  /* if (xcs==4) xcs = 3; */

  Set_XC_Grid(XC_P_switch,xcs,
	      MOD_Grid_D,Null_Grid_D,
	      Null_Grid_D,Null_Grid_D,
	      Vxc_Grid_D[0], Vxc_Grid_D[1],
	      Null_Grid_D, Null_Grid_D,
	      NULL,NULL);

  /****************************************************
             copy Vxc_Grid_D to Vxc_Grid_B
  ****************************************************/

  Ng1 = Max_Grid_Index_D[1] - Min_Grid_Index_D[1] + 1;
  Ng2 = Max_Grid_Index_D[2] - Min_Grid_Index_D[2] + 1;
  Ng3 = Max_Grid_Index_D[3] - Min_Grid_Index_D[3] + 1;

  for (n=0; n<Num_Rcv_Grid_B2D[myid]; n++){
    DN = Index_Rcv_Grid_B2D[myid][n];
    BN = Index_Snd_Grid_B2D[myid][n];

    i = DN/(Ng2*Ng3);
    j = (DN-i*Ng2*Ng3)/Ng3;
    k = DN - i*Ng2*Ng3 - j*Ng3; 

    if ( !(i<=1 || (Ng1-2)<=i || j<=1 || (Ng2-2)<=j || k<=1 || (Ng3-2)<=k)){
      Vxc_Grid_B[0][BN] = Vxc_Grid_D[0][DN];
    }
  }

  /****************************************************
   set RefVxc_Grid from LOSD_Grid 
  ****************************************************/

  if (xcs==3 || xcs==4) {

    for (BN=0; BN<My_NumGridB_AB; BN++){
      RefVxc_Grid_B[BN] = V_XC_PW_LSDA(LOSD_Grid_B[BN]);
    }
  }

  /****************************************************
        calculations of Eva, Eef, EH1, and EXC
  ****************************************************/

  My_EH = 0.0;
  My_EXC = 0.0;
  sum = 0.0;

  for (BN=0; BN<My_NumGridB_AB; BN++){

    sden = LOD_Grid_B[BN];
    aden = LOSD_Grid_B[BN];

    /* EH = \int dn(r) dV_H dr + 2\int s(r) dV_H dr */
    My_EH += (sden+aden)*MO_VHart_Grid_B[BN];

    /* EXC = \sum_{\sigma} n_{\sigma}\epsilon_{xc}
       - n_{spherical}\epsilon_{xc}(n_{spherical})  */
    My_EXC += sden*Vxc_Grid_B[0][BN] - aden*RefVxc_Grid_B[BN];

    /* sum += aden*RefVxc_Grid_B[BN]; */
  } /* BN */

  /* printf("[%d] TotDen = %lf\n",myid,sum*GridVol); */

  /****************************************************
       multiplying GridVol and adding references
  ****************************************************/

  /* EH = \int dn(r) dV_H dr + 2\int s(r) dV_H dr */
  My_EH *= GridVol;

  /*   EXC = \sum_{\sigma} n_{\sigma}\epsilon_{xc}
       - n_{spherical}\epsilon_{xc}(n_{spherical}) */
  My_EXC *= GridVol;

  /* sum *= GridVol; */

  /****************************************************
   MPI:

   EH, EXC
  ****************************************************/

  MPI_Barrier(mpi_comm_level1);
  MPI_Allreduce(&My_EH, &EH, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
  MPI_Allreduce(&My_EXC, &EXC, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

  tmp1 = RefAOUh(qp,sdr);
  tmp2 = RefAOUxc(qp,sdr);

  EH += tmp1;
  EXC += tmp2;

  /* if (myid==Host_ID) printf(" RefUh = %lf (%lf), RefUxc = %lf (%lf)\n",tmp1, RefAOUh(qp,1.0), tmp2, RefAOUxc(qp,1.0)); */

  Usic[0] = EH;
  Usic[1] = EXC;
}


void Calc_EX_EH(double *MOD_Grid_D, double *Null_Grid_D, double *MOSD_Grid_D, double *MO_VHart_Grid_B, int nI, double **qp, double sdr, double Usic[2])
{
  /****************************************************
     EXC = \sum_{\sigma}
           n_{\sigma}(\epsilon_{xc}-\mu_{xc,\sigma})

     EH = 1/2\int n(r) V_H dr
  ****************************************************/

  static int firsttime=1;
  int i,XC_P_switch,xcs;
  int n,n1,n2,n3,Ng1,Ng2,Ng3,j,k;
  int GN,GNs,BN,DN,LN,N2D,n2D,N3[4];
  double EXC,EH;
  double My_EXC,My_EH;
  double sum,tmp1,tmp2;
  int numprocs,myid,tag=999,ID,IDS,IDR;
  double sden,aden;
  double den[2],E_unif[2],V_unif[2];

  MPI_Status stat;
  MPI_Request request;
  /* for OpenMP */
  int OMPID,Nthrds,Nthrds0,Nprocs,Nloop;

  /* MPI */
  MPI_Comm_size(mpi_comm_level1,&numprocs);
  MPI_Comm_rank(mpi_comm_level1,&myid);

  /****************************************************
   set Vxc_Grid
  ****************************************************/

  XC_P_switch = 5;
  xcs = XC_switch;
  /* if (xcs==4) xcs = 3; */

  Set_XC_Grid(XC_P_switch,xcs,
	      MOD_Grid_D,Null_Grid_D,
	      Null_Grid_D,Null_Grid_D,
	      Vxc_Grid_D[0], Vxc_Grid_D[1],
	      Null_Grid_D, Null_Grid_D,
	      NULL,NULL);

  /****************************************************
             copy Vxc_Grid_D to Vxc_Grid_B
  ****************************************************/

  Ng1 = Max_Grid_Index_D[1] - Min_Grid_Index_D[1] + 1;
  Ng2 = Max_Grid_Index_D[2] - Min_Grid_Index_D[2] + 1;
  Ng3 = Max_Grid_Index_D[3] - Min_Grid_Index_D[3] + 1;

  for (n=0; n<Num_Rcv_Grid_B2D[myid]; n++){
    DN = Index_Rcv_Grid_B2D[myid][n];
    BN = Index_Snd_Grid_B2D[myid][n];

    i = DN/(Ng2*Ng3);
    j = (DN-i*Ng2*Ng3)/Ng3;
    k = DN - i*Ng2*Ng3 - j*Ng3; 

    if ( !(i<=1 || (Ng1-2)<=i || j<=1 || (Ng2-2)<=j || k<=1 || (Ng3-2)<=k)){
      Vxc_Grid_B[0][BN] = Vxc_Grid_D[0][DN];
    }
  }

  /****************************************************
   set RefVxc_Grid from LOSD_Grid 
  ****************************************************/

  if (xcs==2) {
    for (BN=0; BN<My_NumGridB_AB; BN++){
      RefVxc_Grid_B[BN] = XC_Ceperly_Alder(LOSD_Grid_B[BN],XC_P_switch);
    }
  }
  else if (xcs==3 || xcs==4) {

    for (BN=0; BN<My_NumGridB_AB; BN++){
      if (LOSD_Grid_B[BN] > 1e-14) {
	den[0] = LOSD_Grid_B[BN];
	den[1] = 0.0;
	XC_EX(1,2.0*den[0],den,E_unif,V_unif);
	RefVxc_Grid_B[BN] = V_unif[0];
      }
      else RefVxc_Grid_B[BN] = 0.0;
    }
  }

  /****************************************************
        calculations of Eva, Eef, EH1, and EXC
  ****************************************************/

  My_EH = 0.0;
  My_EXC = 0.0;
  sum = 0.0;

  for (BN=0; BN<My_NumGridB_AB; BN++){

    sden = LOD_Grid_B[BN];
    aden = LOSD_Grid_B[BN];

    /* EH = \int dn(r) dV_H dr + 2\int s(r) dV_H dr */
    My_EH += (sden+aden)*MO_VHart_Grid_B[BN];

    /* EXC = \sum_{\sigma} n_{\sigma}\epsilon_{xc}
       - n_{spherical}\epsilon_{xc}(n_{spherical})  */
    My_EXC += sden*Vxc_Grid_B[0][BN] - aden*RefVxc_Grid_B[BN];

    /* sum += aden*RefVxc_Grid_B[BN]; */
  } /* BN */

  /* printf("[%d] TotDen = %lf\n",myid,sum*GridVol); */

  /****************************************************
       multiplying GridVol and adding references
  ****************************************************/

  /* EH = \int dn(r) dV_H dr + 2\int s(r) dV_H dr */
  My_EH *= GridVol;

  /*   EXC = \sum_{\sigma} n_{\sigma}\epsilon_{xc}
       - n_{spherical}\epsilon_{xc}(n_{spherical}) */
  My_EXC *= GridVol;

  /* sum *= GridVol; */

  /****************************************************
   MPI:

   EH, EXC
  ****************************************************/

  MPI_Barrier(mpi_comm_level1);
  MPI_Allreduce(&My_EH, &EH, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);
  MPI_Allreduce(&My_EXC, &EXC, 1, MPI_DOUBLE, MPI_SUM, mpi_comm_level1);

  tmp1 = RefAOUh(qp,sdr);
  tmp2 = RefAOUxc(qp,sdr);

  EH += tmp1;
  EXC += tmp2;

  /* if (myid==Host_ID) printf(" RefUh = %lf (%lf), RefUxc = %lf (%lf)\n",tmp1, RefAOUh(qp,1.0), tmp2, RefAOUxc(qp,1.0)); */

  Usic[0] = EH;
  Usic[1] = EXC;
}

double angle_xc(double denR, double r, double dRdr, double *dFdGr, int i, int pi, int ni, double dr, int xc_switch)
{
  double den;
  double yt2,result;
  double dndr;
  double dndz,E[2],F[3];

  /* xc_switch

     2: LSDA-CA
     3: LSDA-PW
     4: GGA-PBE 1st part
     8: GGA-PBE 2nd part

   */

  switch (xc_switch) {

  case 2:
    yt2 = 0.25/PI;
    den = denR*yt2;
    XC_CA_LSDA(den,0.0,E,1);
    result = E[0];
    break;

  case 3:
    yt2 = 0.25/PI;
    den = denR*yt2;
    result = V_XC_PW_LSDA(den);
    break;

  case 4:
    yt2 = 0.25/PI;
    den = denR*yt2;
    dndr = dRdr*yt2;
    dndz = dndr;
    result = V_XC_GGA(den,0.0,0.0,dndz,F);
    dFdGr[i] = F[2];
    break;

  case 8:
    result = (dFdGr[ni] - dFdGr[pi])/dr + 2.0*dFdGr[i]/r;
    break;

  default:
    result = 1.0;
  }

  return result;
}

double V_XC_GGA(double den0, double dndx, double dndy, double dndz, double F[3])
{
  double den_min=1.0e-14; 
  double Exc[2],Vxc;
  double ED[2],GDENS[3][2];
  double DEXDD[2],DECDD[2];
  double DEXDGD[3][2],DECDGD[3][2];

  if (den0<den_min) {
    F[0] = 0.0;
    F[1] = 0.0;
    F[3] = 0.0;
    return 0.0;
  }

  ED[0] = den0;
  ED[1] = 0.0;

  GDENS[0][0] = dndx;
  GDENS[1][0] = dndy;
  GDENS[2][0] = dndz;
  GDENS[0][1] = 0.0;
  GDENS[1][1] = 0.0;
  GDENS[2][1] = 0.0;

  XC_PBE(ED, GDENS, Exc, DEXDD, DECDD, DEXDGD, DECDGD);

  Vxc = DEXDD[0] + DECDD[0];

  F[0] = DEXDGD[0][0] + DECDGD[0][0];
  F[1] = DEXDGD[1][0] + DECDGD[1][0];
  F[2] = DEXDGD[2][0] + DECDGD[2][0];

  return Vxc;
}

double V_XC_PW_LSDA(double den0)
{
  double Vxc;
  double den[2],E_unif[2],V_unif[2];
  double min_den = 1.0e-14;

  if (den0<min_den) return 0.0;

  den[0] = den0;
  den[1] = 0.0;

  XC_PW92C(den,E_unif,V_unif);
  Vxc = V_unif[0];

  XC_EX(1,2.0*den0,den,E_unif,V_unif);
  Vxc += V_unif[0];

  return Vxc;
}


double RefAOUh(double **qp, double sdr)
{
  int xcs = XC_switch;
  double aoh,a0,b0;
  int spe,GL,Mul,i,j;
  double r,r2,pr,rmin,rmax,h,h2;
  double dRdr;
  double den,fden,bden;

  /* calculation */

  spe = WhatSpecies[target_atom];

  rmin = Spe_PAO_RV[spe][0];
  rmax = Spe_Atom_Cut1[spe] + 0.5;
  h = (rmax - rmin)/(double)OneD_Grid;
  h2 = h*2.0;

  /* Uh: inner integration */

  r = rmin;
  r2 = r*r;
  den = 0.0;
  fden = 0.0;
  for (GL=0; GL<=Spe_MaxL_Basis[spe]; GL++) {
    for (Mul=0; Mul<Spe_Num_Basis[spe][GL]; Mul++){
      den += qp[GL][Mul]*DenRF(spe,GL,Mul,r);
      fden += qp[GL][Mul]*DenRF(spe,GL,Mul,r+h);
    }
  }
  den *= sdr;
  fden *= sdr;
  dRdr = (fden - den)/h;

  a0 = 0.0;
  b0 = 0.5*den*r;

  for (i=1; i<OneD_Grid; i++){
    r = rmin + (double)i*h;
    r2 = r*r;
    bden = den;
    den = fden;
    fden = 0.0;
    for (GL=0; GL<=Spe_MaxL_Basis[spe]; GL++) {
      for (Mul=0; Mul<Spe_Num_Basis[spe][GL]; Mul++){
	fden += qp[GL][Mul]*DenRF(spe,GL,Mul,r+h);
      }
    }
    fden *= sdr;
    dRdr = (fden - bden)/h2;

    b0 += den*r;
  }

  r = rmax;
  r2 = r*r;
  bden = den;
  den = fden;
  dRdr = (den - bden)/h;

  b0 += 0.5*den*r;

  /* Uh: outer integration with reintegration of the inner part */

  r = rmin;
  r2 = r*r;
  den =0.0;
  for (GL=0; GL<=Spe_MaxL_Basis[spe]; GL++) {
    for (Mul=0; Mul<Spe_Num_Basis[spe][GL]; Mul++){
      den += qp[GL][Mul]*DenRF(spe,GL,Mul,r);
    }
  }
  den *= sdr;

  aoh = 0.5*(a0 + r*b0)*den*r;

  for (i=1; i<OneD_Grid; i++){
    pr = r;
    r = rmin + (double)i*h;
    r2 = r*r;
    bden = den;
    den = 0.0;
    for (GL=0; GL<=Spe_MaxL_Basis[spe]; GL++) {
      for (Mul=0; Mul<Spe_Num_Basis[spe][GL]; Mul++){
	den += qp[GL][Mul]*DenRF(spe,GL,Mul,r);
      }
    }
    den *= sdr;

    a0 += (bden*pr*pr + den*r2)/2.0;
    b0 -= (bden*pr + den*r)/2.0;
    aoh += (a0 + r*b0)*den*r;

  }

  pr = r;
  r = rmax;
  r2 = r*r;
  bden = den;
  den = 0.0;
  for (GL=0; GL<=Spe_MaxL_Basis[spe]; GL++) {
    for (Mul=0; Mul<Spe_Num_Basis[spe][GL]; Mul++){
      den += qp[GL][Mul]*DenRF(spe,GL,Mul,r);
    }
  }
  den *= sdr;

  a0 += (bden*pr*pr + den*r2)*0.5;
  b0 -= (bden*pr + den*r)*0.5;

  aoh += 0.5*(a0 + r*b0)*den*r;

  aoh *= h*h;
  return aoh;
}


double RefAOUxc(double **qp, double sdr)
{
  int xcs = XC_switch;
  double aoxc;
  int spe,GL,Mul,i,j;
  double r,r2,pr,rmin,rmax,h,h2;
  double dRdr;
  double den,fden,bden;
  double *dFxcdGr;

  if (xcs==4) xcs==3;

  /* allocation */

  dFxcdGr = (double*)malloc(sizeof(double)*(OneD_Grid+1));

  /* calculation */

  spe = WhatSpecies[target_atom];

  rmin = Spe_PAO_RV[spe][0];
  rmax = Spe_Atom_Cut1[spe] + 0.5;
  h = (rmax - rmin)/(double)OneD_Grid;
  h2 = h*2.0;

  /* Uxc: integration */

  r = rmin;
  r2 = r*r;
  den = 0.0;
  fden = 0.0;
  for (GL=0; GL<=Spe_MaxL_Basis[spe]; GL++) {
    for (Mul=0; Mul<Spe_Num_Basis[spe][GL]; Mul++){
      den += qp[GL][Mul]*DenRF(spe,GL,Mul,r);
      fden += qp[GL][Mul]*DenRF(spe,GL,Mul,r+h);
    }
  }
  den *= sdr;
  fden *= sdr;

  aoxc = 0.5*angle_xc(den, r, dRdr, dFxcdGr, 0, 0, 1, h, xcs)*den*r2;

  for (i=1; i<OneD_Grid; i++){
    r = rmin + (double)i*h;
    r2 = r*r;
    den = fden;
    fden = 0.0;
    for (GL=0; GL<=Spe_MaxL_Basis[spe]; GL++) {
      for (Mul=0; Mul<Spe_Num_Basis[spe][GL]; Mul++){
	fden += qp[GL][Mul]*DenRF(spe,GL,Mul,r+h);
      }
    }
    fden *= sdr;

    aoxc += angle_xc(den, r, dRdr, dFxcdGr, i, i-1, i+1, h2, xcs)*den*r2;
  }

  r = rmax;
  r2 = r*r;
  bden = den;
  den = fden;

  aoxc += 0.5*angle_xc(den, r, dRdr, dFxcdGr, i, i-1, i, h, xcs)*den*r;

  aoxc *= h;

  /* freeing */

  free(dFxcdGr);

  return aoxc;

}
