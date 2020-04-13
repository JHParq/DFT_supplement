/**********************************************************************
  DenRF.c:

     RadialF.c is a subroutine to calculate the radial function of
     pseudo atomic orbital density specified by "l" for atomic species
     at R.

  Log of DenRF.c:

     Mar/2020  Released by J.Parq

***********************************************************************/

#include <stdio.h>
#include <math.h>
#include "openmx_common.h"

double DenRF(int spe, int L, int Mul, double R)
{
  /* radial function of atomic density */

  int mp_min,mp_max,m;
  double h1,h2,h3,f1,f2,f3,f4;
  double g1,g2,x1,x2,y1,y2,f;
  double result;
  double tmp;

  mp_min = 0;
  mp_max = Spe_Num_Mesh_PAO[spe] - 1;

  if (R<Spe_PAO_RV[spe][0]){
    if (L==0) {
      tmp = Spe_PAO_RWF[spe][L][Mul][0];
      f = tmp*tmp;
    }
    else      f = 0.0;
  }
  else if (Spe_PAO_RV[spe][Spe_Num_Mesh_PAO[spe]-1]<R){
    f = 0.0;
  }
  else{
    do{
      m = (mp_min + mp_max)/2;
      if (Spe_PAO_RV[spe][m]<R)
        mp_min = m;
      else 
        mp_max = m;
    }
    while((mp_max-mp_min)!=1);
    m = mp_max;

    if (m<2)
      m = 2;
    else if (Spe_Num_Mesh_PAO[spe]<=m)
      m = Spe_Num_Mesh_PAO[spe] - 2;

    /****************************************************
                   Spline like interpolation
    ****************************************************/

    if (m==1){
      h2 = Spe_PAO_RV[spe][m]   - Spe_PAO_RV[spe][m-1];
      h3 = Spe_PAO_RV[spe][m+1] - Spe_PAO_RV[spe][m];

      tmp = Spe_PAO_RWF[spe][L][Mul][m-1];
      f2 = tmp*tmp;
      tmp = Spe_PAO_RWF[spe][L][Mul][m];
      f3 = tmp*tmp;
      tmp = Spe_PAO_RWF[spe][L][Mul][m+1];
      f4 = tmp*tmp;

      h1 = -(h2+h3);
      f1 = f4;
    }
    else if (m==(Spe_Num_Mesh_PAO[spe]-1)){
      h1 = Spe_PAO_RV[spe][m-1] - Spe_PAO_RV[spe][m-2];
      h2 = Spe_PAO_RV[spe][m]   - Spe_PAO_RV[spe][m-1];

      tmp = Spe_PAO_RWF[spe][L][Mul][m-2];
      f1 = tmp*tmp;
      tmp = Spe_PAO_RWF[spe][L][Mul][m-1];
      f2 = tmp*tmp;
      tmp = Spe_PAO_RWF[spe][L][Mul][m];
      f3 = tmp*tmp;

      h3 = -(h1+h2);
      f4 = f1;
    }
    else{
      h1 = Spe_PAO_RV[spe][m-1] - Spe_PAO_RV[spe][m-2];
      h2 = Spe_PAO_RV[spe][m]   - Spe_PAO_RV[spe][m-1];
      h3 = Spe_PAO_RV[spe][m+1] - Spe_PAO_RV[spe][m];

      tmp = Spe_PAO_RWF[spe][L][Mul][m-2];
      f1 = tmp*tmp;
      tmp = Spe_PAO_RWF[spe][L][Mul][m-1];
      f2 = tmp*tmp;
      tmp = Spe_PAO_RWF[spe][L][Mul][m];
      f3 = tmp*tmp;
      tmp = Spe_PAO_RWF[spe][L][Mul][m+1];
      f4 = tmp*tmp;
    }

    /****************************************************
                calculate the value at R
    ****************************************************/

    g1 = ((f3-f2)*h1/h2 + (f2-f1)*h2/h1)/(h1+h2);
    g2 = ((f4-f3)*h2/h3 + (f3-f2)*h3/h2)/(h2+h3);

    x1 = R - Spe_PAO_RV[spe][m-1];
    x2 = R - Spe_PAO_RV[spe][m];
    y1 = x1/h2;
    y2 = x2/h2;

    f =  y2*y2*(3.0*f2 + h2*g1 + (2.0*f2 + h2*g1)*y2)
       + y1*y1*(3.0*f3 - h2*g2 - (2.0*f3 - h2*g2)*y1);
  }
  result = f;
  return result;
}
