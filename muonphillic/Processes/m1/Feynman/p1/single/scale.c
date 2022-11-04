#include<math.h>
#include<stdlib.h>
#include"/home/duncan/UChicago/MuonBeam/calchep_3.8.10/include/nType.h"
#define min(x,y) (x<y? x:y)
#define max(x,y) (x>y? x:y)
extern void ScaleCC(int nsub, REAL*, double (*calcPV)(char,char*,double*), double*,double*,double*,double*,double *);
void ScaleCC(int nsub,REAL*modelVal, double (*calcPV)(char,char*,double*), double *pvect,double *Qren, double *Qpdf1,double *Qpdf2, double*Qshow)
{ double X[9];          
  double ss,s1; int i;
  X[0]=calcPV('M',"\1\2",pvect);
X[1]=X[0];
 *Qren= X[1];
  X[2]=calcPV('M',"\1\2",pvect);
X[3]=X[2];
 *Qpdf1= X[3];
X[5]=X[3];
 *Qpdf2= X[5];
X[7]=X[1];
 *Qshow= X[7];
}
