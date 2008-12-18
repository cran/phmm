#include <math.h>
#define NRANSI
#include "nrutil.h"
#include "phmm.h"

#define FREERETURN {return;}

/*#define FREERETURN {free_matrix(fjac,1,n,1,n);free_vector(fvec,1,n);\
    free_vector(p,1,n);free_ivector(indx,1,n);return;}*/

void mnewt(int n, int ntrial, float *x, float tolx, float tolf,
   int nobs, double *omega, double *sum0, double **sum1, double **z,
   int *delta, double *sum2[n+1][n+1])
{
    void lubksb(float **a, int n, int *indx, float b[]);
    void ludcmp(float **a, int n, int *indx, float *d);
    int k,i,*indx;
    float errx,errf,d,*fvec,**fjac,*p;

    indx=ivector(1,n);
    p=vector(1,n);
    fvec=vector(1,n);
    fjac=matrix(1,n,1,n);
    errx=tolx+1;
    for (k=1;k<=ntrial;k++) {            /*printf("trial no.1");*/
            usrfun(x,n,fvec,fjac,nobs,omega,sum0,sum1,z,delta,sum2);
        if (errx <= tolx) FREERETURN
        errf=0.0;
        for (i=1;i<=n;i++) errf += fabs(fvec[i]);
        if (errf <= tolf) FREERETURN
        for (i=1;i<=n;i++) p[i] = -fvec[i];
        ludcmp(fjac,n,indx,&d);
        lubksb(fjac,n,indx,p);
        errx=0.0;
        for (i=1;i<=n;i++) {
            errx += fabs(p[i]);
            x[i] += p[i];
        }
    }
    FREERETURN
}
#undef FREERETURN
#undef NRANSI

/* BETTER MOVED IN THE SAME FILE WITH   mnewt() */
void usrfun(float *beta, int n, float *score, float **fjac,
   int nobs, double *omega, double *sum0, double *sum1[n], double **z,
   int *delta, double *sum2[n+1][n+1])
{           /* used by mnewt(); updated for frailties */

   int i, j, k;
   double temp;

/** printf("usrfun "); **/
   for (j=1; j<=n; j++)   {
      score[j]=0;
      for (k=1; k<=n; k++)   fjac[j][k]=0;
   }

   for ( i=nobs; i>=1; i-- )  {
      temp=exp(BetaZ(i, beta, n, z))*omega[i];
      sum0[i] = sum0[i+1]+temp;
      for (j=1;j<=n;j++) {
         sum1[j][i] = sum1[j][i+1]+temp*z[j][i];
         if (delta[i])   score[j]+= (z[j][i]-sum1[j][i]/sum0[i]);
      }
      for (j=1;j<=n;j++)
         for (k=1;k<=n;k++)  {
            sum2[j][k][i] = sum2[j][k][i+1]+temp*z[j][i]*z[k][i];
            if (delta[i])
              fjac[j][k]+= ((sum1[j][i]/sum0[i])*(sum1[k][i]/sum0[i]) -
            sum2[j][k][i]/sum0[i]);
         }
   }
}   /* end usrfun() */
