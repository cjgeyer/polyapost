/*Computes the Markov chain of the points in the polytope*/
/*for combined equality and inequality constraints A1x=b1,A2x<=b2. It returns the last value.*/

#include <R.h>
double sum1( double x[], int vectsize )
{
     double s;
     int i;
     s=0.0;
     for (i=0; i<vectsize; i++)
          s=s+x[i]*x[i];
     return(s);
}

void f1(double x[], double y[], int vectSize, double* pMin, double* pMax)
{
  int i;
  double maxDiv=-DBL_MAX; /*I make sure that,initially, any value is bigger than maxDiv*/
  double minDiv=DBL_MAX; /*I make sure that,initially, any value is smaller than minDiv*/
  double div;

  for(i=0; i<vectSize; i++)
  {
    if(fabs(y[i]) > 0.00000001) /* don't do anything if y[i]==0 */
    {
      div = -x[i]/y[i];

      if(y[i] < 0)
      {
        /* update max */
        maxDiv = (div > maxDiv)? div : maxDiv;
      }
      else
      {
        /* update min */
        minDiv = (div < minDiv)? div : minDiv;
      }
    }
  }

  /* return min and max calculated */

  *pMin = minDiv;
  *pMax = maxDiv;
}





/*p is B%*%t(B) where B is the matrix that generates N(A1)
matsize is the dimension of p
a2 constraint matrix for the inequality part of dimension nrow*matsize
b2 the right hand side vector of the inequality part of dimension nrow
result stores the means we generate
*/

#define P(i,j) p[(i)+n*(j)]
#define A2(i,j) a2[(i)+m*(j)]
void gen1( double *p, int matsize, double *a2, int nrow, double *b2, double *initsol,
 double *result,  double *z, double *z1, double *d, double *d1, double *d2,
          double *x1, double* x2)
{
  int m=nrow;
  int n=matsize;

  double s;
  double s1;
  double lambda;
  double lambda1;
  double lambda2;
  int i,j;

  /*generate  nornmal variables*/
  for(i=0;i<n;i++)
  {
    z1[i]=norm_rand();
  }
  /*generate uniform variable on the sphere*/
  s=sqrt(sum1(z1,n));
  for(i=0;i<n;i++)
  {
    z[i]=z1[i]/s;
  }
  /*projection of the uniform variable onto the null space*/
  for(i=0;i<n;i++)
  {
    d1[i]=0.0;
    for(j=0;j<n;j++)
    {
      d1[i]+=P(i,j)*z[j];
    }
  }
  /*generate uniform variable on the nullspace and sphere in a lower dim*/
  s1=sqrt(sum1(d1,n));
  for(i=0;i<n;i++)
  {
    d[i]=d1[i]/s1;
  }
  /*multiply A2 and initsol*/
  for(i=0;i<m;i++)
  {
    x1[i]=0.0;
    for(j=0;j<n;j++)
    {
      x1[i]+=A2(i,j)*initsol[j];
    }
  }
  /*multiply A2 and d*/
  for(i=0;i<m;i++)
  {
    d2[i]=0.0;
    for(j=0;j<n;j++)
    {
      d2[i]+=A2(i,j)*d[j];
    }
  }
  /*generate x0-b2*/
  for(i=0;i<m;i++)
  {
    x2[i]=x1[i]-b2[i];
  }
  /*get the bounds in the direction d*/
  f1(x2,d2,m,&lambda2,&lambda1);
  lambda=lambda1+(lambda2-lambda1)*unif_rand();
  /*generate a  new feasible solution*/
  for(i=0;i<n;i++)
  {
    result[i]=initsol[i]+lambda*d[i];
  }

}

/*generates the Markov chain and returns the last point
*/
void probvect1( double *p, int *matsize, double *a2, int *nrow, double *b2,
double *initsol, int * length, double *estimate)
{
  int n=matsize[0];
  int m=nrow[0];
  int k=length[0];
  int i,j;
  double *result = (double *) R_alloc(n, sizeof(double));

  double *z = (double *) R_alloc(n, sizeof(double));
  double *z1 = (double *) R_alloc(n, sizeof(double));
  double *d = (double *) R_alloc(n, sizeof(double));
  double *d1 = (double *) R_alloc(n, sizeof(double));
  double *d2 = (double *) R_alloc(m, sizeof(double));
  double *x1 = (double *) R_alloc(m, sizeof(double));
  double *x2 = (double *) R_alloc(m, sizeof(double));
   GetRNGstate(); 

  for(i=0;i<k;i++)
  {
     gen1(p,n,a2,m,b2,initsol,result,z,z1,d,d1,d2,x1,x2);
     for(j=0;j<n;j++)
     {
      initsol[j]=result[j];
     }
  }
  for(j=0;j<n;j++)
     {
       estimate[j]=initsol[j];
     }
  PutRNGstate();
}


