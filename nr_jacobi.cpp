// -*- c++ -*-

// #include "colvartypes.h"
#include "nr_jacobi.h"
#include <cmath>


#define ROTATE(a,i,j,k,l) g=a[i][j]; \
  h=a[k][l];                         \
  a[i][j]=g-s*(h+g*tau);             \
  a[k][l]=h+s*(g-h*tau);

#define n 4


namespace NR_Jacobi {

int jacobi(double a[4][4], double d[4], double v[4][4], int *nrot)
{
  int j,iq,ip,i;
  double tresh,theta,tau,t,sm,s,h,g,c;

  double b[n];
  double z[n];

  for (ip=0;ip<n;ip++) {
    for (iq=0;iq<n;iq++) {
      v[ip][iq]=0.0;
    }
    v[ip][ip]=1.0;
  }
  for (ip=0;ip<n;ip++) {
    b[ip]=d[ip]=a[ip][ip];
    z[ip]=0.0;
  }
  *nrot=0;
  for (i=0;i<=50;i++) {
    sm=0.0;
    for (ip=0;ip<n-1;ip++) {
      for (iq=ip+1;iq<n;iq++)
        sm += std::fabs(a[ip][iq]);
    }
    if (sm == 0.0) {
      return COLVARS_OK;
    }
    if (i < 4)
      tresh=0.2*sm/(n*n);
    else
      tresh=0.0;
    for (ip=0;ip<n-1;ip++) {
      for (iq=ip+1;iq<n;iq++) {
        g=100.0*std::fabs(a[ip][iq]);
        if (i > 4 && (double)(std::fabs(d[ip])+g) == (double)std::fabs(d[ip])
            && (double)(std::fabs(d[iq])+g) == (double)std::fabs(d[iq]))
          a[ip][iq]=0.0;
        else if (std::fabs(a[ip][iq]) > tresh) {
          h=d[iq]-d[ip];
          if ((double)(std::fabs(h)+g) == (double)std::fabs(h))
            t=(a[ip][iq])/h;
          else {
            theta=0.5*h/(a[ip][iq]);
            t=1.0/(std::fabs(theta)+std::sqrt(1.0+theta*theta));
            if (theta < 0.0) t = -t;
          }
          c=1.0/std::sqrt(1+t*t);
          s=t*c;
          tau=s/(1.0+c);
          h=t*a[ip][iq];
          z[ip] -= h;
          z[iq] += h;
          d[ip] -= h;
          d[iq] += h;
          a[ip][iq]=0.0;
          for (j=0;j<=ip-1;j++) {
            ROTATE(a,j,ip,j,iq)
              }
          for (j=ip+1;j<=iq-1;j++) {
            ROTATE(a,ip,j,j,iq)
              }
          for (j=iq+1;j<n;j++) {
            ROTATE(a,ip,j,iq,j)
              }
          for (j=0;j<n;j++) {
            ROTATE(v,j,ip,j,iq)
              }
          ++(*nrot);
        }
      }
    }
    for (ip=0;ip<n;ip++) {
      b[ip] += z[ip];
      d[ip]=b[ip];
      z[ip]=0.0;
    }
  }
  return COLVARS_ERROR;
}


int eigsrt(double d[4], double v[4][4])
{
  int k,j,i;
  double p;

  for (i=0;i<n;i++) {
    p=d[k=i];
    for (j=i+1;j<n;j++)
      if (d[j] >= p) p=d[k=j];
    if (k != i) {
      d[k]=d[i];
      d[i]=p;
      for (j=0;j<n;j++) {
        p=v[j][i];
        v[j][i]=v[j][k];
        v[j][k]=p;
      }
    }
  }
  return COLVARS_OK;
}


int transpose(double v[4][4])
{
  double p;
  int i,j;
  for (i=0;i<n;i++) {
    for (j=i+1;j<n;j++) {
      p=v[i][j];
      v[i][j]=v[j][i];
      v[j][i]=p;
    }
  }
  return COLVARS_OK;
}

}

#undef n
#undef ROTATE
