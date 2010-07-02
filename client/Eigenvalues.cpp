/*
 *===============================================
 *  Taken from Numerical Recipies
 *-----------------------------------------------
 *  Modified. Name, Date and a small description!
 *
 *-----------------------------------------------
 *  Todo:
 *
 *===============================================
 */
#include <iostream>
#include <cmath>
#include <cassert>
#include "Eigenvalues.h"

void eigenValues(const int n,double *d, double **H){
      assert(n>0);
      double *e;
      e = new double[n];
      MakeSymmetric(n,H);
      TriDiagonalize(n,d,e,H);
      Diagonalize(n,d,e);
      delete [] e;
}

void MakeSymmetric(const int n, double **H){
	for (int i=0; i<n ;i++){
	  for (int j=i+1; j<n; j++){
	    H[i][j]=H[j][i]=0.5*(H[i][j]+H[j][i]);}}
}

void TriDiagonalize(const int n,double *d,double *e,double **H)
{
      assert(n>0);
      int i,j,k,l;
      double f,g,h,hh,sigma;
      for (i=n-1; i>0; i--){
            l=i-1;
            h=sigma=double(0.0);
            if(l > 0){
                  for(k=0; k<=l; k++){
                        sigma+=fabs(H[i][k]);
                  }
                  if(sigma == double(0.0)){
                        e[i]=H[i][l];
                  }
                  else{
                        for(k=0; k<=l; k++){
                              H[i][k]=H[i][k]/sigma;
                              h+=H[i][k]*H[i][k];
                        }
                        f=H[i][l];
                        g=-Sign(f)*sqrt(h);
                        e[i]=sigma*g;
                        h-=f*g;
                        H[i][l]=f-g;
                        f=double(0.0);
                        for(j=0; j<=l; j++){
                              g=double(0.0);
                              for(k=0; k<=j; k++){
                                    g+=H[j][k]*H[i][k];
                              }
                              for(k=j+1; k<=l; k++){
                                    g+=H[k][j]*H[i][k];
                              }
                              e[j]=g/h;
                              f+=e[j]*H[i][j];}
                        hh=0.5*f/h;
                        for(j=0; j<=l; j++){
                              f=H[i][j];
                              g=e[j]-hh*f;
                              e[j]=g;
                              for(k=0; k<=j; k++){
                                    H[j][k]-=f*e[k]+g*H[i][k];
                              }
                        }
                  }
            }
            else{
                  e[i]=H[i][l];
            }
      }
      e[0]=double(0.0);
      for(i=0; i<n; i++){
            d[i]=H[i][i];}
}

void Diagonalize(const int n,double *d,double *e){
	int i,itr,l,m;
	double b,c,dd,f,g,p,r,s;
	for(i=1; i<n; i++){
	  e[i-1]=e[i];}
	e[n-1]=double(0.0);
	for(l=0; l<n; l++){
	  itr=0;
	  do{
	    for(m=l; m<n-1; m++){
	      dd=Abs(d[m])+Abs(d[m+1]);
	      if( (double)(Abs(e[m])+dd) == dd){
	        break;}}
	    if(m != l){
	      if(itr++ == 30) std::cerr << "Too many iteractions in Diagonalize\n";
	      g=0.5*(d[l+1]-d[l])/e[l];
	      r=scaling(g,double(1.0));
	      g=d[m]-d[l]+e[l]/(g+Sign(g)*r);
	      s=c=double(1.0);
	      p=double(0.0);
	      for(i=m-1;i>=l;i--){
	        f=s*e[i];
	        b=c*e[i];
	        e[i+1]=r=scaling(f,g);
	        if(r == double(0.0)){
	          d[i+1]-=p;
	          e[m]=double(0.0);
	          break;}
	        s=f/r;
	        c=g/r;
	        g=d[i+1]-p;
	        r=(d[i]-g)*s+2.0*c*b;
	        p=s*r;
	        d[i+1]=g+p;
	        g=c*r-b;}
	      if(r == double(0.0) && i >= l) continue;
	      d[l]-=p;
	      e[l]=g;
	      e[m]=double(0.0);}
	}while(m != l);}
}

double scaling(double a,double b){
	double absa,absb;
	absa=Abs(a);
	absb=Abs(b);
	if(absa > absb) return absa*sqrt(1.0+(absb/absa)*(absb/absa));
	else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+(absa/absb)*(absa/absb)));
}

