/*
 *===============================================
 *  EON HelperFunctions.cpp
 *===============================================
 */

#include "HelperFunctions.h"

#include <math.h>
#include <cassert>

// Random number generator

double helper_functions::random(long newSeed){
    static long seed = -1;
    if(newSeed){
        seed=-newSeed;}
    int j;
    long k;
    static long seed2=123456789;
    static long iy=0;
    static long iv[NTAB];
    double temp;
    if (seed <= 0) {
        if (-(seed) < 1) seed=3;
        else seed = -(seed);
        seed2=(seed);
        for (j=NTAB+7;j>=0;j--) {
            k=(seed)/IQ1;
            seed=IA1*(seed-k*IQ1)-k*IR1;
            if (seed < 0) seed += IM1;
            if (j < NTAB) iv[j] = seed;}
        iy=iv[0];}
    k=(seed)/IQ1;
    seed=IA1*(seed-k*IQ1)-k*IR1;
    if (seed < 0) seed += IM1;
    k=seed2/IQ2;
    seed2=IA2*(seed2-k*IQ2)-k*IR2;
    if (seed2 < 0) seed2 += IM2;
    j=int(iy/NDIV);
    iy=iv[j]-seed2;
    iv[j] = seed;
    if (iy < 1) iy += IMM1;
    if ((temp=double(AM*iy)) > RNMX) return RNMX;
    else return temp;
}

double helper_functions::randomDouble(){
    return(random());
}


// Random value in interval
double helper_functions::randomDouble(int max){
    double dmax = double(max);
    return(dmax*randomDouble());
}

double helper_functions::randomDouble(long max){
    double dmax = double(max);
    return(dmax*randomDouble());
}

double helper_functions::randomDouble(double dmax){
    return(dmax*randomDouble());
}

double helper_functions::guaRandom(double avg,double std){
	double r=2,v1,v2,l,result;
    while (r >= 1){
	    v1=2.0*randomDouble()-1.0;
		v2=2.0*randomDouble()-1.0;
	    r=v1*v1+v2*v2;
	}
    l=v1*sqrt(-2.0*log(r)/r);
	result=avg+std*l;
	return(result);
}

// Vector functions.
// All functions 'returning' an array 
// the first argument should be a pointer to the array 
// where the result should be stored.
double helper_functions::dot(const double *v1, const double *v2, long size){
    double result=0;
    for(int i=0;i<size;i++)
        result = result+v1[i]*v2[i];
    return result;
}

double helper_functions::length(const double *v1, long numFreeCoord){
    return(sqrt(dot(v1, v1, numFreeCoord)));
}

void helper_functions::add(double *result, const double *v1, 
                           const double *v2, long size){
    for(int i=0;i<size;i++) {
        result[i] = v1[i]+v2[i];
    };
    return;
}

void helper_functions::subtract(double *result, const double *v1, 
                                const double * v2, long size){
    for(int i=0;i<size;i++)
        result[i] = v1[i]-v2[i];
    return;
}

void helper_functions::multiplyScalar(double *result, const double *v1, 
                                      double scalar, long size) {
    for(int i=0;i<size;i++) {
        result[i] = v1[i]*scalar;
    };
    return;
}

void helper_functions::divideScalar(double *result, const double *v1, 
                                    double scalar, long size){
    for(int i=0;i<size;i++)
        result[i] = v1[i]/scalar;
    return;
}

void helper_functions::copyRightIntoLeft(double *result, const double *v1, 
                                         long size){
    for(int i=0;i<size;i++)
        result[i] = v1[i];
    return;
}

void helper_functions::normalize(double *v1, long size){
    double const norm=length(v1, size);
    //XXX: dirty hack while waiting for merge
    if(norm == 0.0){throw 14323;}
    divideScalar(v1, v1, norm, size);
    return;
}

// Make v1 orthogonal to v2 
Matrix<double, Eigen::Dynamic, 3> helper_functions::makeOrthogonal(const Matrix<double, Eigen::Dynamic, 3> v1, const Matrix<double, Eigen::Dynamic, 3> v2){
    return v1 - (v1.cwise()*v2).sum() * v2.normalized();
}

// result contains v1 projection on v2 
void helper_functions::makeProjection(double *result, const double *v1, 
                                      const double *v2, long size){
    double *tempListDouble;
    double tempDouble;
    tempListDouble = new double[size];
    
    tempDouble = dot(v1, v2, size);
    multiplyScalar(result, v2, tempDouble, size);
    
    delete [] tempListDouble;
    return;
}
