/*
 *===============================================
 *  Created by Andreas Pedersen on 10/17/06.
 *-----------------------------------------------
 *  Modified. Name, Date and a small description!
 *
 *-----------------------------------------------
 *  Todo:
 *
 *===============================================
 */
#include <cassert>
#include "HelperFunctions.h"

// Random values in interval
double helper_functions::randomDouble(){
    // RAND_MAX is in the library math.h
    return((double)rand()/((double)RAND_MAX+1));
}

double helper_functions::randomDouble(int max){
    double dmax;
    dmax = double(max);
    // RAND_MAX is in the library math.h
    return(dmax*(double)rand()/((double)RAND_MAX+1));
}

double helper_functions::randomDouble(long max){
    double dmax;
    dmax = double(max);
    // RAND_MAX is in the library math.h
    return(dmax*(double)rand()/((double)RAND_MAX+1));
}

double helper_functions::randomDouble(double dmax){
    // RAND_MAX is in the library math.h
    return(dmax*(double)rand()/((double)RAND_MAX+1));
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
        assert(isfinite(v1[i]));
        assert(isfinite(v2[i]));
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
    assert(isfinite(scalar));
    for(int i=0;i<size;i++) {
        assert(isfinite(v1[i]));
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
    assert(norm != 0.0);
    divideScalar(v1, v1, norm, size);
    return;
}

// Make v1 orthogonal to v2 
void helper_functions::makeOrthogonal(double *result, const double *v1, 
                                      const double *v2, long size){
    double *tempListDouble;
    double tempDouble;
    tempListDouble = new double[size];
    
    tempDouble = dot(v2, v1, size);
    multiplyScalar(tempListDouble, v2, tempDouble, size);
    subtract(result, v1, tempListDouble, size);
    
    delete [] tempListDouble;
    return;
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
