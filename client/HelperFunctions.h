/*
 *===============================================
 *  HelperFunctions.h
 *  eon2
 *-----------------------------------------------
 *  Created by Andreas Pedersen on 10/17/06.
 *-----------------------------------------------
 *  Modified. Name, Date and a small description!
 *
 *-----------------------------------------------
 *  Todo:
 *
 *===============================================
 */
#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H

#include <stdlib.h>
#include <math.h>

// remember to delete
#include <iostream>


/** Collection of supporting functions. Functions that handle arrays of doubles as vectors and different random number generators*/
namespace helper_functions {

    double randomDouble();///< Random value between 0 and 1.
    double randomDouble(int max);///< Random value between 0 and \a max.
    double randomDouble(long max);///< Random value between 0 and \a max.
    double randomDouble(double max);///< Random value between 0 and \a max.

    /** Calculate the dot product of the two (\a v1 \a v2) double arrays being passed in of length \a size.
    @param[in]   *v1    Array of double.
    @param[in]   *v2    Array of double.
    @param[in]   size   Number of elements in one of the arrays.*/
    double dot(const double *v1, const double *v2, long size);

    /** Calculate the length of (\a v1) array being passed in of length \a size.
    @param[in]   *v1    Array of double.
    @param[in]   size   Number of elements in array.*/
    double length(const double *v1, long size);

    /** Calculate the sum of the two (\a v1 \a v2) double arrays being passed in of length \a size in \a result.
    @param[out]  *result  Array of double.
    @param[in]   *v1      Array of double.
    @param[in]   *v2      Array of double.
    @param[in]   size     Number of elements in one of the arrays.*/
    void add(double *result, const double *v1, const double *v2, long size);

    /** Calculate the difference of the two (\a v1 \a v2) double arrays being passed in of length \a size in \a result.
    @param[out]  *result  Array of double.
    @param[in]   *v1      Array of double.
    @param[in]   *v2      Array of double.
    @param[in]   size     Number of elements in one of the arrays.*/
    void subtract(double *result, const double *v1, const double *v2, long size);

    /** Calculate the product of the (\a v1) double array being passed in of length \a size and scalar in \a result.
    @param[out]  *result  Array of double.
    @param[in]   *v1      Array of double.
    @param[in]   scalar   Value of double.
    @param[in]   size     Number of elements in one of the arrays.*/
    void multiplyScalar(double *result, const double *v1, double scalar, long size);

    /** Calculate the ratio of the (\a v1) double array being passed in of length \a size and scalar in \a result.
    @param[out]  *result  Array of double.
    @param[in]   *v1      Array of double.
    @param[in]   scalar   Value of double.
    @param[in]   size     Number of elements in one of the arrays.*/
    void divideScalar(double *result, const double *v1, double scalar, long size);

    /** Copy \a v2 into \a v1, double arrays of length \a size.
    @param[out]  *result  Array of double.
    @param[in]   *v1      Array of double.
    @param[in]   size     Number of elements in one of the arrays.*/
    void copyRightIntoLeft(double *result, const double *v1, long size);

    /** Normalize the double array \a v1 of length \a size.
    @param[in]   *v1      Array of double.
    @param[out]  *v1      Array of double.
    @param[in]   size     Number of elements in one of the arrays.*/
    void normalize(double *v1, long size);

    /** Compute the orthogonal part of \a v1 to \a v2, of length \a size and store in \a result.
    @param[out]  *result  Array of double.
    @param[in]   *v1      Array of double.
    @param[in]   *v2      Array of double.
    @param[in]   size     Number of elements in one of the arrays.*/
    void makeOrthogonal(double *result, const double *v1, const double *v2, long size);

    /** Compute the projection of \a v1 on \a v2, of length \a size and store in \a result.
    @param[out]  *result  Array of double.
    @param[in]   *v1      Array of double.
    @param[in]   *v2      Array of double.
    @param[in]   size     Number of elements in one of the arrays.*/
    void makeProjection(double *result, const double *v1, const double *v2, long size);
}
#endif
