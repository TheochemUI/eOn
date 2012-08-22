//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "HelperFunctions.h"
#include "Log.h"

#include <math.h>
#include <cassert>
#include <iostream>
#include <time.h>
#include <string.h>

#ifndef WIN32
    #include <sys/time.h>
    #include <sys/resource.h>
#endif

// Random number generator

double helper_functions::random(long newSeed)
{
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

long helper_functions::randomInt(int lower, int upper) {
    return lround( (upper-lower)*randomDouble()+lower );
}

double helper_functions::gaussRandom(double avg,double std){
    double r=2,v1,v2,l,result;
    while (r >= 1){
        v1 = 2.0*randomDouble()-1.0;
        v2 = 2.0*randomDouble()-1.0;
        r = v1*v1+v2*v2;
    }
    l = v1*sqrt(-2.0*log(r)/r);
    result = avg+std*l;
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

void helper_functions::add(double *result, const double *v1, const double *v2, long size){
    for(int i=0;i<size;i++) {
        result[i] = v1[i]+v2[i];
    };
    return;
}

void helper_functions::subtract(double *result, const double *v1, const double * v2, long size){
    for(int i=0;i<size;i++)
        result[i] = v1[i]-v2[i];
    return;
}

void helper_functions::multiplyScalar(double *result, const double *v1, double scalar, long size){
    for(int i=0;i<size;i++) {
        result[i] = v1[i]*scalar;
    };
    return;
}

void helper_functions::divideScalar(double *result, const double *v1, double scalar, long size){
    for(int i=0;i<size;i++)
        result[i] = v1[i]/scalar;
    return;
}

void helper_functions::copyRightIntoLeft(double *result, const double *v1, long size){
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
AtomMatrix helper_functions::makeOrthogonal(const AtomMatrix v1, const AtomMatrix v2){
    return v1 - (v1.cwise()*v2).sum() * v2.normalized();
}

// result contains v1 projection on v2
void helper_functions::makeProjection(double *result, const double *v1, const double *v2, long size){
    double *tempListDouble;
    double tempDouble;
    tempListDouble = new double[size];

    tempDouble = dot(v1, v2, size);
    multiplyScalar(result, v2, tempDouble, size);

    delete [] tempListDouble;
    return;
}

bool helper_functions::rotationMatch(const Matter *m1, const Matter *m2, const double max_diff)
{
    AtomMatrix r1 = m1->getPositions();
    AtomMatrix r2 = m2->getPositions();
    
    // Align centroids
    Eigen::VectorXd c1(3);
    Eigen::VectorXd c2(3);

    c1[0] = r1.col(0).sum();
    c1[1] = r1.col(1).sum();
    c1[2] = r1.col(2).sum();
    c2[0] = r2.col(0).sum();
    c2[1] = r2.col(1).sum();
    c2[2] = r2.col(2).sum();
    c1/=r1.rows();
    c2/=r2.rows();

    for(int i = 0; i < r1.rows(); i++)
    {
        r1(i,0) -= c1[0];
        r1(i,1) -= c1[1];
        r1(i,2) -= c1[2];

        r2(i,0) -= c2[0];
        r2(i,1) -= c2[1];
        r2(i,2) -= c2[2];
    }

    //Determine optimal rotation
    //Horn, J. Opt. Soc. Am. A, 1987
    Eigen::Matrix3d m = r1.transpose() * r2;

    double sxx = m(0,0);
    double sxy = m(0,1);
    double sxz = m(0,2);
    double syx = m(1,0);
    double syy = m(1,1);
    double syz = m(1,2);
    double szx = m(2,0);
    double szy = m(2,1);
    double szz = m(2,2);

    Eigen::Matrix4d n;
    n.setZero();
    n(0,1) = syz-szy;
    n(0,2) = szx-sxz;
    n(0,3) = sxy-syx;

    n(1,2) = sxy+syx;
    n(1,3) = szx+sxz;

    n(2,3) = syz+szy;

    n += n.transpose().eval();

    n(0,0) = sxx + syy + szz;
    n(1,1) = sxx-syy-szz;
    n(2,2) = -sxx + syy -szz;
    n(3,3) = -sxx -syy + szz;

    Eigen::SelfAdjointEigenSolver<Matrix4d> es(n);
    Eigen::Vector4d maxv = es.eigenvectors().col(3);

    Eigen::Matrix3d R;

    double aa = maxv[0]*maxv[0];
    double bb = maxv[1]*maxv[1];
    double cc = maxv[2]*maxv[2];
    double dd = maxv[3]*maxv[3];
    double ab = maxv[0]*maxv[1];
    double ac = maxv[0]*maxv[2];
    double ad = maxv[0]*maxv[3];
    double bc = maxv[1]*maxv[2];
    double bd = maxv[1]*maxv[3];
    double cd = maxv[2]*maxv[3];

    R(0,0) = aa + bb - cc - dd;
    R(0,1) = 2*(bc-ad);
    R(0,2) = 2*(bd+ac);
    R(1,0) = 2*(bc+ad);
    R(1,1) = aa - bb + cc - dd;
    R(1,2) = 2*(cd-ab); 
    R(2,0) = 2*(bd-ac); 
    R(2,1) = 2*(cd+ab); 
    R(2,2) = aa - bb - cc + dd;

    //Eigen is transposed relative to numpy
    r2 = r2 * R;

    for(int i=0; i<r1.rows(); i++)
    {
        double diff = (r2.row(i) - r1.row(i)).norm();
        if( diff > max_diff)
        {
            return false;
        }
    }
    return true;
}

double helper_functions::maxAtomMotion(const AtomMatrix v1) {
    double max = 0.0;
    for(int i = 0; i < v1.rows(); i++) {
        double norm = v1.row(i).norm();
        if(max < norm) {
            max = norm;
        }
    }    
    return max;
}

double helper_functions::maxAtomMotionV(const VectorXd v1) {
    double max = 0.0;
    for(int i=0;i<v1.rows();i+=3) {
        double norm = v1.segment<3>(i).norm();
        if(max < norm) {
            max = norm;
        }
    }    
    return max;
}

long helper_functions::numAtomsMoved(const AtomMatrix v1, double cutoff) {
    long num = 0;
    for(int i = 0; i < v1.rows(); i++) {
        double norm = v1.row(i).norm();
        if(norm >= cutoff) {
            num += 1;
        }
    }    
    return num;
}

AtomMatrix helper_functions::maxAtomMotionApplied(const AtomMatrix v1, double maxMotion)
{
    /*
    Function ensures (by scaling) that there is no single element of the AtomMatrix which is larger than maxMotion.
    */
    AtomMatrix v2(v1);
    
    double max = maxAtomMotion(v1);
    //double max = v1.norm();
    if(max > maxMotion)
    {
        v2 *= maxMotion/max;
    }
    return v2;
}

VectorXd helper_functions::maxAtomMotionAppliedV(const VectorXd v1, double maxMotion)
{
    VectorXd v2(v1);

    double max = maxAtomMotionV(v1);
    if (max > maxMotion)
    {
        v2 *= maxMotion/max;
    }
    return v2;
}

void helper_functions::getTime(double *real, double *user, double *sys)
{
    #ifdef WIN32
        *real = (double)time(NULL);
        if(user != NULL)
        {
            *user = 0.0;
        }
        if(sys != NULL)
        {
            *sys = 0.0;
        }
    #else
        struct timeval time;
        gettimeofday(&time, NULL);
        *real = (double)time.tv_sec + (double)time.tv_usec/1000000.0;
        struct rusage r_usage;
        if (getrusage(RUSAGE_SELF, &r_usage)!=0) 
        {
            fprintf(stderr, "problem getting usage info: %s\n", strerror(errno));
        }
        if(user != NULL)
        {
            *user = (double)r_usage.ru_utime.tv_sec + (double)r_usage.ru_utime.tv_usec/1000000.0;
        }
        if(sys != NULL)
        {
            *sys = (double)r_usage.ru_stime.tv_sec + (double)r_usage.ru_stime.tv_usec/1000000.0;
        }
    #endif
}

bool helper_functions::existsFile(string filename)
{
    FILE *fh;
    fh = fopen(filename.c_str(), "rb");
    if (fh == NULL) {
        return 0;
    }else fclose(fh);
    return 1;
}

string helper_functions::getRelevantFile(string filename)
{
    string filenameRelevant;
    string filenamePrefix;
    string filenamePostfix;

    // check if the _checkpoint version of the file is present
    int i = filename.rfind(".");
    filenamePrefix.assign(filename, 0, i);
    filenamePostfix.assign(filename, i, filename.size());
    filenameRelevant = filenamePrefix + "_checkpoint" + filenamePostfix;
    if(existsFile(filenameRelevant)) {
        return filenameRelevant;
    }
    // check if the _passed version of the file is present
    filenameRelevant = filenamePrefix + "_passed" + filenamePostfix;
    if(existsFile(filenameRelevant)) {
        return filenameRelevant;
    }
    // otherwise return original filename
    return filename;
}

AtomMatrix helper_functions::loadMode(FILE *modeFile, int nAtoms)
{
    AtomMatrix mode;
    mode.resize(nAtoms, 3);
    mode.setZero();
    for (int i=0; i < nAtoms; i++)
    {
        fscanf(modeFile, "%lf %lf %lf", &mode(i,0), &mode(i,1), &mode(i,2));
    }
    return mode;
}

AtomMatrix helper_functions::loadMode(string filename, int nAtoms)
{
    FILE *modeFile;
    modeFile = fopen(filename.c_str(), "rb");
    if (!modeFile) {
        cerr << "File " << filename << " was not found.\n";
        log("Stop\n");
        exit(1);
    }
    AtomMatrix mode = loadMode(modeFile, nAtoms);
    fclose(modeFile);
    return mode;
}

void helper_functions::saveMode(FILE *modeFile, Matter *matter, AtomMatrix mode)
{
    long const nAtoms = matter->numberOfAtoms();
    for (long i=0; i < nAtoms; ++i) {
        if (matter->getFixed(i)) {
            fprintf(modeFile, "0 0 0\n");
        }
        else {
            fprintf(modeFile, "%lf\t%lf \t%lf\n", mode(i,0), mode(i,1), mode(i,2));
        }
    }
    return;
}

std::vector<int> helper_functions::split_string_int(std::string s, std::string delim)
{
    char *pch;
    char *str;
    str = (char*) malloc(sizeof(char)*(s.length()+1));
    s.copy(str, s.length(), 0);
    str[s.length()] = '\0';
    cout << "COrrect: " << s << endl;
    printf("HI: %s\n", str);
    pch = strtok(str, delim.c_str());

    std::vector<int> list;

    while (pch != NULL) {
        char *endptr;
        int value = (int)strtol(pch, &endptr, 10); 
        if (strcmp(pch, endptr) == 0) {
            //Error reading
            std::vector<int> emptylist;
            free(str);
            return emptylist;
        }
        list.push_back(value);
        pch = strtok(NULL, delim.c_str());
    }
    free(str);
    return list;
}
