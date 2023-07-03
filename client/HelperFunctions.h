#ifndef HELPER_FUNCTIONS_H
#define HELPER_FUNCTIONS_H

#include "Eigen.h"
#include "Matter.h"
#include <string>
#include <vector>

// Random number generator constants

#define IM 2147483647
#define AM (1.0 / IM)
#define NTAB 32
#define NDIV (1 + (IM - 1) / NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)
#define IM1 2147483563
#define IM2 2147483399
#define IMM1 (IM1 - 1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791

/* Collection of supporting functions that handle arrays of doubles as vectors
 * and different random number generators */
namespace helper_functions {

double
random(long newSeed = 0);     // random number generator from numerical recipies
double randomDouble();        // random value between 0 and 1
double randomDouble(int max); // random value between 0 and max
double randomDouble(long max);   // random value between 0 and max
double randomDouble(double max); // random value between 0 and max
long randomInt(int lower, int upper);
double gaussRandom(double avg,
                   double std); // Gaussion random number with avg and std;
double dot(const double *v1, const double *v2, long size); // dot product
double length(const double *v1, long size); // length of vector v1
void add(double *result, const double *v1, const double *v2,
         long size); // v1 + v2
void subtract(double *result, const double *v1, const double *v2,
              long size); // v1 - v2
void multiplyScalar(double *result, const double *v1, double scalar,
                    long size); // scalar * v1
void divideScalar(double *result, const double *v1, double scalar,
                  long size); // (1/scalar) * v1
void copyRightIntoLeft(double *result, const double *v1,
                       long size);     // copy v2 into v1
void normalize(double *v1, long size); // v1 / |v1|
AtomMatrix makeOrthogonal(
    const AtomMatrix v1,
    const AtomMatrix v2); // return orthogonal component of v1 from v2
void makeProjection(double *result, const double *v1, const double *v2,
                    long size); // result = projection of v1 on v2
RotationMatrix rotationExtract(const AtomMatrix r1, const AtomMatrix r2);
bool rotationMatch(const Matter &m1, const Matter &m2, const double max_diff);
void rotationRemove(const AtomMatrix r1, std::shared_ptr<Matter> m2);
void rotationRemove(const std::shared_ptr<Matter> m1,
                    std::shared_ptr<Matter> m2);
void translationRemove(Matter &m1, const AtomMatrix r1);
void translationRemove(Matter &m1, const Matter &m2);
double maxAtomMotion(const AtomMatrix v1);
double maxAtomMotionV(const VectorXd v1);
long numAtomsMoved(const AtomMatrix v1, double cutoff);
AtomMatrix maxAtomMotionApplied(const AtomMatrix v1, double maxMotion);
VectorXd maxAtomMotionAppliedV(const VectorXd v1, double maxMotion);
AtomMatrix maxMotionApplied(const AtomMatrix v1, double maxMotion);
VectorXd maxMotionAppliedV(const VectorXd v1, double maxMotion);
void getTime(double *real, double *user, double *sys);
bool existsFile(string filename); // does filename exist
string
getRelevantFile(string filename); // return filename containing _checkpoint or
                                  // _passed if such a file exists
VectorXd loadMasses(string filename, int nAtoms);
AtomMatrix loadMode(FILE *modeFile, int nAtoms);
AtomMatrix loadMode(string filename, int nAtoms);
void saveMode(FILE *modeFile, std::shared_ptr<Matter> matter, AtomMatrix mode);
std::vector<int> split_string_int(std::string s, std::string delim);

bool identical(const Matter &m1, const Matter &m2,
               const double distanceDifference);
bool sortedR(const Matter &m1, const Matter &m2,
             const double distanceDifference);
void pushApart(std::shared_ptr<Matter> m1, double minDistance);

} // namespace helper_functions
#endif
