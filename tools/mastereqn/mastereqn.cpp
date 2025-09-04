#include <iostream>
#include <stdio.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <Eigen/MPRealSupport>
#include <Eigen/Dense>
#include <Eigen/Sparse>

using namespace std;
using namespace mpfr;
using namespace Eigen;

typedef Matrix<mpreal,Dynamic,Dynamic> MatrixXmp;
typedef Matrix<mpreal,Dynamic,1> VectorXmp;


void getTime(double *real, double *user, double *sys)
{
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
}

MatrixXmp readRateMatrix(string filename) {
    FILE *fh = fopen(filename.c_str(), "r");

    char buff[128];
    fgets(buff, 128, fh);
    int N = atoi(buff);
    MatrixXmp rateMatrix = MatrixXmp(N,N);

    for (int i=0;i<N;i++) {
        for (int j=0;j<N;j++) {
            double rate;
            int records = fscanf(fh, "%lg", &rate);
            if (records != 1) {
                fprintf(stderr, "error reading file: row %i col %i\n", i, j);
                abort();
            }
            rateMatrix(i,j) = rate;
        }
    }

    fclose(fh);
    return rateMatrix;
}

int main()
{
    // set precision to 256 bits (double has only 53 bits)
    mpreal::set_default_prec(256);
    double realTime0, realTime1;

    getTime(&realTime0, NULL, NULL);
    fprintf(stderr, "reading in rate_matrix\n");
    MatrixXmp rateMatrix = readRateMatrix("rate_matrix");
    getTime(&realTime1, NULL, NULL);
    fprintf(stderr, "took %.3f seconds\n", realTime1-realTime0);
    int N = rateMatrix.rows();
    fprintf(stderr, "read in %ix%i rate matrix\n", N, N);

    fprintf(stderr, "solving generalized eigen problem\n");
    getTime(&realTime0, NULL, NULL);
    EigenSolver<MatrixXmp> eigensolver(rateMatrix);
    getTime(&realTime1, NULL, NULL);
    fprintf(stderr, "took %.3f seconds\n", realTime1-realTime0);
    if (eigensolver.info() != Success) {
        fprintf(stderr, "eigensolver failed\n");
        abort();
    }
    VectorXmp ew = eigensolver.eigenvalues().real();
    MatrixXmp ev = eigensolver.eigenvectors().real();

    VectorXmp p0 = VectorXmp::Zero(N);
    p0(0) = 1.0;
    fprintf(stderr, "converting initial probability vector into eigenvector basis\n");
    getTime(&realTime0, NULL, NULL);
    VectorXmp c0 = ev.colPivHouseholderQr().solve(p0);
    getTime(&realTime1, NULL, NULL);
    fprintf(stderr, "took %.3f seconds\n", realTime1-realTime0);

    printf("set logscale x\n");
    printf("plot \"-\" w l\n");

    int final_state = 9;
    for (int p=15;p>=0;p--) {
        mpreal t = pow(10,-p);
        mpreal prob = 0.0;
        for (int i=0;i<N;i++) {
            prob += c0(i)*exp(-ew(i)*t)*ev(final_state,i);
        }
        cout << t << " " << prob << endl;
    }

    return 0;
}
