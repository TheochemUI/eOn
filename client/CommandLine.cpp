#include "CommandLine.h"
#include "ConjugateGradients.h"
#include "Matter.h"
#include "Parameters.h"
#include "Log.h"
#include "version.h"

#include <cstdlib>
#include <string>
#include <unistd.h>

using namespace std;

void singlePoint(Parameters *parameters, Matter *matter)
{
    printf("Energy:         %.10f\n", matter->getPotentialEnergy());
    printf("Max atom force: %.10g\n", matter->maxForce());
}

void minimize(Parameters *parameters, Matter *matter, string confileout)
{
    matter->relax(false, false);
    if (confileout.length() > 0) {
        printf("saving relaxed structure to %s\n", confileout.c_str());
    }else{
        printf("no output file specified, not saving\n");
    }
    matter->matter2con(confileout);
}

void usage(void)
{
    fprintf(stderr, "Usage: eonclient [options] inputConfile [outputConfile]\n");
    char fmtStr[] = "  -%-2s %s\n";

    fprintf(stderr, "Job Type:\n");
    fprintf(stderr, fmtStr, "v", "print version information");
    fprintf(stderr, fmtStr, "m", "Minimization of inputConfile saves to outputConfile");
    fprintf(stderr, fmtStr, "s", "Single point energy of inputConfile");
    fprintf(stderr, fmtStr, "c", "Compare structures of inputConfile to outputConfile");
    fprintf(stderr, fmtStr, "o", "Optimization method [default: qm]");
    fprintf(stderr, fmtStr, "f", "Convergence force [default: 0.001]");
    fprintf(stderr, fmtStr, "t", "Distance tolerance [default: 0.1]");

    fprintf(stderr, "Required Options:\n");
    fprintf(stderr, fmtStr, "p", "The potential (e.g. qsc, lj, eam_al)");
}

void commandLine(int argc, char **argv)
{
    // no getopt on windows
    #ifndef WIN32
    int c;
    bool sflag = false, mflag = false, pflag = false, cflag = false;
    double optConvergedForce = 0.001;

    string potential;
    string confile;
    string optimizer("cg");

    Parameters *parameters = new Parameters;

    while ((c=getopt(argc,argv,"chsmp:f:o:t:v")) != -1) {
        switch (c) {
            case 'c':
                cflag = true;
                break;
            case 's':
                sflag = true;
                break;
            case 'm':
                mflag = true;
                break;
            case 'p':
                pflag = true;
                potential = optarg;
                break;
            case 't':
                parameters->distanceDifference = atof(optarg);
                break;
            case 'o':
                optimizer = optarg;
                break;
            case 'f':
                cout <<optarg<<endl;
                optConvergedForce = atof(optarg);
                break;
            case 'h':
                usage();
                exit(0);
            case 'v':
                printf("eonclient version r%s\n", VERSION);
                printf("          compiled %s\n", BUILD_DATE);
                exit(0);
            case '?':
                if (optopt == 'p')
                    fprintf(stderr, "Option -%c requires an argument.\n", optopt);
                else
                    fprintf(stderr, "Unknown option `-%c'.\n", optopt);
                usage();
                exit(2);
        }
    }

    if (sflag && mflag) {
        fprintf(stderr, "Cannot specify both minimization and single point\n");
        exit(2);
    }

    if (!pflag && (sflag || mflag)) {
        fprintf(stderr, "Must specify a potential\n");
        exit(2);
    }else if (!cflag) {
        for (string::size_type i = 0; i < potential.length(); ++i) {
          potential[i] = tolower(potential[i]);
        }
    }

    int extraArgs = argc - optind;

    if (extraArgs < 1) {
        fprintf(stderr, "Only one non-option argument is allowed: the con file\n");
        exit(2);
    }else{
        confile = argv[optind];
    }

    if (!cflag) {
        parameters->potential = potential;
    }
    parameters->optMethod = optimizer;
    parameters->optConvergedForce = optConvergedForce;

    log_init(parameters, (char*)"client.log");

    Matter *matter = new Matter(parameters);
    Matter *matter2 = new Matter(parameters);
    matter->con2matter(confile);

    string confileout;
    if (extraArgs == 2) {
        confileout = argv[optind+1];
        if (cflag) matter2->con2matter(confileout);
    }

    if (sflag) {
        singlePoint(parameters, matter);
    }else if (mflag) {
        minimize(parameters, matter, confileout);
    }else if (cflag) {
        parameters->checkRotation = true;
        if (matter->compare(matter2, true)) {
            printf("structures match\n");
        }else{
            printf("structures do not match\n");
        }
    }

    delete parameters;
    delete matter;
    delete matter2;
    #endif
}
