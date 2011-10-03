#include "CommandLine.h"
#include "ConjugateGradients.h"
#include "Matter.h"
#include "Parameters.h"
#include "Log.h"

#include <cstdlib>
#include <string>

using namespace std;

void singlePoint(Parameters *parameters, Matter *matter)
{
    printf("Energy:         %f\n", matter->getPotentialEnergy());
    printf("Max atom force: %g\n", matter->maxForce());
}

void minimize(Parameters *parameters, Matter *matter, string confileout)
{
    ConjugateGradients cg(matter, parameters);
    cg.setOutput(1);
    cg.fullRelax();

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
    fprintf(stderr, fmtStr, "m", "Minimization of inputConfile saves to outputConfile");
    fprintf(stderr, fmtStr, "s", "Single point energy of inputConfile");

    fprintf(stderr, "Required Options:\n");
    fprintf(stderr, fmtStr, "p", "The potential (e.g. qsc, lj, eam_al)");
}

void commandLine(int argc, char **argv)
{
    // no getopt on windows
    #ifndef WIN32
    int c;
    bool sflag=false,mflag=false, pflag=false;

    string potential;
    string confile;

    while ((c=getopt(argc,argv,"hsmp:")) != -1) {
        switch (c) {
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
            case 'h':
                usage();
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

    if (!pflag) {
        fprintf(stderr, "Must specify a potential\n");
        exit(2);
    }else{
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

    Parameters *parameters = new Parameters;
    parameters->potential = potential;

    log_init(parameters, "client.log");

    Matter *matter = new Matter(parameters);
    matter->con2matter(confile);

    string confileout;
    if (extraArgs == 2) {
        confileout = argv[optind+1];
    }

    if (sflag) {
        singlePoint(parameters, matter);
    }else if (mflag) {
        minimize(parameters, matter, confileout);
    }

    delete parameters;
    delete matter;
    #endif
}
