/*
 *===============================================
 *  Modified. Name, Date and a small description!
 *
 *-----------------------------------------------
 *  Todo:
 *
 *===============================================
 */
#include <ctime>
#include <cstring>
#include <cassert>
#ifdef NDEBUG
      #include <boinc/boinc_api.h>
      #include <boinc/diagnostics.h>     // boinc_init_diagnostics()
      #include <boinc/filesys.h>         // boinc_fopen(), etc...
#else
      #include "false_boinc.h"
#endif
#include "common/Parameters.h"

#include "ClientEON.h"
#include "LowestEigenmodeInterface.h"  


using namespace constants;
using namespace client_eon;

int main(int argc, char **argv) 
{
    unsigned const size=255;
    char parameters_passed[size], reactant_passed[size], displacement_passed[size], mode_passed[size];
    if (argc > 1 ) {
        strncpy(parameters_passed, argv[1], sizeof(parameters_passed));
        strncpy(reactant_passed, argv[2], sizeof(reactant_passed));
        strncpy(displacement_passed, argv[3], sizeof(displacement_passed));
        strncpy(mode_passed, argv[4], sizeof(mode_passed));
    }
    else {
        strncpy(parameters_passed, PARMS_PASSED_FILE_NAME.c_str(), sizeof(parameters_passed));
        strncpy(reactant_passed, REAC_PASSED_FILE_NAME.c_str(), sizeof(reactant_passed));
        strncpy(displacement_passed, DISPLACEMENT_PASSED_FILE_NAME.c_str(), sizeof(displacement_passed));
        strncpy(mode_passed, MODE_PASSED_FILE_NAME.c_str(), sizeof(mode_passed));
    };

    // BOINC is started up
    // before initializing BOINC itself, intialize diagnostics, so as
    // to get stderr output to the file stderr.txt, and thence back home.
    //boinc_init_diagnostics(BOINC_DIAG_REDIRECTSTDERR|
    //                       BOINC_DIAG_REDIRECTSTDOUT|
    //                       BOINC_DIAG_MEMORYLEAKCHECKENABLED|
    //                       BOINC_DIAG_DUMPCALLSTACKENABLED| 
    //                       BOINC_DIAG_TRACETOSTDERR);
    
    rc = boinc_init();
    if(rc) 
    {
        exit(rc);
    }
    
    // Checks if the user has requested information
    if ((argc > 1) && (*argv[1] == '-')) 
    {    
        printRequestedInfo(argv[1]);
        boinc_finish(0);
        return 0;
    }


    // Loads runtime parameters and relaxes the initial configuration.
    loadDataAndRelax(parameters_passed, reactant_passed);
    
    // If we are not only minimizing, perform saddle searches.
    if(!parameters.getMinimizeOnly())
    {
        // If the server performed the displacement, load it and perform the saddle search.
        if (parameters.getRefineSP()) 
        {
            // While we can read in the displacement and mode files, perform a saddle search on them.
            while (loadDisplacementAndMode(displacement_passed, mode_passed) == true) 
            {
                doSaddleSearch();
            }
        }
        // Otherwise, displace and search ourselves.
        else
        {
            doSaddleSearch();
        }
    }
    else
    {
        saveData();
    }
    // BOINC applications must exit via boinc_finish(rc), not merely exit */
    boinc_finish(0);
    return 0;
}

void client_eon::doSaddleSearch(void)
{
    bool continueDimerSearch=true;
    int currentTry=0;
    // ----- start dimer search ----- 
    do
    {
        fprintf(stdout, "\nNew search started.\n");

        initializeNewSearch();
        // state is equal to getStateInit() if a saddle point is determined
        state = saddlePoint.locate(min1, min2);
        
        // Barrier determination part
        if(state == getStateInit() && connectedToReactantState())
        {
            determineBarriers();
            barriersOK = barriersWithinWindow();
        }
        if(!parameters.getPrefactorsTag() && barriersOK)
            state = getStateGood();

        // Prefactor determination part
        if(state==getStateInit() && parameters.getPrefactorsTag())
        {
            prefactors.compute(prefactorsValues);
            prefactorsOK = prefactorsWithinWindow();
        }
        if(parameters.getPrefactorsTag() && prefactorsOK)
            state = getStateGood();

        currentTry++;

        if (parameters.getRefineSP()) {
            continueDimerSearch = false;
            continue;
        }
            
        if ( (state == getStateGood()) ||
             // getStateSaddleSearchNoConvexRegion() should be treated special 
             // as the saddle point search does not starts since the initial 
             // displacement is too small
             (state == getStateSaddleSearchNoConvexRegion()) ) {
            continueDimerSearch = false;
        }

        if (currentTry > parameters.getNrOfTriesToDetermineSaddlePoints_SP()) {
            printf("Did maximum number (%i) of saddle searches.", currentTry);
            continueDimerSearch = false;
        }
        printEndState(state);

    } while (continueDimerSearch);
    // ------ end simulation ------ 

    saveData();
    parameters.printOutput();    
}

static FILE * openFile(char const name[])
{
    char resolved[512];
    int rc = boinc_resolve_filename(name, resolved, sizeof(resolved));
    if (rc) {
        fprintf(stderr, "Cannot resolve intput file: %s, RC=%d\n", name, rc);
        boinc_finish(rc);    // back to BOINC core
    }
    FILE * file = boinc_fopen(resolved, "r");
    if (file == NULL){
        fprintf(stderr, "Cannot open file: %s\n", resolved);
        boinc_finish(rc);    // back to BOINC core
    }
    return file;
}

static void closeFile(FILE * file, char const name[])
{
        if (ferror(file)) {
                fprintf(stderr, "Error while reading %s\n", name);
                boinc_finish(rc);
        };
        fclose(file);
}

void client_eon::loadDataAndRelax(char const parameters_passed[], char const reactant_passed[])
{
    FILE * file;
    file = openFile(parameters_passed);
    parameters.load(file);
    closeFile(file, parameters_passed);
    
    // Initialize random generator
    if(parameters.getRandomSeed() < 0) 
    {
        unsigned i = time(NULL);
        parameters.setRandomSeed(i);
        srand(i);
    }
    else
    {
        srand(parameters.getRandomSeed());
    }
    printf("Random seed is: %ld\n", parameters.getRandomSeed());

    // The paramters have been loaded && the Matter objects can be initialized
    initial = new Matter(&parameters);
    displacement = new Matter(&parameters);
    saddle = new Matter(&parameters);
    min1 = new Matter(&parameters);
    min2 = new Matter(&parameters);

    file = openFile(reactant_passed);
    int ok = initial->con2matter(file);
    if (not ok) 
    {
        fprintf(stderr, "Cannot read file %s\n", reactant_passed);
        std::exit(1);
    };
    closeFile(file, reactant_passed);
    
    // Relax the passed configuration.
    // RT: if we are only minimizing AND are minimizing the box, use QMBox
    if (parameters.getMinimizeOnly() && parameters.getMinimizeBox())
    {
        QMBox qmboxInitial(initial, &parameters);
        qmboxInitial.fullRelax();
    }
    // RT: Otherwise, minimize as usual.
    else
    {    
        ConjugateGradients cgInitial(initial, &parameters);
        cgInitial.fullRelax();
    }

    return;
}


bool client_eon::loadDisplacementAndMode(char const displacement_passed[], char const mode_passed[])
{
    static int modePosition=0;
    static int displacementPosition=0;
    static int count=0;
    static int displacementFileLength=0;
    FILE * file;
    char buf[1024];

    //read in displacement con file
    file=openFile(displacement_passed);
    if (count == 0) {
        fseek(file, 0, SEEK_END);
        displacementFileLength = ftell(file);
        fseek(file, 0, SEEK_SET);
    }else{
        if (displacementFileLength <= displacementPosition) {
            return false;
        }
        fseek(file, displacementPosition, SEEK_SET);
    }
    saddle->con2matter(file);
    //read in displacement con file
    file=openFile(displacement_passed);
    if (count == 0) {
        fseek(file, 0, SEEK_END);
        displacementFileLength = ftell(file);
        fseek(file, 0, SEEK_SET);
    }else{
        if (displacementFileLength <= displacementPosition) {
            return false;
        }
        fseek(file, displacementPosition, SEEK_SET);
    }
    displacement->con2matter(file);
    displacementPosition = ftell(file);
    closeFile(file, displacement_passed);

    //read in mode file
    file=openFile(mode_passed);
    fseek(file, modePosition, SEEK_SET);
    long nall=0, nfree=0;
   
    //fgets(buf, 1024, file);
    //printf(buf);
    //fseek(file, 0, SEEK_END);
    //printf("length %i\n", ftell(file));
    //fseek(file, modePosition, SEEK_SET);
    //printf("seeked to %i\n", modePosition);
    fscanf(file, "%ld %ld", &nall, &nfree);
    assert(nall == 3*saddle->numberOfAtoms());
    assert(nfree == 3*saddle->numberOfFreeAtoms());
    client_eon::eigenvector.size=nfree;
    double * vector=new double[nall];
    client_eon::eigenvector.vector=vector;
    for (int i=0, j=0; i < nall; ++i) {
        //if (not saddle->getFixed(i/3)) {
            fscanf(file, "%lf", &vector[j]);
            ++j;
        //};
    };
    modePosition = ftell(file);
    closeFile(file, mode_passed);

    count++;
    return true;
}


void client_eon::initializeNewSearch()
{
    barriersValues[0] = barriersValues[1] = 0;
    prefactorsValues[0] = prefactorsValues[1] = 0;
    barriersOK = prefactorsOK = false;
    if (not parameters.getRefineSP()) {
        // Matter objects
        *saddle = *min1 = *min2 = *initial;
    }
    else {        
        *min1 = *min2 = *initial;
    }
    // SaddlePoint object
    saddlePoint.initialize(initial, saddle, &parameters);
    // If prefactor has to be determined Prefactors object
    if(parameters.getPrefactorsTag())
    {
        prefactors.initialize(saddle, min1, min2, &parameters);
    }
    return;
}


bool client_eon::connectedToReactantState(){
    bool result = false;
    Matter matterTemp(&parameters);
    
    // Ensure that not both minima correspond to the initial state
    if((*initial==*min1) && (!(*initial==*min2)))
        result = true;
    // If min2 corresponds to initial state swap min1 && min2
    else if(!(*initial==*min1) && ((*initial==*min2))){
        matterTemp = *min1;
        *min1 = *min2;
        *min2 = matterTemp;
        result = true;
    }
    else{
        result = false;
        state = getStateNotConnected();
    }
    return result;
}


void client_eon::determineBarriers()
{ 
    barriersValues[0] = saddle->potentialEnergy()-min1->potentialEnergy();
    barriersValues[1] = saddle->potentialEnergy()-min2->potentialEnergy();
    return;
}


bool client_eon::barriersWithinWindow()
{
    bool result = true;

    if((parameters.getMaxEnergy_SP() < barriersValues[0]) || (parameters.getMaxEnergy_SP() < barriersValues[1]))
    {
        result = false;
        state = getStateBadBarrier();
    }
    return result;
}


bool client_eon::prefactorsWithinWindow(){
    bool result = true;
    
    if((getPrefactorMax()<prefactorsValues[0]) || 
       (prefactorsValues[0]<getPrefactorMin())){
        result = false;
        state = getStateBadPrefactor();
    }
    if((getPrefactorMax()<prefactorsValues[1]) || 
       (prefactorsValues[1]<getPrefactorMin())){
        result = false;
        state = getStateBadPrefactor();
    }
    return result;
}


void client_eon::saveData(){
    FILE* fileResults;
    FILE* fileReactant;
    FILE* fileSaddle;
    FILE* fileProduct;
    FILE* fileMode;

    char resolvedNameResults[STRING_SIZE];
    char resolvedNameReactant[STRING_SIZE];
    char resolvedNameSaddle[STRING_SIZE];
    char resolvedNameProduct[STRING_SIZE];
    char resolvedNameMode[STRING_SIZE];

    parameters.setPotentialEnergySP(saddle->potentialEnergy());
    parameters.setPotentialEnergyMin1(min1->potentialEnergy());
    parameters.setPotentialEnergyMin2(min2->potentialEnergy());
    
    parameters.setPrefactorReac_Prod(prefactorsValues[0]);
    parameters.setPrefactorProd_Reac(prefactorsValues[1]);
    
    parameters.setBarrierReac_Prod(barriersValues[0]);
    parameters.setBarrierProd_Reac(barriersValues[1]);

    parameters.addForceCalls(parameters.getForceCallsSaddlePoint()+
                             parameters.getForceCallsPrefactors());
    
    parameters.setTerminationReason(state);
    parameters.setDisplacementSaddleDistance(displacement->per_atom_norm(*saddle));
    
    // Input need to be "resolved" from their logical name
    // for the application to the actual path on the client's disk
    rc = boinc_resolve_filename(RESULTS_FILE_NAME.c_str(), 
                                resolvedNameResults, 
                                sizeof(resolvedNameResults));
    if (rc){
        fprintf(stderr, "Cannot resolve output file name. RC=%d\n", rc);
        boinc_finish(rc);    // back to BOINC core
    }
    rc = boinc_resolve_filename(REAC_FILE_NAME.c_str(), 
                                resolvedNameReactant, 
                                sizeof(resolvedNameReactant));
    if (rc){
        fprintf(stderr, "Cannot resolve output file name. RC=%d\n", rc);
        boinc_finish(rc);    // back to BOINC core
    }
    rc = boinc_resolve_filename(SEND_SP_CONF_FILE_NAME.c_str(), 
                                resolvedNameSaddle, 
                                sizeof(resolvedNameSaddle));
    if (rc){
        fprintf(stderr, "Cannot resolve output file name. RC=%d\n", rc);
        boinc_finish(rc);    // back to BOINC core
    }
    rc = boinc_resolve_filename(SEND_PROD_FILE_NAME.c_str(), 
                                resolvedNameProduct, 
                                sizeof(resolvedNameProduct));
    if (rc){
        fprintf(stderr, "Cannot resolve output file name. RC=%d\n", rc);
        boinc_finish(rc);    // back to BOINC core
    }
    rc = boinc_resolve_filename(MODE_FILE_NAME.c_str(), 
                                resolvedNameMode, 
                                sizeof(resolvedNameMode));
    if (rc){
        fprintf(stderr, "Cannot resolve output file name. RC=%d\n", rc);
        boinc_finish(rc);    // back to BOINC core
    }
        
    fileResults = boinc_fopen(resolvedNameResults, "a");
    fileReactant = boinc_fopen(resolvedNameReactant, "a");
    fileSaddle = boinc_fopen(resolvedNameSaddle, "a");
    fileProduct = boinc_fopen(resolvedNameProduct, "a");
    fileMode = boinc_fopen(resolvedNameMode, "a");

    parameters.saveOutput(fileResults);

    if(parameters.getMinimizeOnly())
    {
        initial->matter2con(fileReactant);
    }
    else
    {
        min1->matter2con(fileReactant);
    }
    saddle->matter2con(fileSaddle);    
    min2->matter2con(fileProduct);
    long const nfree=saddlePoint.getnFreeCoord();
    long const nAtoms=saddle->numberOfAtoms();
    fprintf(fileMode, "%ld %ld\n", nAtoms*3, nfree);
    double const *const eigenMode=saddlePoint.getEigenMode();
    for (long i=0, j=0; i < nAtoms; ++i) {
        if (saddle->getFixed(i)) {
            fprintf(fileMode, "0 0 0\n");
        }
        else {
            fprintf(fileMode, "%lf\t%lf \t%lf\n", eigenMode[j], eigenMode[j+1], eigenMode[j+2]);
            j+=3;
        };
    };
    fclose(fileResults);               
    fclose(fileReactant);  
    fclose(fileSaddle);               
    fclose(fileProduct);
    fclose(fileMode);
    return;
}

void client_eon::printEndState(long state){
    fprintf(stdout, "Final state: ");
    if(state == getStateGood())
        fprintf(stdout, "Succesful.\n");

    else if(state == getStateSaddleSearchNoConvexRegion())
        fprintf(stdout, "Initial displacement, not able to reach convex region.\n");
   
    else if(state == getStateSaddleSearchTerminatedBarrier())
        fprintf(stdout, "Saddle search, too high barrier.\n");
        
    else if(state == getStateSaddleSearchTerminatedConcaveIterations()) 
        fprintf(stdout, "Saddle search, a too long series of iterations in concave region.\n");

    else if(state == getStateSaddleSearchTerminatedTotalIterations())
        fprintf(stdout, "Saddle search, too many iterations in saddle point search.\n");

    else if(state == getStateNotConnected())
        fprintf(stdout, "Minima, saddle is not connected to initial state.\n");

    else if(state == getStateBadPrefactor())
            fprintf(stdout, "Prefactors, not within window as defined in Constants\n");

    else if(state == getStateBadBarrier())
        fprintf(stdout, "Energy barriers, not within window as defined in Constants\n");

    else
        fprintf(stdout, "Unknown state!!!\n");
    return;
}


// Just to print information to the user
void client_eon::printRequestedInfo(char *argv) {
    if (strcmp(argv, "-version") == 0 || strcmp(argv, "--version") == 0) {
        fprintf(stdout, "Version: %s\n",VERSION_INFO_EON2.c_str());
    }
    else if (strcmp(argv, "-print_parameters") == 0 || 
             strcmp(argv, "--print_parameters") == 0) {
        fprintf(stdout, "\n-------------------\n");
        fprintf(stdout, "Parameters used when executed.\n");
        fprintf(stdout, "If \"Bus error\" parameter file is not in directory (-make_parameters creates file).\n\n");
        FILE* fileParameters;
        fileParameters = fopen(PARMS_PASSED_FILE_NAME.c_str(), "r");            
        parameters.load(fileParameters);
        parameters.printInput();
        fprintf(stdout, "-------------------\n\n");
    }
    else if (strcmp(argv, "-make_parameters") == 0 || 
             strcmp(argv, "--make_parameters") == 0) {
        fprintf(stdout, "\n-------------------\n");
        fprintf(stdout, "Create needed parameter file.\n");
        FILE* fileParameters;
        fileParameters = fopen(PARMS_PASSED_FILE_NAME.c_str(), "w");            
        parameters.saveInput(fileParameters);
        fprintf(stdout, "-------------------\n\n");
    }
    else if (strcmp(argv, "-files") == 0 || strcmp(argv, "--files") == 0) {
        fprintf(stdout, "\n-------------------\n");
        fprintf(stdout, "The following files have to exist in the same directory as the this executable!\n\n"); 
        fprintf(stdout, "Input reactant configuration:           %s\n",
                REAC_PASSED_FILE_NAME.c_str());            
        fprintf(stdout, "Input (-make_parameters creates file):  %s\n",
                PARMS_PASSED_FILE_NAME.c_str());            
        fprintf(stdout, "Output (can be empty):                  %s\n",
                REAC_FILE_NAME.c_str());
        fprintf(stdout, "Output (can be empty):                  %s\n",
                SEND_SP_CONF_FILE_NAME.c_str());
        fprintf(stdout, "Output (can be empty):                  %s\n",
                SEND_PROD_FILE_NAME.c_str()); 
        fprintf(stdout, "Output (can be empty):                  %s\n",
                RESULTS_FILE_NAME.c_str());
        fprintf(stdout, "-------------------\n");
        fprintf(stdout, "Files created by BOINC.\n\n");
        fprintf(stdout, "stderr.txt\n");
        fprintf(stdout, "init_data.xml\n");
        fprintf(stdout, "boinc_finish_called\n");
        fprintf(stdout, "-------------------\n\n");
    }
    else if (strcmp(argv, "-help") == 0 || strcmp(argv, "--help") == 0 || 
             strcmp(argv, "-h") == 0 || strcmp(argv, "--h") == 0) {
        fprintf(stdout, "\n-------------------\n");
        fprintf(stdout, "-help               prints this.\n");
        fprintf(stdout, "-version            prints version number.\n");
        fprintf(stdout, "-print_parameters   prints parameters used during execution.\n");
        fprintf(stdout, "-make_parameters    makes parameter file needed. Overwrites existing file!\n");
        fprintf(stdout, "-files              prints the files needed for execution.\n");
        fprintf(stdout, "-------------------\n\n");
    }
    else{
        fprintf(stdout, "\n-------------------\n");
        fprintf(stdout, "Unknown command! (-help prints known commands)\n");
        fprintf(stdout, "-------------------\n\n");
    }        
    return;
}
