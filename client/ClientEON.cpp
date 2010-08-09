/*
 *===============================================
 *  Modified. Name, Date and a small description!
 *
 *-----------------------------------------------
 *  Todo:
 *
 *===============================================
 */

#include "ClientEON.h"

using namespace constants;
using namespace client_eon;

int main(int argc, char **argv) 
{
	// BOINC is started
	rc = boinc_init();
	if(rc){
		exit(rc);
	}
	char parameters_passed[STRING_SIZE], reactant_passed[STRING_SIZE], displacement_passed[STRING_SIZE], mode_passed[STRING_SIZE];
    if (argc > 1 ) {
		if (*argv[1] == '-'){    
			if (strcmp(argv[1], "-test") == 0 || strcmp(argv[1], "--test") == 0) {				
				forcesOfConfig();
			}
			else{
				printRequestedInfo(argv[1]);   //// NEEDS TO BE UPDATED
			}
			// To prevent segmentation fault
			boinc_finish(0);
			return 0;

        }
        else { 
			// Uses data passed from user
			strncpy(parameters_passed, argv[1], sizeof(parameters_passed));
			strncpy(reactant_passed, argv[2], sizeof(reactant_passed));
			strncpy(displacement_passed, argv[3], sizeof(displacement_passed));
			strncpy(mode_passed, argv[4], sizeof(mode_passed));
		}//
		
    }
    else {
		// Uses data from Constant file
        strncpy(parameters_passed, PARMS_PASSED_FILE_NAME.c_str(), sizeof(parameters_passed));
        strncpy(reactant_passed, REAC_PASSED_FILE_NAME.c_str(), sizeof(reactant_passed));
        strncpy(displacement_passed, DISPLACEMENT_PASSED_FILE_NAME.c_str(), sizeof(displacement_passed));
        strncpy(mode_passed, MODE_PASSED_FILE_NAME.c_str(), sizeof(mode_passed));
    };
    
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
                try{doSaddleSearch();}
                catch(int except){}
                parameters.resetForceCalls();
                parameters.resetForceCallsSaddlePoint();
                parameters.resetForceCallsPrefactors();
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
	// To prevent segmentation fault
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


static FILE * openFile(char const name[], char const readWriteAppend[])
{
    char resolved[STRING_SIZE];
    int rc = boinc_resolve_filename(name, resolved, sizeof(resolved));
    if (rc) {
        fprintf(stderr, "Cannot resolve intput file: %s, RC=%d\n", name, rc);
        boinc_finish(rc);    // back to BOINC core
    }
    FILE * file = boinc_fopen(resolved, readWriteAppend);
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
	FILE *fileParameters, *fileReactant;
    fileParameters = openFile(parameters_passed, READ.c_str());
    parameters.load(fileParameters);
    closeFile(fileParameters, parameters_passed);
    
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

    // The parameters have been loaded && the Matter objects can be initialized
    initial = new Matter(&parameters);
    displacement = new Matter(&parameters);
    saddle = new Matter(&parameters);
    min1 = new Matter(&parameters);
    min2 = new Matter(&parameters);

    fileReactant = openFile(reactant_passed, READ.c_str());
    int ok = initial->con2matter(fileReactant);
    if (not ok) 
    {
        fprintf(stderr, "Cannot read file %s\n", reactant_passed);
        std::exit(1);
    };
    closeFile(fileReactant, reactant_passed);
    
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
	//Dangerous to reuse the same pointer name for different files!
    FILE * file;

    //read in displacement con file
    file=openFile(displacement_passed, READ.c_str());
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
    file=openFile(displacement_passed, READ.c_str());
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
    file=openFile(mode_passed, READ.c_str());
    fseek(file, modePosition, SEEK_SET);
    
    saddlePoint.loadMode(file);   /// new.
    
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
    
    // Ensure that both minima do not correspond to the initial state
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
    FILE *fileResults, *fileReactant, *fileSaddle, *fileProduct, *fileMode;

    if (parameters.getMinimizeOnly()) {
        //If we are only minimizing set the min1 pointer
        //to the matter object we just minimized so that
        //the energy can be saved to the results.dat file.
        min1 = initial;
    }

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

	fileResults = openFile(RESULTS_FILE_NAME.c_str(), APPEND.c_str());
    parameters.saveOutput(fileResults);
	closeFile(fileResults, RESULTS_FILE_NAME.c_str());

	fileReactant = openFile(REAC_FILE_NAME.c_str(), APPEND.c_str());
    if(parameters.getMinimizeOnly())
    {
        initial->matter2con(fileReactant);
    }
    else
    {
        min1->matter2con(fileReactant);
        // Only save the saddle in the case that we did a saddle search.
        fileMode = openFile(MODE_FILE_NAME.c_str(), APPEND.c_str());
        saddlePoint.saveMode(fileMode);
        closeFile(fileMode, MODE_FILE_NAME.c_str());
    }
	closeFile(fileReactant, REAC_FILE_NAME.c_str());  

	fileSaddle = openFile(SEND_SP_CONF_FILE_NAME.c_str(), APPEND.c_str());
    saddle->matter2con(fileSaddle);
	closeFile(fileSaddle, SEND_SP_CONF_FILE_NAME.c_str());               

	fileProduct = openFile(SEND_PROD_FILE_NAME.c_str(), APPEND.c_str());
	min2->matter2con(fileProduct);
    closeFile(fileProduct, SEND_PROD_FILE_NAME.c_str());

    return;
}


void client_eon::printEndState(long state){
    fprintf(stdout, "Final state: ");
    if(state == getStateGood())
        fprintf(stdout, "Succesful.\n");

    else if(state == getStateSaddleSearchNoConvexRegion())
        fprintf(stdout, "Initial displacement, not able to reach convex region.\n");
   
    else if(state == getStateSaddleSearchTerminatedBarrier())
        fprintf(stdout, "Saddle search, barrier too high.\n");
        
    else if(state == getStateSaddleSearchTerminatedConcaveIterations()) 
        fprintf(stdout, "Saddle search, too many iterations in concave region.\n");

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
        fprintf(stdout, "The following files have to exist in the same directory as this executable\n\n"); 
        fprintf(stdout, "Input reactant configuration:           %s\n",
                REAC_PASSED_FILE_NAME.c_str());            
        fprintf(stdout, "Input (-make_parameters creates file):  %s\n",
                PARMS_PASSED_FILE_NAME.c_str());            
        fprintf(stdout, "Input displacement done by server:      %s\n",
				DISPLACEMENT_PASSED_FILE_NAME.c_str());
		fprintf(stdout, "Input initial mode:                     %s\n",
				MODE_PASSED_FILE_NAME.c_str());
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
        fprintf(stdout, "-help               prints this\n");
        fprintf(stdout, "-version            prints version number\n");
        fprintf(stdout, "-print_parameters   prints parameters used during execution\n");
        fprintf(stdout, "-make_parameters    makes parameter file needed. Overwrites existing file!\n");
        fprintf(stdout, "-files              prints the files needed for execution\n");
        fprintf(stdout, "-test               prints potential energy and forces of reactant.test to potentialInfo.txt\n");
        fprintf(stdout, "-------------------\n\n");
    }
    else{
        fprintf(stdout, "\n-------------------\n");
        fprintf(stdout, "Unknown command! (-help prints known commands)\n");
        fprintf(stdout, "-------------------\n\n");
    }        
    return;
}

void client_eon::forcesOfConfig(){
	
	FILE *fileParameters, *fileReactant;
	
    fileParameters = fopen("parameters.test", "r");            
    parameters.load(fileParameters);
    fclose(fileParameters);
    
    testConfig = new Matter(&parameters);
    
    fileReactant = fopen("reactant.test", "r");
    //printf("test->con2matter(fileReactant);\n");
    int ok = testConfig->con2matter(fileReactant);   //load the reactant
    if (not ok) 
    {
        fprintf(stderr, "Cannot read reactant.test\n");
        std::exit(1);
    };
    fclose(fileReactant);
    //printf("Reactant loaded\n");
    ConjugateGradients cgInitial(testConfig, &parameters);
    //printf("Conjugate gradients initialized\n");
    
    double potEnergy;
    potEnergy = testConfig->potentialEnergy();
    //printf("Potential energy %lf\n", potEnergy);
    
    long int numOfAtoms;    
    numOfAtoms = testConfig->numberOfAtoms();
    
    //printf("Number of atoms %li\n",numOfAtoms);
    double initialForces[3*numOfAtoms];    // three coordinates for each atom
    testConfig->getForces(initialForces);
	FILE * file;
    file = fopen("potentialInfo.txt","w");
    fprintf(file, "Information about potential energy and forces of a\ncurrent configuration in reactant_passed.con\n\n\n\n");
    fprintf(file, "Potential tag: %li\n\n\n", parameters.getPotentialTag());
    fprintf(file, "Potential energy: %lf\n\n", potEnergy);
    fprintf(file, "Force coordinates:\n");
    fprintf(file, "Atom     %-13.13s%-13.13s%-13.13s\n", "x:","y:","z:");
    for (int i=0; i<numOfAtoms;){
		fprintf(file, "%4i%13.6lf%13.6lf%13.6lf\n",i,initialForces[3*i],initialForces[3*i+1],initialForces[3*i+2]);
		i++;
		}
		
    fclose(file);
	return;
	}
