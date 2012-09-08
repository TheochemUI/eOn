#include "GlobalOptimizationJob.h"
#include "GlobalOptimization.h"
//#include "MinimizationJob.h"
#include "Dynamics.h"
#include "HelperFunctions.h"
#include "Log.h"
using namespace helper_functions;

static FILE *monfile;
static FILE *earrfile;
//****************************************************************************************
GlobalOptimizationJob::GlobalOptimizationJob(Parameters *params)
{
	parameters = params;
	nlmin=0;
	//etoler=parameters->globalOptimizationEtoler;
	//decisionMethod="NPEW";
	ediff=1.E-1;
	ekin=5.E-2;
	beta1=parameters->globalOptimizationBeta;
	beta2=parameters->globalOptimizationBeta;
	beta3=1.0/parameters->globalOptimizationBeta;
	alpha1=1.0/parameters->globalOptimizationAlpha;
	alpha2=parameters->globalOptimizationAlpha;
	mdmin=parameters->globalOptimizationMdmin;
	monfile=fopen("monitoring.dat","w");
	earrfile=fopen("earr.dat","w");
	firstStep=true;
	fcallsMove=0;
	fcallsRelax=0;
}
//****************************************************************************************
GlobalOptimizationJob::~GlobalOptimizationJob(void)
{
}
//****************************************************************************************
std::vector<std::string> GlobalOptimizationJob::run(void)
{
    //int status;
	GlobalOptimization globopt = GlobalOptimization(parameters);
    string reactant_passed = helper_functions::getRelevantFile(parameters->conFilename);
    std::vector<std::string> returnFiles;
    //returnFiles.push_back(reactant_output);
    Matter *matter_curr = new Matter(parameters);
    Matter *matter_hopp = new Matter(parameters);
    matter_curr->con2matter(reactant_passed);
    bool converged;
	long nstep = parameters->globalOptimizationSteps;
	AtomMatrix rat_t(matter_curr->numberOfAtoms(),3);
	//double epot_hopp;
	//std::vector<double> earr;
	//std::vector<Matter> allmatter;
    printf("\nBeginning minima hopping of %s\n", reactant_passed.c_str());
	//long fcalls;
	//printf("fcalls= %ld\n",matter_curr->getForceCalls());
	//printf("epot= %24.15E\n",matter_curr->getPotentialEnergy());
	//printf("fcalls= %ld\n",matter_curr->getForceCalls());
	//printf("epot= %24.15E\n",matter_curr->getPotentialEnergy());
	//printf("fcalls= %ld\n",matter_curr->getForceCalls());
    converged = matter_curr->relax(false, parameters->writeMovies, 
                                parameters->checkpoint, "min", "matter_curr");
    printf("converged %s \n",(converged)?"TRUE":"FALSE");
	//nlmin=0;
	//if(nlmin==0) 
	earr.push_back(matter_curr->getPotentialEnergy());
	matter_hopp[0]=matter_curr[0];
	GlobalOptimizationJob::report(matter_hopp);
	for(long istep=1;istep<=nstep;istep++) {
		//hoppingStep attempts to hopp into a new state followed by a minimization
		hoppingStep(istep,matter_curr,matter_hopp);
		//decisionStep decides to accept or reject this step if it escaped
		decisionStep(matter_curr,matter_hopp);
		//analyzing what happened in this step if it escaped
		analyze(matter_curr,matter_hopp);
		//reporting useful information about this hopp
		report(matter_hopp);
	}
	for(size_t i=0;i<earr.size();i++) {
		double earrim1;
		if(i==0) {
			earrim1=earr[0];
		}
		else {
			earrim1=earr[i-1];
			
		}
		fprintf(earrfile,"%5lu  %15.5f  %15.5f  %15.5f  \n",i+1,earr[i],earr[i]-earr[0],earr[i]-earrim1);
	}
	//globopt.run();
	fclose(monfile);
	fclose(earrfile);
	return returnFiles;
} //end of GlobalOptimizationJob::run
//****************************************************************************************
void GlobalOptimizationJob::analyze(Matter *matter_curr,Matter *matter_hopp) {
	if(escapeResult=="failure") {
		//matter_hopp[0]=matter_curr[0];
		hoppingResult="same";
		return;
	}
	double epot=matter_hopp->getPotentialEnergy();
	size_t jlo=hunt(epot);
	//printf("REZA: %lu \n",jlo);
	if(abs(epot-earr[jlo])<parameters->energyDifference) {
		hoppingResult="already_visited";
	}
	else {
		hoppingResult="new";
	}
	if(decisionResult=="accepted") {
		double epot_hopp=matter_hopp->getPotentialEnergy();
		matter_curr[0]=matter_hopp[0];
		size_t jlo=hunt(matter_hopp->getPotentialEnergy());
		if(hoppingResult=="new" && jlo==0 && epot_hopp<earr[0]) {
			printf("new lowest: nlmin,epot_hopp,dE %7d %15.5f %10.5f\n",1,
				epot_hopp,epot_hopp-earr[0]);
		}
		//printf("%15.5f  %15.5f  %15.5f  \n",epot_hopp,earr[0],earr[1]);
		//exit(1);
		insert(matter_curr);
	}
	else if(decisionResult=="rejected") {
		//matter_hopp[0]=matter_curr[0];
	}
	else {
		printf("ERROR: new minimum is neither accepted nor rejected: client stops.\n");
		exit(0);
	}
}
//****************************************************************************************
void GlobalOptimizationJob::examineEscape(Matter *matter_curr,Matter *matter_hopp) {
	//fprintf(monfile,"%15.5f  %15.5f  %15.5f  \n",matter_hopp->getPotentialEnergy(),
	//		matter_curr->getPotentialEnergy(),matter_hopp->getPotentialEnergy()-matter_curr->getPotentialEnergy());
	double epot, epot_hopp;
	epot=matter_curr->getPotentialEnergy();
	epot_hopp=matter_hopp->getPotentialEnergy();
	if(abs(epot_hopp-epot)<parameters->energyDifference) {
		escapeResult="failure";
	}
	else {
		escapeResult="success";
		//fprintf(monfile,"%15.5f  %15.5f  %15.5f           \n",epot_hopp,ediff,ekin);
	}
}
//****************************************************************************************
void GlobalOptimizationJob::applyMoveFeedbackMD(void) {
	if(firstStep) {
		firstStep=false;
		return;
	}
	if     (hoppingResult=="same")            {ekin*=beta1;}
	else if(hoppingResult=="already_visited") {ekin*=beta2;}
	else if(hoppingResult=="new")             {ekin*=beta3;}
	else {
		printf("ERROR: client does not know what to do with ekin. \n");
		printf("ERROR: client stops in applyMoveFeedbackMD.\n");
		exit(0);
	}
}
//****************************************************************************************
void GlobalOptimizationJob::applyDecisionFeedback(void) {
	if(decisionResult=="accepted") {ediff*=alpha1;}
	else {ediff*=alpha2;}
}
//****************************************************************************************
void GlobalOptimizationJob::report(Matter *matter_hopp) {
	char C1, C2;
	double ekin_p;
	if(hoppingResult=="same") {
		C1='S';
		C2='-';
		ekin_p=ekin*beta1;
	}
	else if(hoppingResult=="already_visited") { 
		C1='O';
		ekin_p=ekin*beta2;
	}
	else if(hoppingResult=="new") { 
		C1='N';
		ekin_p=ekin*beta3;
	}
	else {
		C1='-';
		ekin_p=ekin;
	}
	if(decisionResult=="accepted") {
		C2='A';
	}
	else if(decisionResult=="rejected") {
		C2='R';
	}
	else {
		C2='-';
	}
	double epot_hopp=matter_hopp->getPotentialEnergy();
	double temp = (2.0*ekin_p/8.6173857E-5);
	double dt=parameters->mdTimeStep/0.09823;
	fprintf(monfile,"%15.5f  %15.5f  %11lu  %12.2f         %c%c  %5ld  %5ld\n",
		epot_hopp,ediff,(size_t)temp,dt,C1,C2,fcallsMove,fcallsRelax);
}
//****************************************************************************************
void GlobalOptimizationJob::decisionStep(Matter *matter_curr,Matter *matter_hopp) {
	decisionResult="unknown";
	examineEscape(matter_curr,matter_hopp);
	if(escapeResult=="failure") {
		//matter_hopp[0]=matter_curr[0];
		return;
	}
	if(parameters->globalOptimizationDecisionMethod=="npew") {
		acceptRejectNPEW(matter_curr,matter_hopp);
		//GlobalOptimizationJob::update_minhopp_param(matter_hopp);
	}
	else if(parameters->globalOptimizationDecisionMethod=="boltzmann") {
		printf("WARNING: it will be implemented by Sam\n");
		exit(0);
	}
	else {
		printf("ERROR: accept/reject method not specified. client stops.\n");
		exit(0);
	}
	applyDecisionFeedback();
}
//****************************************************************************************
//NPEW: non-probablistic energy window (used in minima hopping method)
void GlobalOptimizationJob::acceptRejectNPEW(Matter *matter_curr,Matter *matter_hopp) {
	if(matter_hopp->getPotentialEnergy()<matter_curr->getPotentialEnergy()+ediff) {
		decisionResult="accepted";
	}
	else {
		decisionResult="rejected";
	}
}
//****************************************************************************************
void GlobalOptimizationJob::hoppingStep(long istep,Matter *matter_curr,Matter *matter_hopp) {
    bool converged;
	matter_hopp[0]=matter_curr[0];
	long fcalls1=matter_hopp->getForceCalls();
	if(parameters->globalOptimizationMoveMethod=="md") {
		applyMoveFeedbackMD();
		mdescape(matter_hopp);
	}
	else if(parameters->globalOptimizationMoveMethod=="random") {
		printf("WARNING: it will be implemented by Sam\n");
		exit(0);
	}
	long fcalls2=matter_hopp->getForceCalls();
	hoppingResult="unknown";
    converged = matter_hopp->relax(false, parameters->writeMovies, 
		parameters->checkpoint, "min", "matter_hopp");
    printf("converged %s \n",(converged)?"TRUE":"FALSE");
	long fcalls3=matter_hopp->getForceCalls();
	fcallsMove=fcalls2-fcalls1;
	fcallsRelax=fcalls3-fcalls2;
}
//****************************************************************************************
void GlobalOptimizationJob::mdescape(Matter *matter)
{
	int nmd;
	double ekinc, epot, etot, epot0, etot0;
	Dynamics *dyn = new Dynamics(matter,parameters);
	velopt(matter);
	epot=matter->getPotentialEnergy();
	ekinc=matter->getKineticEnergy();
	etot=ekinc+epot;
	epot0=epot;
	etot0=etot;
	size_t nummax=0, nummin=0;
	double enmin1=0.0, enmin2=0.0, en0000=0.0, econs_max=-1.E100, econs_min=1.E100, devcon;
	bool md_presumably_escaped=false;
    printf("MD  %5d  %20.10E  %15.5E  %15.5E  \n",0,epot-epot0,ekinc,etot-etot0);
	nmd=1000;
	for(int imd=1;imd<=nmd;imd++) {
		enmin2=enmin1;
		enmin1=en0000;
		dyn->velocityVerlet();
		epot=matter->getPotentialEnergy();
		ekinc=matter->getKineticEnergy();
		etot=ekinc+epot;
		en0000=epot-epot0;
		if(enmin1>enmin2 && enmin1>en0000)  nummax=nummax+1;
		if(enmin1<enmin2 && enmin1<en0000)  nummin=nummin+1;
    	printf("MD  %5d  %15.5f  %15.5f  %12.2E  %4lu  %4lu\n",imd,epot-epot0,ekinc,etot-etot0,nummax,nummin);
		econs_max=max(econs_max,ekinc+epot);
		econs_min=min(econs_min,ekinc+epot);
		if(nummin>=(size_t) mdmin) {
			if(nummax!=nummin) printf("WARNING: iproc,nummin,nummax %4lu %4lu",nummin,nummax);
			md_presumably_escaped=true;
			break;
		}
	} //end of loop over imd
	devcon=econs_max-econs_min;
	if(md_presumably_escaped) {
		devcon=devcon/(double)(matter->numberOfFreeAtoms()*3);
		if(devcon/ekin<2.E-3) {
			parameters->mdTimeStep*=1.1;
		}
		else {
			parameters->mdTimeStep/=1.1;
		}
	}
	else {
		printf("TOO MANY MD STEPS  \n");
		parameters->mdTimeStep*=2.0;
	}
}
//****************************************************************************************
void GlobalOptimizationJob::velopt(Matter *matter)
{
	AtomMatrix vat(matter->numberOfAtoms(),3);
	double tt1, tt2, tt3, vtot[3]; //, ekin_t;
	int iat;
	//matter->numberOfAtoms();
	vtot[0]=0.0; vtot[1]=0.0; vtot[2]=0.0;
	for(iat=0;iat<matter->numberOfAtoms();iat++) {
		tt1=randomDouble();
		tt2=randomDouble();
		tt3=randomDouble();
		vat(iat,0)=(tt1-0.5)*2.0;
		vat(iat,1)=(tt2-0.5)*2.0;
		vat(iat,2)=(tt3-0.5)*2.0;
		vtot[0]+=vat(iat,0);
		vtot[1]+=vat(iat,1);
		vtot[2]+=vat(iat,2);
	}
    printf("Linear momentum  %15.5E  %15.5E  %15.5E  \n",vtot[0],vtot[1],vtot[2]);
	vtot[0]/=matter->numberOfAtoms();
	vtot[1]/=matter->numberOfAtoms();
	vtot[2]/=matter->numberOfAtoms();
	for(iat=0;iat<matter->numberOfAtoms();iat++) {
		vat(iat,0)-=vtot[0];
		vat(iat,1)-=vtot[1];
		vat(iat,2)-=vtot[2];
	}
	matter->setVelocities(vat);
	long nFreeCoords = matter->numberOfFreeAtoms()*3;
	double kinE = matter->getKineticEnergy();
	double kinT = (2.0*kinE/nFreeCoords/8.6173857E-5);
	double temperature = (2.0*ekin/8.6173857E-5);
	matter->setVelocities(vat*sqrt(temperature/kinT));
}
//****************************************************************************************
//void ::rescaleVelocity()
//{
//    AtomMatrix velocity = matter->getVelocities();
//    double kinE = matter->getKineticEnergy();
//    double kinT = (2.0*kinE/nFreeCoords/kb);
//    matter->setVelocities(velocity*sqrt(temperature/kinT));
//}
//****************************************************************************************
void GlobalOptimizationJob::insert(Matter *matter)
{
	double epot;
	//vector<double> epot_hopp;
	size_t jlo, jlo_insert;
	//vector<double>::iterator it;
	epot=matter->getPotentialEnergy();
	jlo=hunt(epot);
    printf("JLO= %5lu  %10.5f  \n",jlo,abs(epot-earr[jlo]));
	//it=earr.begin()+jlo;
	//epot_hopp.push_back(epot);
	if(!(abs(epot-earr[jlo])<parameters->energyDifference)) {
		//earr.insert(it,epot_hopp.begin(),epot_hopp.end());
		jlo_insert=jlo;
		if(epot>earr[jlo]) jlo_insert++;
		earr.insert(earr.begin()+jlo_insert,1,epot);
	}
}
//****************************************************************************************
size_t GlobalOptimizationJob::hunt(double epot) {
    //epot is in interval [earr(jlo),earr(jlo+1)[ ; earr(0)=-Infinity ; earr(n+1) = Infinity
    size_t jlo;
	double de;
	for(jlo=0;jlo<earr.size();jlo++) if(epot<earr[jlo]) break;
    if(jlo==earr.size()) jlo--;
	de=abs(epot-earr[jlo]);
    if(jlo>0) if(abs(epot-earr[jlo-1])<de) jlo--; //{jlo--;de=abs(epot-earr[jlo]);}
    //if(jlo!=earr.size()-1) if(abs(epot-earr[jlo+1])<parameters->energyDifference) jlo++;
	return jlo;
}
//****************************************************************************************
//****************************************************************************************
//****************************************************************************************
//****************************************************************************************
