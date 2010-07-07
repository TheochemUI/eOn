#include <math.h> 
#include <iostream>

//#include "../../system_unit.h" // unit converters
#include "PotentialsInterface.h"

//	Variables
	long *celllist_old;
	long *celllist_new;
	long *neigh_list;
	long fcalled;
	double Dm, alphaM, Rm, beta1,beta2;
	double r_cut;
	double rc[3]={5,5,5}; //2 is arbitrary number. rc represents the optimal length for each cell in cell list
    double func_param[9]={67.2169,-253.032,392.956,-328.003,165.763,-59.8235,18.0797,-2.00292,-.0102076};
    double u0;
    double cuttOffR;
    double psi;
    
    double cuttOffU;
    
    // Functions
	// constructor
	
    // Just to satify interface
    void initialize() {
		fcalled=0;
		Dm=1;
		alphaM=3.0205380362464;
		Rm=2.65;
		beta1=6.6137657075868;
		beta2= 6.0000000000000;
		
		

		r_cut=5.5;
		}
    void cleanMemory(){
		delete celllist_old;
		delete celllist_new;
		delete neigh_list;   
	} 
    
    void force(long N, const double *R, const long *atomicNrs, double *F, double *U, const double *box);
    void setParameters(double r0Recieved, double u0Recieved, double psiRecieved);


