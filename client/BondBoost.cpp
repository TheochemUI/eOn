#include "BondBoost.h"
#include <math.h>

BondBoost::BondBoost(Matter *matt,Parameters *params)
{
    matter = matt;    
    parameters = params;
    nAtoms = matter->numberOfAtoms();

};


BondBoost::~BondBoost()
{
    return;
};


void BondBoost::initial(){

  long i, j,count = 0;
  bool flag  = 1;
  
  nBBs = 0;
  nReg = 0;
  nBAs = 3; //nBoostAtoms
  nRAs = nAtoms - nBAs;//nRestAtoms
  nTABs = nBAs*(nBAs-1)/2+nBAs*nRAs; //number of Bonds involved with Tagged Atoms
  BAList = new long[nBAs];//BoostAtomsList
  RAList = new long[nRAs];//RestAtomsList
  TABAList = new long[2*nTABs];//CorrespondingAtomsList of TABLList
  TABLList.setZero(nTABs,1);
  
  for(i = 0;i < nBAs; i++){
      BAList[i]=i;
  }

  for(i = 0 ;i < nAtoms; i++){
     flag = 1;
     for (j = 0; j < nBAs; j++){
         if(i == BAList[j]){
             flag = 0;
         }
     }
     if(flag == 1){
         RAList[count] = i;   
         count ++;
     }
  } 
  
  if(count != nRAs){
      printf("Error: nRestAtoms does not equal to counted number!\n");
  }

// TEST PRINT
/*
  for(i=0 ;i< nBAs;i++){
      printf("nBoostAtoms %ld: %ld\n",i,BAList[i]);
  }
  for(i=0; i< nRAs;i++){
      printf("nRestAtoms %ld: %ld\n",i,RAList[i]);
  }
*/  
}

void BondBoost::boost(){

    long RMDS = parameters->mdSteps-3;
    Matrix<double, Eigen::Dynamic, 1> TABL_tmp(nTABs,1);
    bool flag = 0;
    
    if(nReg < RMDS){
        flag = 0;
    }
    else{ 
        flag = 1; 
    }
    
    if( flag == 0 ){
        TABL_tmp=Rmdsteps();
        TABLList = TABLList + (1.0/RMDS)*TABL_tmp;
        nReg ++;
        // TEST PRINT
        for(long i=0;i<nTABs;i++){
              printf("Distance beweem Atoms %ld and %ld is %lf\n",TABAList[2*i],TABAList[2*i+1],TABL_tmp(i,0));
        }
        
        printf("\n");
        for(long i=0;i<nTABs;i++){
              printf("Distance beweem Atoms %ld and %ld is %lf\n",TABAList[2*i],TABAList[2*i+1],TABLList(i,0));
        }
        // TEST END
    }
    else{ 
           if( nReg == RMDS+1){
               printf("First steps = RMDS+1\n");
               nBBs = BondSelect();
               for (long i=0;i<nBBs;i++){
                   printf("Equilibrium Distance between Atoms %ld and %ld is %lf\n",BBAList[2*i],BBAList[2*i+1],EBBLList(i,0));
               }
           }

           Epsr_Q = new double[nBBs];           
           CBBLList.setZero(nBBs,1);
           Booststeps();
           nReg ++;
    }

       
  printf("hello,let's boost! nReg = %ld\n",nReg);
}

void BondBoost::Booststeps(){
     long i,j,Mi;
     long AtomI_1,AtomI_2;
     double QRR, PRR, Epsr_MAX, A_EPS_M, Sum_V, Boost_Fact, DVMAX;
     double Dforce, Fact_1, Fact_2,Mforce;
     double Ri[3] = {0.0} , R = 0.0;

     Matrix<double, Eigen::Dynamic, 3> OriForce(nAtoms,3);
     Matrix<double, Eigen::Dynamic, 3> Free(nAtoms,3);
     Matrix<double, Eigen::Dynamic, 3> AddForce(nBBs,3);
     Matrix<double, Eigen::Dynamic, 3> TADF(nAtoms,3);
     Matrix<double, Eigen::Dynamic, 3> NewForce(nAtoms,3);
 
     OriForce.setZero();
     TADF.setZero();
     Free.setZero();
     AddForce.setZero();
     NewForce.setZero();

     QRR = 0.3;
     PRR = 0.95;
     Epsr_MAX = 0.0;
     A_EPS_M = 0.0;
     Boost_Fact = 0.0;
     Sum_V = 0.0;
     DVMAX = 1.0;
     Dforce = 0.0;
     Mforce = 0.0;
     Fact_1 = 0.0;
     Fact_2 = 0.0;

     printf("haha::BoostSteps\n");

     for(i=0;i<nBBs;i++){
         AtomI_1 = BBAList[2*i];
         AtomI_2 = BBAList[2*i+1];
         CBBLList(i,0) = matter->distance(AtomI_1,AtomI_2);
     }

     for (i=0;i<nBBs;i++){
         Epsr_Q[i]=(CBBLList(i,0)-EBBLList(i,0))/EBBLList(i,0)/QRR;
         if (abs(Epsr_Q[i]) >= Epsr_MAX ){
             Epsr_MAX = abs(Epsr_Q[i]);
         }
     }
     
     A_EPS_M=(1.0-Epsr_MAX*Epsr_MAX)*(1.0-Epsr_MAX*Epsr_MAX)/(1.0-PRR*PRR*Epsr_MAX*Epsr_MAX);
/*
     for (i=0;i<nBBs;i++){
         printf("Boost::Distance between Atoms %ld and %ld is %lf; EPSR_Q=%lf\n",BBAList[2*i],BBAList[2*i+1],CBBLList(i,0),Epsr_Q[i]);
     }

     printf("Boost::Epsr_MAX= %lf\n",Epsr_MAX);
     printf("Boost::A_EPS_M= %lf\n",A_EPS_M);
*/   
     if(Epsr_MAX < 1.0){
        for(i=0;i<nBBs;i++){
	    Sum_V +=DVMAX*(1.0-Epsr_Q[i]*Epsr_Q[i])/double(nBBs);
	}
	Boost_Fact = A_EPS_M*Sum_V;
	//SPtime
 	//SDtime
	//Ave_Boost_Fact
   
        for(i=0;i<nBBs;i++){
	    if( abs(Epsr_Q[i]) < Epsr_MAX ){
		Fact_1 = 2.0*A_EPS_M*DVMAX*Epsr_Q[i]/QRR/EBBLList(i,0)/double(nBBs);
		Dforce = Fact_1;
             
		AtomI_1 = BBAList[2*i];
		AtomI_2 = BBAList[2*i+1];
		
 		R=matter->distance(AtomI_1,AtomI_2);

		for(j=0;j<3;j++){
                    Ri[j] = matter->pdistance(AtomI_1,AtomI_2,j);   
	            AddForce(i,j)=Ri[j]/R*Dforce;
                    TADF(AtomI_1,j)=TADF(AtomI_1,j)+AddForce(i,j);
                    TADF(AtomI_2,j)=TADF(AtomI_2,j)-AddForce(i,j);
		}
              //  printf("%ld Bond:: Rx=%lf;Ry=%lf;Rz=%lf;Rsum=%lf;R=%lf\n",i,Ri[0],Ri[1],Ri[2],sqrt(Ri[0]*Ri[0]+Ri[1]*Ri[1]+Ri[2]*Ri[2]),R);

	    }	
            else if( abs(Epsr_Q[i]) == Epsr_MAX ){
  		Fact_1 = 2.0*A_EPS_M*DVMAX*Epsr_Q[i]/QRR/EBBLList(i,0)/double(nBBs);
 		Fact_2 = 2.0*(1.0-(Epsr_Q[i])*Epsr_Q[i])*Epsr_Q[i]*(2.0*(1.0-PRR*PRR*Epsr_Q[i]*Epsr_Q[i])-PRR*PRR*(1.0-(Epsr_Q[i])*Epsr_Q[i]))/QRR/EBBLList(i,0)/(1.0-PRR*PRR*Epsr_Q[i]*Epsr_Q[i])/(1.0-PRR*PRR*Epsr_Q[i]*Epsr_Q[i]); 
        	Dforce=Fact_1+Sum_V*Fact_2;

		AtomI_1 = BBAList[2*i];
                AtomI_2 = BBAList[2*i+1];

                Mi=i;
		Mforce=Dforce;
                R=matter->distance(AtomI_1,AtomI_2);

                for(j=0;j<3;j++){
                    Ri[j] = matter->pdistance(AtomI_1,AtomI_2,j);
                    AddForce(i,j)=Ri[j]/R*Dforce;
                    TADF(AtomI_1,j)=TADF(AtomI_1,j)+AddForce(i,j);
                    TADF(AtomI_2,j)=TADF(AtomI_2,j)-AddForce(i,j);
		}
	    }
/*
            if(i != Mi){
	        printf("DeltaF(Epsr_Q = %lf ), %lf  %lf  %lf; Atom %ld and Atom %ld\n", Epsr_Q[i],AddForce(i,0),AddForce(i,1),AddForce(i,2),BBAList[2*i],BBAList[2*i+1]);
	    }else{
		printf("DeltaF(Epsr_MAX = %lf) = %lf, %lf  %lf  %lf; Atom %ld and Atom %ld\n", Epsr_MAX,Mforce,AddForce(i,0),AddForce(i,1),AddForce(i,2),BBAList[2*i],BBAList[2*i+1]);
	   }
*/
	}
         OriForce = matter->getForces();
        NewForce = OriForce + TADF;
        Free = matter->getFree();
        NewForce = NewForce.cwise() *  Free;
	matter->setForces(NewForce);
       
        printf("TADF::\n");
	for(i=0;i<nAtoms;i++){
		printf("%lf   %lf   %lf\n",TADF(i,0),TADF(i,1),TADF(i,2));
        }
 
        printf("OriForces::\n");
        for(i=0;i<nAtoms;i++){
                printf("%lf   %lf   %lf\n",OriForce(i,0),OriForce(i,1),OriForce(i,2));
        }
  
        printf("NewForces::\n");
        for(i=0;i<nAtoms;i++){
                printf("%lf   %lf   %lf\n",NewForce(i,0),NewForce(i,1),NewForce(i,2));
        }
        OriForce = matter->getForces();
        printf("NewSettedForces::\n");
        for(i=0;i<nAtoms;i++){
                printf("%lf   %lf   %lf\n",OriForce(i,0),OriForce(i,1),OriForce(i,2));
        }

       



     }
     else if(Epsr_MAX >= 1.0){
   	Boost_Fact = 0.0;
	//SPtime = 
	//Ave_Boost_Fact =
	//
     }
}


Matrix<double, Eigen::Dynamic,1> BondBoost::Rmdsteps() {
     Matrix<double, Eigen::Dynamic, 1> TABL_t(nTABs,1);
     long count = 0, i, j;

     for(i=0;i<nBAs;i++){
         for(j=i+1;j<nBAs;j++){
             //printf("distance = %lf\n", matter->distance(BAList[i],BAList[j]));
             TABL_t(count,0) = matter->distance(BAList[i],BAList[j]);
             TABAList[2*count] = BAList[i];
             TABAList[2*count+1] = BAList[j];
             count ++;
         }
     }

     for(i=0;i<nBAs;i++){
         for(j=0;j<nRAs;j++){
             //printf("distance = %lf\n", matter->distance(BAList[i],RAList[j]));
             TABL_t(count,0) = matter->distance(BAList[i],RAList[j]);
             TABAList[2*count] = BAList[i];
             TABAList[2*count+1] = RAList[j];
             count ++;
         }
     }
//     printf("count=%ld\n",count);
     if(count != nTABs){
         printf("Total Involved Bond number does not equal counted number\n");
     }

//     for(i=0;i<count;i++){
//         printf("Distance beweem Atoms %ld and %ld is %lf\n",TABAList[2*i],TABAList[2*i+1],TABL_t(i,0));
//    }
     return TABL_t;
}

long BondBoost::BondSelect(){
    long count = 0,i,nBBs_tmp = 0;
    double Qcutoff = 3.0;

    for(i=0;i<nTABs;i++){
        if(TABLList(i,0) <= Qcutoff){
            count ++;
        }
    }
    nBBs_tmp = count;
    EBBLList.setZero(nBBs_tmp,1);
    BBAList = new long[2*nBBs_tmp];
    count = 0;
    for(i=0;i<nTABs;i++){
        if(TABLList(i,0) <= Qcutoff){
           EBBLList(count,0) = TABLList(i,0);
           BBAList[2*count] = TABAList[2*i];
           BBAList[2*count+1] = TABAList[2*i+1];
           count ++;
        }
    }
    return nBBs_tmp;
}



