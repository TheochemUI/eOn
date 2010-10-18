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

  long i, j,count;
  bool flag  = 1;
  
  nBAs = 3; //nBoostAtoms
  nRAs = nAtoms - nBAs;//nRestAtoms
  nTABs = nBAs*(nBAs-1)/2+nBAs*nRAs; //number of Bonds involved with Tagged Atoms
  printf("nTABs=%ld\n",nTABs);
  count = 0;
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


  for(i=0 ;i< nBAs;i++){
      printf("nBoostAtoms %ld: %ld\n",i,BAList[i]);
  }
  for(i=0; i< nRAs;i++){
      printf("nRestAtoms %ld: %ld\n",i,RAList[i]);
  }
}

void BondBoost::boost(){
    long nReg = 0,RMDS = parameters->mdSteps;
    Matrix<double, Eigen::Dynamic, 1> TABL_tmp(nTABs,1);
    bool flag = 0;
    
    if(nReg <= RMDS){
        flag = 0;
    }
    else{ 
        flag = 1; 
    }
    
    if( flag == 0 ){
        TABL_tmp=Rmdsteps();
        TABLList = TABLList + (1.0/RMDS)*TABL_tmp;

        for(long i=0;i<nTABs;i++){
              printf("Distance beweem Atoms %ld and %ld is %lf\n",TABAList[2*i],TABAList[2*i+1],TABL_tmp(i,0));
        }
        
        printf("\n");
        for(long i=0;i<nTABs;i++){
              printf("Distance beweem Atoms %ld and %ld is %lf\n",TABAList[2*i],TABAList[2*i+1],TABLList(i,0));
        }

       nReg ++ ;   
    }
  printf("hello,let's boost!\n");
}


Matrix<double, Eigen::Dynamic,1> BondBoost::Rmdsteps() {
     Matrix<double, Eigen::Dynamic, 1> TABL_t(nTABs,1);
     long count = 0, i, j;

     for(i=0;i<nBAs;i++){
         for(j=i+1;j<nBAs;j++){
             printf("distance = %lf\n", matter->distance(BAList[i],BAList[j]));
             TABL_t(count,0) = matter->distance(BAList[i],BAList[j]);
             TABAList[2*count] = BAList[i];
             TABAList[2*count+1] = BAList[j];
             count ++;
         }
     }

     for(i=0;i<nBAs;i++){
         for(j=0;j<nRAs;j++){
             printf("distance = %lf\n", matter->distance(BAList[i],RAList[j]));
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
     

