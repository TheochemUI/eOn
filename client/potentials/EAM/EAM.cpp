    #include <math.h>
    #include <stdio.h>
    #include <string.h>
    #include <limits.h>
    #include <float.h> 
	#include <stdlib.h> 
	#include "EAM.h"


FILE *file;

int main();
void intialize();
void cleanmemory();
double force(long N, double *R, int *atomicNrs, double *F, double *U, const double *box);
double calc_force(long N, double *R, int *atomicNrs, double *F, double *U, const double *box);
void new_celllist(long N, const double *box, long *num_axis, long *cell_length, long *celllist_new, long num_cells, double *Rnew);
void cell_to_neighbor(long N, long num_of_cells, long *num_axis, long *cell_length, long *celllist_new, long *neigh_list);
int update_cell_list(long N, long num_cells,long *num_axis, long *cell_length, long *celllist_old, double *Rnew);//returns 0 if unchanged, >0 if changed - represents number of atoms that changed lists
double density (long N, long atom, double *R, int *atomicNrs, const double *box);//calculates local density of single atom

int main()
{
	//double Dm_min=0;
	//double Rm_min=0;
	//double en_min=2000;
	//double min_force=1000;
	
		initialize();
//	Rm=1.25;
   FILE *fforce=fopen("force_radius.xyz", "w");
	//for (Rm=.1;Rm<5.;Rm+=.005){
	//long N=2;
	//double R[]={0,0,0,6,0,0};
	long N=26;
	double R[26*3]={3.13796967E-06	,	3.1379696E-06	,	1.33342665952255	,
2.78708157269247	,	3.1379696E-06	,	1.33342665952256	,
5.57416000741527	,	3.1379696E-06	,	1.33342665952256	,
8.36123844213807	,	3.13796961E-06	,	1.33342665952256	,
11.1483168768609	,	3.13796961E-06	,	1.33342665952255	,
3.13796966E-06	,	2.78708157253175	,	1.33342665952255	,
2.78708157269247	,	2.78708157253175	,	1.33342665952255	,
5.57416000741528	,	2.78708157253175	,	1.33342665952255	,
8.36123844213808	,	2.78708157253175	,	1.33342665952255	,
11.1483168768609	,	2.78708157253175	,	1.33342665952254	,
3.13796966E-06	,	5.5741600070939	,	1.33342665952256	,
2.78708157269246	,	5.5741600070939	,	1.33342665952256	,
5.57416000741528	,	5.5741600070939	,	1.33342665952256	,
8.36123844213807	,	5.5741600070939	,	1.33342665952256	,
11.1483168768609	,	5.5741600070939	,	1.33342665952255	,
3.13796965E-06	,	8.36123844165604	,	1.33342665952257	,
2.78708157269246	,	8.36123844165604	,	1.33342665952256	,
5.57416000741527	,	8.36123844165604	,	1.33342665952255	,
8.36123844213808	,	8.36123844165605	,	1.33342665952256	,
11.1483168768609	,	8.36123844165604	,	1.33342665952256	,
3.13796966E-06	,	11.1483168762182	,	1.33342665952257	,
2.78708157269246	,	11.1483168762182	,	1.33342665952256	,
5.57416000741527	,	11.1483168762182	,	1.33342665952257	,
8.36123844213808	,	11.1483168762182	,	1.33342665952255	,
11.1483168768609	,	11.1483168762182	,	1.33342665952256	,
3.18062078872634	,	3.18062081747821	,	3.10568162746665};
	//double *R=new double[N*3];
	int *atomicNrs=new int[N];
	double box[3]={15,15,15};
	
	long k1=0;
	//long k2=0;
	//long k3=0;
	//long k4=0;
	//for (k1=0;k1<4;k1++)
		//for (k2=0;k2<4;k2++)
			//for (k3=0;k3<4;k3++)
				//{
				//R[k4]=k1;
				//R[k4+1]=k2;
				//R[k4+2]=k3;
				//k4+=3;
			//}
	for (k1=0;k1<N;k1++)
		atomicNrs[k1]=1;

	
	double *F=new double[3*N];
	double *U;
	
	char filename[20]="";
	sprintf(filename, "%ldatoms.xyz", N);
	
    file=fopen(filename, "w");
    



	double t=0;
	U=&t;
	
	double max_for=0;
	
	for (k1=0;k1<10000;k1++){
		*U=0;
	max_for=force(N, R, atomicNrs, F, U, box);
	//fprintf(fenergy, "%f %f %f\n", R[3], *U, F[3]);

	printf("run: %ld energy: %f\n", fcalled, *U);
	long ii=0;
	//R[3]-=.01;

	for (ii=25*3;ii<N*3;ii++)
		R[ii]+=F[ii];
	
	
}
printf("force on single: %.20e", sqrt(pow(F[75],2)+pow(F[76],2)+pow(F[77],2)));
printf("\nend location of lone atom: %f %f %f", R[75], R[76], R[77]);

	fprintf(fforce,"%f %.25f\n", Rm, sqrt(pow(F[75],2)+pow(F[76],2)+pow(F[77],2)));
	//if(max_for<min_force){
		//en_min=*U;
		//min_force=max_for;
		//Dm_min=Dm;
		//Rm_min=Rm;
	//}
		
	//delete R;
	delete atomicNrs;
	//delete box;
//}
//printf("Dm: %f, Rm: %f", Dm_min, Rm_min);
	return 1;
}

void cleanmemory()
{}

double force(long N, double *R, int *atomicNrs, double *F, double *U, const double *box)//calculate here long num_cells, long *num_axis, long *cell_length, //become global variables -long *celllist_old, long *celllist_new, long *neigh_list, long fcalled)
{	
	fcalled++;
long *num_axis=new long[3];
	long *cell_length=new long[3];
	//num_axis contains the number of cell lengths on each axis
	 long i=0;
	
	 for (i=0;i<3;i++){
	num_axis[i]=(long)(box[i]/rc[i])+1;
	}

	long num_cells;
	for (i=0;i<3;i++)
		cell_length[i]=(long)(box[i]/(num_axis[i]-1));
	//for (i=0;i<3;i++)
		//printf("celllength: %ld", cell_length[i]);

	num_cells=num_axis[0]*num_axis[1]*num_axis[2];
	
	if(fcalled==1)
	{
		celllist_old=new long[num_cells*(N+1)];
		celllist_new=new long[num_cells*(N+1)];
		neigh_list=new long[N*(N+1)];
	}
fprintf(file,"%ld\n\n", N);
long i3=0;
for (i3=0;i3<N;i3++)
    fprintf(file, "%s %f %f %f\n","He",R[3*i3], R[3*i3+1], R[3*i3+2]);
		
	*U=0;
	long kk=0;
	for (kk=0;kk<3*N;kk++)
		F[kk]=0;
	double max_frc_mag=0;
//	long neigh_list[N][N+1];


	//fcalled++;



   long xmin=LONG_MAX;
   long ymin=LONG_MAX;
   long zmin=LONG_MAX;


   //shift of the R values so that all coordinates are positive (makes creating cell table easier)
   double *Rtemp=new double[3*N];
   double *Rnew=new double[3*N];
   double *Rold=new double[3*N];
   double *boxtemp=new double[3];


   for (i=0;i<3*N;i+=3)
   {
	   Rtemp[i]=R[i];
	   Rtemp[i+1]=R[i+1];
	   Rtemp[i+2]=R[i+2];
	   if(R[i]<xmin)xmin=R[i];
	   if(R[i+1]<ymin)ymin=R[i+1];
	   if(R[i+2]<zmin)zmin=R[i+2];
	 }

   for (i=0;i<3*N;i+=3)
   {
	   Rtemp[i]+=abs(xmin);
	   Rtemp[i+1]+=abs(ymin);
	   Rtemp[i+2]+=abs(zmin);
   }
   boxtemp[0]=box[0];
   boxtemp[1]=box[1];
   boxtemp[2]=box[2];

//enforce periodic boundary conditions
for (i=0;i<3*N;i++)
    while(Rtemp[i]>boxtemp[i%3]) Rtemp[i]-=boxtemp[i%3];

	for (i=0;i<3*N;i++)
		Rnew[i]=Rtemp[i];
		
	if (fcalled==1)
	{
		new_celllist(N,boxtemp, num_axis, cell_length, celllist_new, num_cells, Rnew);
		cell_to_neighbor(N, num_cells, num_axis, cell_length, celllist_new, neigh_list);


	//long j;
	//printf("\n cell list new: \n");
	//for (j=0;j<num_cells;j++){
		//for (i=0;i<N+1;i++)
			//printf("%ld ", celllist_new[j*(N+1)+i]);
		//printf("\n");
	//}
	printf("\n\n end list");
	//for (i=0;i<N;i++){
		//for (j=0;j<N+1;j++)
			//printf("%ld ",neigh_list[i*(N+1)+j]);
		//printf("\n");
	//}
	
	}
	else
	{
		//for (i=0;i<3*N;i++)
			//Rold[i]=Rnew[i];

		if(update_cell_list(N, num_cells,num_axis, cell_length, celllist_old, Rnew)>0){
		new_celllist(N,boxtemp, num_axis, cell_length, celllist_new, num_cells, Rnew);	
			cell_to_neighbor(N, num_cells, num_axis, cell_length, celllist_new, neigh_list);
		
		
			//long j;
	//printf("\n updated cell list new: \n");
	//for (j=0;j<num_cells;j++){
		//printf("cell %ld: ", j);
		//for (i=0;i<N+1;i++)
			//printf("%ld ", celllist_new[j*(N+1)+i]);
		//printf("\n");
	//}
	//printf("\n\n");
}
}

max_frc_mag=calc_force(N, Rnew, atomicNrs, F, U, boxtemp);

for (i=0;i<num_cells*(N+1);i++)
	celllist_old[i]=celllist_new[i];

//double force_coef=.001;

for (i=0;i<3*N;i++){
   // Rnew[i]+=force_coef*F[i];
    Rold[i]=Rnew[i];

    if(i%3==0)    Rold[i]-=abs(xmin);
    else if (i%3==1)Rold[i]-=abs(ymin);
    else if(i%3==2)	Rold[i]-=abs(zmin);
    Rnew[i]=0;
    R[i]=Rold[i];
    }
// delete neigh_list;
 delete Rtemp;
 delete Rnew;
 delete Rold;
 delete boxtemp;
	return max_frc_mag;
 }
 
double calc_force(long N, double *R, int *atomicNrs, double *F, double *U, const double *box)
{
	//printf("%f %f %f %f", Dm, alphaM, Rm, betaR);
	long i=0;
    long j=0;
    long k=0;
    
    double max_frc_mag=0;

    for (i=0;i<3*N;i++) F[i]=0;

    *U=0;
    //double *vector_force= new double[3];
	double *vector_force=new double[3*N];
    for (i=0;i<N;i++)
    {

		long kk=0;
		for (kk=0;kk<3*N;kk++)
			vector_force[kk]=0;
			
		double dens=0;
		double add=0;
		double phi_r=0;

        for (j=0;j<neigh_list[i*(N+1)+N];j++)
        {
			long neigh= neigh_list[i*(N+1)+j];
			if(neigh>i){
	        double disx=R[3*i]-R[3*neigh];
            double disy=R[3*i+1]-R[3*neigh+1];
            double disz=R[3*i+2]-R[3*neigh+2];
			double dirs[]={disx, disy, disz};
			long u=0;
			for (u=0;u<3;u++)
				if(dirs[u]>box[u]/2) dirs[u]=box[u]-dirs[u];
				else if (dirs[u]< -box[u]/2)dirs[u]=-box[u]-dirs[u];

			double r=0;
			for (k=0;k<3;k++)
				r+= dirs[k]*dirs[k];
			r=sqrt(r);

			//if(r>rcut){ 
				//double slope=(rcut_energy-rend_energy)/(rcut-rend_cut);
			////	printf("\nslope: %f\n", slope);
				//phi_r+=slope*r+rcut_energy;
				//double comp = Dm*pow(1-exp(-alphaM*(r-Rm)),2)-Dm;
			////	printf("inside phi_r & compare: %f %f\n", phi_r, comp);
				////add=e_fit+phi_r*.5;
			//	//}
			//else{
         //   printf("disx: %f disy:%f disz: %f total: %f\n",dirs[0], dirs[1], dirs[2], r);
			  phi_r+= Dm*pow(1-exp(-alphaM*(r-Rm)),2)-Dm;//////////////change equations
            		//if(fcalled%100==0)
           // printf("phi_r: %f\n", phi_r);
     
			double curdens=0;
          curdens=density(N, i, R, atomicNrs, box);
          dens+=curdens;

           double mag_force = 2*alphaM*Dm*(exp(alphaM*r)-exp(alphaM*Rm))*(exp(alphaM*Rm-2*alphaM*r)); //old equation w/o density
			
			double mag_force_den=(537.7352*pow(curdens,7)-1771.224*pow(curdens,6)+2357.736*pow(curdens,5)-1640.015*pow(curdens,4)
								+663.052*pow(curdens,3)-179.4705*pow(curdens,2)+36.1594*curdens-2.00292);
		//	mag_force_den*=-pow(r,5)*((beta1*r-6)*exp(beta1*r)+1024*(beta2*r-3))*exp(-2*beta2*r);
			mag_force_den*=6*pow(r,5)*(512*exp(beta1*r)+exp(2*beta2*r))*exp(-beta1*r-2*beta2*r);
		//	printf("\nmag_force_den %f\n", mag_force_den);
			
			mag_force+=mag_force_den;	

			if(mag_force>max_frc_mag)max_frc_mag=mag_force;

            for (k=0;k<3;k++){
				vector_force[3*i+k]+=dirs[k]/r*-mag_force;
				vector_force[3*neigh+k]-=dirs[k]/r*-mag_force;
			}
        
        }
       }
        //if(dens>0)
        //printf("\ndensity: %f\n", dens);
            ////embedding function
            double density_after_func=0;
            long power=8;
            for (power=8;power>-0;power--)
				density_after_func+=func_param[8-power]*pow(dens, power);
			//if(dens>0)
			//printf("\ndensity after function: %f\n", density_after_func); 
            ///////////////////////
        add=phi_r+density_after_func;
			//printf("energy: %f\n phi: %f\n", add, phi_r);
 	         *U+=add;
 	    		//printf("\nforce vector atom %ld:", i);
        for (k=0;k<3*N;k++){
			double vec=vector_force[k];
			F[k]+=vec;
			//printf("%f ", vec);
		}
		//printf("\n");
		

		
    }
    delete vector_force;
    	//for (i=0;i<N*3;i++){
			//if(i%3==0)printf("\n %ld vector force: ", i/3);
						//printf(" %f", F[i]);
			//}


    return max_frc_mag;


}

 void new_celllist(long N, const double *box, long *num_axis, long *cell_length, long *celllist_new, long num_cells, double *Rnew)// creates a new cell list for the first time force is called
 {
	// printf("new celllist");
	long i=0;
	long j=0;

	////create a 2d array representing the cell list. the last column of each row tells now many atoms are in that cell
	//long cell_list[num_cells][N+1];//instantiates the cell_list.
	
	long *cell_list=new long[num_cells*(N+1)];
	for (i=0;i<num_cells;i++)
		for (j=0;j<N+1;j++)
			if (j%N==0&&j>0) cell_list[i*(N+1)+j]=0;
			else cell_list[i*(N+1)+j]=-1;

	//loop to put each atom into appropriate cell
	for (i=0;i<N;i++)
	{
		long *coor=new long[3];
		for (j=0;j<3;j++)
		coor[j]=(long)((Rnew[3*i+j])/cell_length[j]); //this calculates, for each axis, how many lengths away from 0,0,0 the atom is

		int cell=coor[0]*num_axis[1]*num_axis[2] + coor[1]*num_axis[2] + coor[2]; //this calculates which cell number the atom is in

		cell_list[cell*(N+1) + cell_list[cell*(N+1)+N]]=i;
		cell_list[cell*(N+1)+N]++;
		
		delete coor;
		}

		for (i=0;i<num_cells;i++)
			for (j=0;j<N+1;j++)
				celllist_new[i*(N+1)+j]=cell_list[i*(N+1)+j];
	delete cell_list;
}
void cell_to_neighbor(long N, long num_of_cells, long *num_axis, long *cell_length, long *celllist_new, long *neigh_list)
{
	num_of_cells=num_axis[0]*num_axis[1]*num_axis[2];
	long i2=0;

	long k=0;
	for (i2=0;i2<N;i2++)
		for (k=0;k<N+1;k++)
			neigh_list[i2*(N+1)+k]=-1;


	long i=0;
	long j=0,j1=0, j2=0;

	for (j=0;j<num_axis[0];j++)
		for (j1=0;j1<num_axis[1];j1++)
			for (j2=0;j2<num_axis[2];j2++)
				{
					long *neighbors=new long[N+1];//again, last index contains number of atoms in array.
					for(i=0;i<N+1;i++)
						neighbors[i]=-1;
					//neighbors[N]=0;
					//printf("%ld", neighbors[N]);
					long cur_index= j*num_axis[1]*num_axis[2]+ j1*num_axis[2]+j2;

						//	printf("%ld: ",cur_index);
					long d1=-1, d2=-1, d3=-1;

					//making a copy of cell list so that no cell is counted twice as a neighbor of a dif cell in small atom systems
					long *cell_list_copy=new long[num_of_cells*(N+1)+N+1];
					for (d1=0;d1<num_of_cells;d1++)
						for (d2=0;d2<N+1;d2++)
							cell_list_copy[d1*(N+1)+d2]=celllist_new[d1*(N+1)+d2];

						if(cell_list_copy[cur_index*(N+1)+N]==0)continue;

					for (d1=-1;d1<2;d1++)
						for (d2=-1;d2<2;d2++)
							for(d3=-1;d3<2;d3++)
							{
								long *pos=new long[3];
								pos[0]=j-d1; pos[1]=j1-d2; pos[2]=j2-d3;
								//check if pos is within bounds
								long y=0;
								for (y=0;y<3;y++)
									{
										if (pos[y]<0) pos[y]=num_axis[y]-1;
										else if (pos[y]>=num_axis[y]) pos[y]=0;
										}

								long neigh_index= pos[0]*num_axis[1]*num_axis[2]+ pos[1]*num_axis[2]+pos[2];

								long num=cell_list_copy[neigh_index*(N+1)+N];
								for (y=0;y<num;y++)
								{

									long check=0;
									k=0;

									for (k=0;k<N;k++)
										if(neighbors[check]==cell_list_copy[neigh_index*(N+1)+y]){check=-1;}
									if(check!=-1){
									neighbors[N]++;
									neighbors[neighbors[N]]=cell_list_copy[neigh_index*(N+1)+y];
								
									}
								}

							delete pos;
							}
						//	long y;
						//if(celllist_new[cur_index][N]>0){
						////	printf("\n%ld neighbors: ", cur_index);
							//for (y=0;y<N+1;y++)
								//printf("%ld ",neighbors[y]);
							//printf("\n");}

					for (i=0;i<celllist_new[cur_index*(N+1)+N];i++){
						long cur=celllist_new[cur_index*(N+1)+i];  //atom for which neighbors are being added to vlist
						long temp=0;
						neigh_list[cur*(N+1)+N]=0;
						for (temp=0;temp<=neighbors[N];temp++)
							if (cur!=neighbors[temp])
							{

								neigh_list[cur*(N+1)+ neigh_list[cur*(N+1)+N]]=neighbors[temp];
								neigh_list[cur*(N+1)+N]++;

							}

						}
						delete neighbors;
						delete cell_list_copy;
				}

}
int update_cell_list(long N, long num_cells,long *num_axis, long *cell_length,long *celllist_old, double *Rnew)//returns 0 if unchanged, >0 if changed - represents number of atoms that changed lists
{
	long i=0;
	long j=0;
	//long k=0;
	int changed=0;
	//long *tempcell=new long[num_cells*N+1];
	//celllist_new=tempcell; ////think about implememnting differently

	//for (i=0;i<num_cells;i++)
		//for (j=0;j<=N;j++)
			//celllist_new[i*(N+1)+j]=celllist_old[i*(N+1)+j];

	long *table=new long[N];//will hold the cell number for each atom in the old cell list

	for (i=0;i<num_cells;i++)
		for (j=0;j<celllist_old[i*(N+1)+N];j++)
			table[celllist_old[i*(N+1)+j]]=i;


	for (i=0;i<N;i++)
	{
		long *coor=new long[3];
		for (j=0;j<3;j++)
		coor[j]=(long)((Rnew[3*i+j])/cell_length[j]); //this calculates, for each axis, how many lengths away from 0,0,0 the atom is
		long cell=coor[0]*num_axis[1]*num_axis[2] + coor[1]*num_axis[2] + coor[2]; //this calculates which cell number the atom is in
		if(cell!=table[i])
		{
			changed++;
			//printf("\ncell new: %ld, cell old: %ld\n", cell, table[i]);
			//long oldcell=table[i];
			//for (k=0;k<celllist_new[oldcell*(N+1)+N];k++)
				//if(celllist_new[oldcell*(N+1)+k]==i)
					//{
						//celllist_new[oldcell*(N+1)+k]=celllist_new[oldcell*(N+1)+celllist_new[oldcell*(N+1)+N]];
						//celllist_new[oldcell*(N+1)+celllist_new[oldcell*(N+1)+N]]=-1;
						//celllist_new[oldcell*(N+1)+N]--;
					//}
		}
		delete coor;

	}
//	delete tempcell;
	delete table;
	return changed;

}
double density (long N, long atom, double *R, int *atomicNrs, const double *box)//calculates local density of single atom
{
	double D=0.0;



	long j=atom;
	long i=0;
	long k;
	for (i=0;i<neigh_list[j*(N+1)+N];i++)
	{
			long neigh=neigh_list[j*(N+1)+i];
			double disx=R[3*neigh]-R[3*j];
            double disy=R[3*neigh+1]-R[3*j+1];
            double disz=R[3*neigh+2]-R[3*j+2];

            double dirs[]={disx, disy, disz};
            for (k=0;k<3;k++)
				if(dirs[k]>box[k]/2) dirs[k]=box[k]-dirs[k];
				else if(dirs[k]<-box[k]/2) dirs[k]=box[k]+dirs[k];

			double r=0;
			for (k=0;k<3;k++)
				r+= dirs[k]*dirs[k];
			r=sqrt(r);
			if(r>r_cut) return 0;
			
			D+=pow(r,6)*(exp(-beta1 *r)+512*exp(-2*beta2*r));
			//printf("\ndensity: %f\n",D); 
	//		delete dirs;
	}
	return D;

}





