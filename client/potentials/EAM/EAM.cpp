#include <math.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <float.h>
#include <stdlib.h>
#include "EAM.h"


EAM::EAM()
{
    return;
}

void EAM::initialize()
{
    initialized = false;
    Dm = 1;
    alphaM = 3.0205380362464;
    Rm = 2.65;
    beta1 = 6.6137657075868;
    beta2 = 6.0000000000000;
    r_cut = 5.5;
    rc = new double[3]; 
    rc[0] = rc[1] = rc[2] = 5.0;// 5 is arbitrary number. rc represents the optimal length for each cell in cell list
    double temp[] = {67.2169,-253.032,392.956,-328.003,165.763,-59.8235,18.0797,-2.00292,-.0102076};
    func_param = new double[9];
    for(int i = 0; i < 9; i++){
        func_param[i] = temp[i];
    }
}

void EAM::cleanMemory()
{
    delete celllist_old;
    delete celllist_new;
    delete neigh_list;
}


// Calculate here long num_cells, long *num_axis, long *cell_length, //become global variables -long *celllist_old, long *celllist_new, long *neigh_list, long fcalled)
void EAM::force(long N, const double *R, const long *atomicNrs, double *F, double *U, const double *box)
{
    long *num_axis=new long[3];
    long *cell_length=new long[3];
    //num_axis contains the number of cell lengths on each axis
    long i=0;

    for (i=0;i<3;i++)
    {
        num_axis[i]=(long)(box[i]/rc[i])+1;
    }

    long num_cells;
    for (i=0;i<3;i++)
        cell_length[i]=(long)(box[i]/(num_axis[i]-1));
    //for (i=0;i<3;i++)
    //printf("celllength: %ld", cell_length[i]);

    num_cells=num_axis[0]*num_axis[1]*num_axis[2];

    if(!initialized)
    {
        celllist_old=new long[num_cells*(N+1)];
        celllist_new=new long[num_cells*(N+1)];
        neigh_list=new long[N*(N+1)];
    }

    *U=0;
    long kk=0;
    for (kk=0;kk<3*N;kk++)
        F[kk]=0;
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

    if (!initialized)
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
//        printf("\n\n end list");
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

        if(update_cell_list(N, num_cells,num_axis, cell_length, celllist_old, Rnew)>0)
        {
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

    calc_force(N, Rnew, atomicNrs, F, U, boxtemp);

    for (i=0;i<num_cells*(N+1);i++)
        celllist_old[i]=celllist_new[i];

    //double force_coef=.001;

//    for (i=0;i<3*N;i++)
//    {
//        // Rnew[i]+=force_coef*F[i];
//        Rold[i]=Rnew[i];

//        if(i%3==0)    Rold[i]-=abs(xmin);
//        else if (i%3==1)Rold[i]-=abs(ymin);
//        else if(i%3==2) Rold[i]-=abs(zmin);
//        Rnew[i]=0;
//        R[i]=Rold[i];
//    }
    // delete neigh_list;
    delete Rtemp;
    delete Rnew;
    delete Rold;
    delete boxtemp;
    initialized = true;
}


void EAM::calc_force(long N, double *R, const long *atomicNrs, double *F, double *U, double *box)
{
    //printf("%f %f %f %f", Dm, alphaM, Rm, betaR);
    long i=0;
    long j=0;
    long k=0;

//    double max_frc_mag=0;

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
            if(neigh>i)
            {
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
                                 //////////////change equations
                phi_r+= Dm*pow(1-exp(-alphaM*(r-Rm)),2)-Dm;
                //if(fcalled%100==0)
                // printf("phi_r: %f\n", phi_r);

                double curdens=0;
                curdens=density(N, i, R, atomicNrs, box);
                dens+=curdens;

                                 //old equation w/o density
                double mag_force = 2*alphaM*Dm*(exp(alphaM*r)-exp(alphaM*Rm))*(exp(alphaM*Rm-2*alphaM*r));

                double mag_force_den=(537.7352*pow(curdens,7)-1771.224*pow(curdens,6)+2357.736*pow(curdens,5)-1640.015*pow(curdens,4)
                    +663.052*pow(curdens,3)-179.4705*pow(curdens,2)+36.1594*curdens-2.00292);
                //	mag_force_den*=-pow(r,5)*((beta1*r-6)*exp(beta1*r)+1024*(beta2*r-3))*exp(-2*beta2*r);
                mag_force_den*=6*pow(r,5)*(512*exp(beta1*r)+exp(2*beta2*r))*exp(-beta1*r-2*beta2*r);
                //	printf("\nmag_force_den %f\n", mag_force_den);

                mag_force+=mag_force_den;

//                if(mag_force>max_frc_mag)max_frc_mag=mag_force;

                for (k=0;k<3;k++)
                {
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
        for (k=0;k<3*N;k++)
        {
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

}


// Creates a new cell list for the first time force is called
void EAM::new_celllist(long N, const double *box, long *num_axis, long *cell_length, long *celllist_new, long num_cells, double *Rnew)
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
                                 //this calculates, for each axis, how many lengths away from 0,0,0 the atom is
            coor[j]=(long)((Rnew[3*i+j])/cell_length[j]);

                                 //this calculates which cell number the atom is in
        int cell=coor[0]*num_axis[1]*num_axis[2] + coor[1]*num_axis[2] + coor[2];

        cell_list[cell*(N+1) + cell_list[cell*(N+1)+N]]=i;
        cell_list[cell*(N+1)+N]++;

        delete coor;
    }

    for (i=0;i<num_cells;i++)
        for (j=0;j<N+1;j++)
            celllist_new[i*(N+1)+j]=cell_list[i*(N+1)+j];
    delete cell_list;
}


void EAM::cell_to_neighbor(long N, long num_of_cells, long *num_axis, long *cell_length, long *celllist_new, long *neigh_list)
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
                                 //again, last index contains number of atoms in array.
                long *neighbors=new long[N+1];
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
                                    if(check!=-1)
                                {
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

                for (i=0;i<celllist_new[cur_index*(N+1)+N];i++)
                {
                                 //atom for which neighbors are being added to vlist
                    long cur=celllist_new[cur_index*(N+1)+i];
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


// Returns 0 if unchanged, >0 if changed - represents number of atoms that changed lists
int EAM::update_cell_list(long N, long num_cells,long *num_axis, long *cell_length,long *celllist_old, double *Rnew)
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

    long *table=new long[N];     //will hold the cell number for each atom in the old cell list

    for (i=0;i<num_cells;i++)
        for (j=0;j<celllist_old[i*(N+1)+N];j++)
            table[celllist_old[i*(N+1)+j]]=i;

    for (i=0;i<N;i++)
    {
        long *coor=new long[3];
        for (j=0;j<3;j++)
                                 //this calculates, for each axis, how many lengths away from 0,0,0 the atom is
            coor[j]=(long)((Rnew[3*i+j])/cell_length[j]);
                                 //this calculates which cell number the atom is in
        long cell=coor[0]*num_axis[1]*num_axis[2] + coor[1]*num_axis[2] + coor[2];
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


// Calculates local density of single atom
double EAM::density (long N, long atom, double *R, const long *atomicNrs, const double *box)
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
