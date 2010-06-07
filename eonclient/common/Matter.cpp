//Jean-Claude C. Berthet, University of Iceland, 2006
/*
 *===============================================
 *  Modified. Name, Date and a small description!
 *
 *-----------------------------------------------
 *  Todo:
 *
 *===============================================
 */
#include <cmath>
#include<cstdlib>
#include<cstdio>
#include<cstring>
#include<string>
#include<iostream>
#include <cassert>
#include "Matter.h"

using namespace std;
using namespace system_unit;
using namespace constants;

namespace {
       const char *elementArray[] = {"Unknown", "H","He","Li","Be","B","C","N","O",
           "F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc",
           "Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se",
           "Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag",
           "Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
           "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta",
           "W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
           "Fr","Ra","Ac","Th","Pa","U", NULL};
      /// guess the atom type from the atomic mass,
        std::string mass2atom(double atomicmass) {
            return elementArray[int(atomicmass+.5)];
        }
      const int MAXC=10;/// Maximum number of components for functions matter2con and con2matter
      int symbol2atomicNumber(char const * symbol)
      {
            int i=0;

            while (elementArray[i] != NULL) {
                if (strcmp(symbol, elementArray[i]) == 0) {
                    return i;
                }
                i++;
            }
      }
      
      char const * atomicNumber2symbol(int n)
      {
            return elementArray[n];
      }
}

Matter::Matter(Parameters *parameters)
{
    initialiseDataMembers(parameters);
}

Matter::Matter(Parameters *parameters, const long int nAtoms)
{
    resize(nAtoms);// prepare memory for nAtoms
    initialiseDataMembers(parameters);
}

Matter::Matter(const Matter& matter)
{
    operator=(matter);
}

Matter::~Matter()
{
    clearMemory();
}

const Matter& Matter::operator=(const Matter& matter)
{
      resize(matter.numberOfAtoms());
      nAtoms_ = matter.nAtoms_;
      long i;
      for(i=0; i<3*nAtoms_; i++) {
            positions_[i] = matter.positions_[i];
            forces_[i] = matter.forces_[i];
      };
      if(matter.velocities_) {
            for(i=0; i<3*nAtoms_; i++)
                  velocities_[i] = matter.velocities_[i];
      };
      for(i=0; i<nAtoms_; i++) {
            masses_[i] = matter.masses_[i];
            atomicNrs_[i] = matter.atomicNrs_[i];
            isFixed_[i] = matter.isFixed_[i];
      };
      cellBoundaries_[0] = matter.cellBoundaries_[0];
      cellBoundaries_[1] = matter.cellBoundaries_[1];
      cellBoundaries_[2] = matter.cellBoundaries_[2];
      
      constraints_ = matter.constraints_;
      usePeriodicBoundaries_ = matter.usePeriodicBoundaries_;
      
      potentialEnergy_ = matter.potentialEnergy_;
      computePotential_ = matter.computePotential_;
      
      strcpy(headerCon1_,matter.headerCon1_);
      strcpy(headerCon2_,matter.headerCon2_);
      strcpy(headerCon4_,matter.headerCon4_);
      strcpy(headerCon5_,matter.headerCon5_);
      strcpy(headerCon6_,matter.headerCon6_);
      
      return *this;
}

bool Matter::operator==(const Matter& matter) {//To compare two matter objects, if all differences in positions are bellow getMaxDifferencePos, they are considered equal
    bool result = false;
    long nCoord, i = 0;
    double *pos, diffR, diffRX, diffRY, diffRZ;
    
    nCoord = 3*nAtoms_;
    pos = new double[nCoord];   
    
    if(matter.numberOfAtoms()==nAtoms_){
        // Note that the state of result is changed 
        // in order to be able to break out of the loop 
        // when the first large difference appear
        result = true;
        matter.getPositions(pos);
        
        for(i=0;i<nAtoms_;i++){
            diffRX = positions_[ 3*i ]-pos[ 3*i ];
            diffRY = positions_[3*i+1]-pos[3*i+1];
            diffRZ = positions_[3*i+2]-pos[3*i+2];
            // floor = largest integer value less than argument
            diffRX=diffRX-cellBoundaries_[0]*floor(diffRX/cellBoundaries_[0]+0.5);  
            diffRY=diffRY-cellBoundaries_[1]*floor(diffRY/cellBoundaries_[1]+0.5);
            diffRZ=diffRZ-cellBoundaries_[2]*floor(diffRZ/cellBoundaries_[2]+0.5);
            
            diffR = sqrt(diffRX*diffRX+diffRY*diffRY+diffRZ*diffRZ);
            
            if(getMaxDifferencePos()<diffR){                
                result = false;
                break;
            }
        }
    }
    delete [] pos;
    return result;
}



double Matter::distanceTo(const Matter& matter) 
{
    /* RT: Returns the distance to the given matter object. */
    
    double *pos;
    pos = new double[3 * nAtoms_];   
    
    long i = 0;
    double dX, dY, dZ;
    double sum = 0.0;
    
    if(matter.numberOfAtoms() == nAtoms_)
    {
        matter.getPositions(pos);
        for(i = 0; i < nAtoms_; i++)
        {
            dX = positions_[3 * i + 0] - pos[3 * i + 0];
            dY = positions_[3 * i + 1] - pos[3 * i + 1];
            dZ = positions_[3 * i + 2] - pos[3 * i + 2];
            dX = dX - cellBoundaries_[0] * floor(dX / cellBoundaries_[0] + 0.5);  
            dY = dY - cellBoundaries_[1] * floor(dY / cellBoundaries_[1] + 0.5);
            dZ = dZ - cellBoundaries_[2] * floor(dZ / cellBoundaries_[2] + 0.5);
            sum += dX * dX + dY * dY + dZ * dZ;
        }
    }
    delete [] pos;
    return sqrt(sum);
}



double Matter::per_atom_norm(const Matter& matter) 
{
    /* RT: Returns the maximum distance between two atoms in the Matter objects. */
    
    double *pos;
    pos = new double[3 * nAtoms_];   
    
    long i = 0;
    double dX, dY, dZ;
    double max_distance = 0.0;
    double distance = 0.0;
    
    if(matter.numberOfAtoms() == nAtoms_)
    {
        matter.getPositions(pos);
        for(i = 0; i < nAtoms_; i++)
        {
            dX = positions_[3 * i + 0] - pos[3 * i + 0];
            dY = positions_[3 * i + 1] - pos[3 * i + 1];
            dZ = positions_[3 * i + 2] - pos[3 * i + 2];
            dX = dX - cellBoundaries_[0] * floor(dX / cellBoundaries_[0] + 0.5);  
            dY = dY - cellBoundaries_[1] * floor(dY / cellBoundaries_[1] + 0.5);
            dZ = dZ - cellBoundaries_[2] * floor(dZ / cellBoundaries_[2] + 0.5);
            distance = sqrt(dX * dX + dY * dY + dZ * dZ);
            if(distance > max_distance)
            {
                max_distance = distance;
            }
        }
    }
    delete [] pos;
    return max_distance;
}



void Matter::resize(const long int nAtoms)
{
      clearMemory();
      potential_ = new Potentials(parameters_);
      
      if(nAtoms>0) {
            nAtoms_ = nAtoms;
            positions_ = new double[3*nAtoms];
            velocities_ = new double[3*nAtoms];
            forces_ = new double[3*nAtoms];
            masses_ = new double[nAtoms];
            atomicNrs_ = new long[nAtoms];
            isFixed_ = new int[nAtoms];
            for(long i=0; i<nAtoms; i++) {
                  positions_[i*3]=0.0;  positions_[i*3+1]=0.0;  positions_[i*3+2]=0.0;
                  velocities_[i*3]=0.0;   velocities_[i*3+1]=0.0;   velocities_[i*3+2]=0.0;
                  forces_[i*3]=0.0;     forces_[i*3+1]=0.0;     forces_[i*3+2]=0.0;            
                  masses_[i] = 0.0;
                  atomicNrs_[i] = 0;
                  isFixed_[i] = false;
            };
      };
}

long int Matter::numberOfAtoms() const {return(nAtoms_);}// return the number of atoms

double Matter::getBoundary(int axis) const
{
    return(cellBoundaries_[axis]);
}

void Matter::setBoundary(int axis, double length)
{
    cellBoundaries_[axis]=length;
    if(usePeriodicBoundaries_){
        applyPeriodicBoundary(axis);
    }
    computePotential_=true;
}

void Matter::activatePeriodicBoundaries()
{
    usePeriodicBoundaries_=true;
    applyPeriodicBoundary(0);
    applyPeriodicBoundary(1);
    applyPeriodicBoundary(2);
}

void Matter::desactivatePeriodicBoundaries()
{
    usePeriodicBoundaries_=false;
}

void Matter::setConstraints(Constraints constraints)
{
    constraints_=constraints;
}

double Matter::getPosition(long int indexAtom, int axis) const {
    return positions_[3*indexAtom+axis];
}

void Matter::setPosition(long int indexAtom, int axis, double position)
{
    positions_[3*indexAtom+axis]=position;
    if(usePeriodicBoundaries_){ 
        applyPeriodicBoundary(indexAtom, axis);
    }
    computePotential_=true;
}

void Matter::getPositions(double pos[]) const {//return coordinates of free atoms in array 'pos'
    for(long int i=0; i<nAtoms_; i++) {
        pos[ 3*i ]=positions_[ 3*i ];
        pos[3*i+1]=positions_[3*i+1];
        pos[3*i+2]=positions_[3*i+2];
    };
}

void Matter::setPositions(const double pos[]) {//Update Matter with the new positions of the free atoms given in array 'pos'
    for(long int i=0; i<nAtoms_; i++) {
        positions_[ 3*i ]=pos[ 3*i ];
        positions_[3*i+1]=pos[3*i+1];
        positions_[3*i+2]=pos[3*i+2];
    };    
    if(usePeriodicBoundaries_){
        applyPeriodicBoundary();
    }
    computePotential_=true;
}

void Matter::getForces(double forces[]) const {// return forces applied on all atoms in array 'force' 
    computePotential();
    for(long int i=0; i<nAtoms_; i++) {
        if(!isFixed_[i]) {
            forces[ 3*i ]=forces_[ 3*i ];
            forces[3*i+1]=forces_[3*i+1];
            forces[3*i+2]=forces_[3*i+2];
        }
        else{
            forces[ 3*i ]=0;
            forces[3*i+1]=0;
            forces[3*i+2]=0;
        };
    };
}

double Matter::distance(long index1, long index2) const{// return distance between the atoms with index1 and index2
    double diffRX, diffRY, diffRZ;
    diffRX = positions_[ 3*index1 ]-positions_[ 3*index2 ];
    diffRY = positions_[3*index1+1]-positions_[3*index2+1];
    diffRZ = positions_[3*index1+2]-positions_[3*index2+2];
    
    // floor = largest integer value less than argument
    diffRX = diffRX-cellBoundaries_[0]*floor(diffRX/cellBoundaries_[0]+0.5);  
    diffRY = diffRY-cellBoundaries_[1]*floor(diffRY/cellBoundaries_[1]+0.5);
    diffRZ = diffRZ-cellBoundaries_[2]*floor(diffRZ/cellBoundaries_[2]+0.5);
    
    return(sqrt(diffRX*diffRX+diffRY*diffRY+diffRZ*diffRZ));
}

double Matter::distance(const Matter& matter, long index) const{// return the distance atom with index has moved between the current Matter object and the Matter object passed as argument
    double diffRX, diffRY, diffRZ;
    double *pos;
    pos = new double[3*nAtoms_];
    
    matter.getPositions(pos);
    
    diffRX = positions_[ 3*index ]-pos[ 3*index ];
    diffRY = positions_[3*index+1]-pos[3*index+1];
    diffRZ = positions_[3*index+2]-pos[3*index+2];
    
    // floor = largest integer value less than argument
    diffRX = diffRX-cellBoundaries_[0]*floor(diffRX/cellBoundaries_[0]+0.5);  
    diffRY = diffRY-cellBoundaries_[1]*floor(diffRY/cellBoundaries_[1]+0.5);
    diffRZ = diffRZ-cellBoundaries_[2]*floor(diffRZ/cellBoundaries_[2]+0.5);
    
    delete [] pos;
    
    return(sqrt(diffRX*diffRX+diffRY*diffRY+diffRZ*diffRZ));
}

double Matter::getMass(long int indexAtom) const
{
    return(masses_[indexAtom]);
}

void Matter::setMass(long int indexAtom, double mass)
{
    masses_[indexAtom]=mass;
    computePotential_=true;
}

long Matter::getAtomicNr(long int indexAtom) const
{
    return(atomicNrs_[indexAtom]);
}

void Matter::setAtomicNr(long int indexAtom, long atomicNr)
{
    atomicNrs_[indexAtom]=atomicNr;
    computePotential_=true;
}

int Matter::getFixed(long int indexAtom) const {
    return(isFixed_[indexAtom]);
}

void Matter::setFixed(long int indexAtom, int isFixed)
{
    isFixed_[indexAtom]=isFixed;
}

double Matter::potentialEnergy() const
{
      if(nAtoms_>0) {
            computePotential();
            return potentialEnergy_;
      } 
      else 
            return 0.0;
}

double Matter::kineticEnergy() const
{
    double K=0;
    if(velocities_) {
        long int n3Atoms=3*nAtoms_;
        for(long int i=0; i<n3Atoms; i++) {
            if(!isFixed_[i/3]) K+=velocities_[i]*velocities_[i]*masses_[i/3]*0.5;
        };
    };
    return K;
}

double Matter::mechanicalEnergy() const
{
    return potentialEnergy()+kineticEnergy();
}

long int Matter::numberOfFreeAtoms() const {
    long int nfree=0;
    for(long int i=0; i< nAtoms_; i++) {
            if(!isFixed_[i]) nfree++;//count the number of free atoms
    };
    return(nfree);
}

void Matter::getFreePositions(double pos[]) const {//returns coordinates of free atoms in array 'pos'
    long int j=0;
    for(long int i=0; i<nAtoms_; i++) {
        if(!isFixed_[i]) {
            pos[3*j  ]=positions_[3*i  ];
            pos[3*j+1]=positions_[3*i+1];
            pos[3*j+2]=positions_[3*i+2];
            j++;
        };
    };
}

void Matter::setFreePositions(const double pos[]) {//Update Matter with the new positions of the free atoms given in array 'pos'
    long int j=0;
    for(long int i=0; i<nAtoms_; i++) {
        if(!isFixed_[i]) {
            assert(isnormal(pos[3*j]));
            positions_[3*i  ]=pos[3*j  ];
            assert(isnormal(pos[3*j+1]));
            positions_[3*i+1]=pos[3*j+1];
            assert(isnormal(pos[3*j+2]));
            positions_[3*i+2]=pos[3*j+2];
            j++;
        };
    };
    if(usePeriodicBoundaries_){
            applyPeriodicBoundary();
    }
    computePotential_=true;
}

bool Matter::getFreeVelocities(double velocities[]) const{
    long int j=0;
    if(velocities_) {
        for(long int i=0; i<nAtoms_; i++) {
            if(!isFixed_[i]) {
                velocities[3*j  ]=velocities_[3*i  ];
                velocities[3*j+1]=velocities_[3*i+1];
                velocities[3*j+2]=velocities_[3*i+2];
                j++;
            };
        };
        return true;
    } else {
        return false;
    };
}

void Matter::setFreeVelocities(const double velocities[]) {
    long int j=0;
    if(!velocities_) {
        velocities_=new double[3*nAtoms_];
        for(long int i=0; i<nAtoms_; i++) {
            if(!isFixed_[i]) {
                velocities_[  3*i   ]=velocities[  3*j  ];
                velocities_[3*i+1]=velocities[3*j+1];
                velocities_[3*i+2]=velocities[3*j+2];
                j++;
            } else {
                velocities_[ 3*i ]=0;
                velocities_[3*i+1]=0;
                velocities_[3*i+2]=0;
            };
        };  
    } else {
        for(long int i=0; i<nAtoms_; i++) {
            if(!isFixed_[i]) {
                velocities_[ 3*i ]=velocities[ 3*j ];
                velocities_[3*i+1]=velocities[3*j+1];
                velocities_[3*i+2]=velocities[3*j+2];
                j++;
            };
        };
    };
}

void Matter::getFreeAccelerations(double accelerations[]) const {
    long int j=0;
    computePotential();
    for(long int i=0; i<nAtoms_; i++) {
        if(!isFixed_[i]) {
            accelerations[3*j  ]=forces_[3*i  ]/masses_[i];
            accelerations[3*j+1]=forces_[3*i+1]/masses_[i];
            accelerations[3*j+2]=forces_[3*i+2]/masses_[i];
            j++;
        };
    };
}

void Matter::getFreeForces(double forces[]) const {// return forces applied on free atoms in array 'force' 
    long int j=0;
    computePotential();
    for(long int i=0; i<nAtoms_; i++) {
        if(!isFixed_[i]) {
            forces[3*j  ]=forces_[3*i  ];
            forces[3*j+1]=forces_[3*i+1];
            forces[3*j+2]=forces_[3*i+2];
            j++;
        };
    };
}

void Matter::updateForces(double positions[], double velocities[], double forces[]) {
    if(constraints_) {// When constraints are applied.
        //The old coordinates are stored in positionBefore_ (for reference). As constraint algorithms requires to old coordinates in order to apply to correction.
        // see members constraint() & UpdateCons()
        if(positionsBefore_==0) positionsBefore_=new double [3*nAtoms_];
        for(int i=0; i<3*nAtoms_; i++) positionsBefore_[i]=positions_[i];
    };
    if(positions) setFreePositions(positions);
    if(velocities) setFreeVelocities(velocities);
    applyConstraints();
    if(positions) getFreePositions(positions);
    if(velocities) getFreeVelocities(velocities);
    if(forces) getFreeForces(forces);
}

void Matter::updateAccelerations(double positions[], double velocities[],double accelerations[]) {
    if(constraints_) {// When constraints are applied.
        //The old coordinates are stored in positionBefore_ (for reference). As constraint algorithms requires to old coordinates in order to apply to correction.
        // see members constraint() & UpdateCons()
        if(positionsBefore_==0) positionsBefore_=new double [3*nAtoms_];
        for(int i=0; i<3*nAtoms_; i++) positionsBefore_[i]=positions_[i];
    };
    if(positions) setFreePositions(positions);
    if(velocities) setFreeVelocities(velocities);
    applyConstraints();
    if(positions) getFreePositions(positions);
    if(velocities) getFreeVelocities(velocities);
    if(accelerations) getFreeAccelerations(accelerations);
}

void Matter::applyConstraints() {
    if(constraints_) {
        // The 2nd and 3rd arguments are normally the position at time t and t+dt respectively. positions_ contains the position at time t+dt. This feature is not yet fully implemented.
        constraints_(nAtoms_, positionsBefore_, positions_, velocities_, cellBoundaries_);
    };
}

void Matter::getFreeMasses(double m[]) const {
    long int j=0;
    for(long int i=0; i<nAtoms_; i++) {
        if(!isFixed_[i]) {
            m[j]=masses_[i];
            j++;
        };
    };
}

long Matter::getForceCalls() const{
    return(forceCalls_);
}

void Matter::resetForceCalls(){
    forceCalls_ = 0;
    return;
}

//Print atomic coordinate to a .xyz file
void Matter::matter2xyz(std::string filename, bool append /*Append if file already exist*/) const {
    FILE * file;
    long int i;
    filename+=".xyz";
    #ifdef _FILESYS_
        std::string resfile;
        boinc_resolve_filename_s(filename.c_str(), resfile);
        if (append) file=boinc_fopen(resfile.c_str(),"a");
        else file=boinc_fopen(resfile.c_str(),"w");
    #else
        if (append) file=fopen(filename.c_str(),"a");
        else file=fopen(filename.c_str(),"w");
    #endif
    if(file==0) {
        cerr << "Can't create file " << filename << endl;
        exit(1);
    };
    fprintf(file,"%ld\nGenerated by EON\n", numberOfAtoms());
    for(i=0;i<numberOfAtoms();i++) {
        fprintf(file,"%s\t%f\t%f\t%f\n", atomicNumber2symbol(getAtomicNr(i)), getPosition(i, 0)/ANGSTROM, getPosition(i, 1)/ANGSTROM, getPosition(i, 2)/ANGSTROM);
    };
    fclose(file);
}

// Print atomic coordinates to a .con file
bool Matter::matter2con(std::string filename) const {
    bool state;
    FILE *file;
    filename += ".con";
    file = fopen(filename.c_str(),"w");     
    state = matter2con(file);
    fclose(file); 
    return(state);
}

bool Matter::matter2con(FILE *file) const
{
      long int i;
      int j;
      long int Nfix=0; // Nfix to store the number of fixed atoms
      int Ncomponent=0; // Used to store the number of components (eg water: two components H and O)
      int first[MAXC]; // To store the position of the first atom of each component plus at the end the total number of atoms
      double mass[MAXC];
      long atomicNrs[MAXC];
      first[0]=0;
      if(numberOfAtoms()>0) {
            if(getFixed(0)) Nfix=1;//count the number of fixed atoms
            mass[0]=getMass(0);
            atomicNrs[0]=getAtomicNr(0);
      };
      j=0;
      for(i=1;i<numberOfAtoms();i++) {
            if(getFixed(i)) Nfix++;//count the number of fixed atoms
            if(getAtomicNr(i) != atomicNrs[j]) {// Check if there is a second component
                  j++;
                  if(j>=MAXC) {
                        std::cerr << "Does not support more than " << MAXC << " components and the atoms must be ordered by component.\n";
                        return false;
                  };
                  mass[j]=getMass(i);
                  atomicNrs[j]=getAtomicNr(i);
                  first[j]=i;
            };
      };
      first[j+1]=numberOfAtoms();
      Ncomponent=j+1;
      
      fputs(headerCon1_, file);
      fputs(headerCon2_, file);
      fprintf(file, "%f\t%f\t%f\n", getBoundary(0)/ANGSTROM, getBoundary(1)/ANGSTROM, getBoundary(2)/ANGSTROM);
      fputs(headerCon4_, file);
      fputs(headerCon5_, file);
      fputs(headerCon6_, file);
      
      fprintf(file, "%d\n", Ncomponent);
      for(j=0; j<Ncomponent; j++) {
            fprintf(file, "%d ", first[j+1]-first[j]);
      }
      fprintf(file, "\n");  
      for(j=0; j<Ncomponent; j++) {
            mass[j]/=G_PER_MOL;
            fprintf(file, "%f ", mass[j]);
      };
      fprintf(file, "\n");
      for(j=0; j<Ncomponent; j++) {
            fprintf(file, "%s\n", atomicNumber2symbol(atomicNrs[j]));
            fprintf(file, "Coordinates of Component %d\n", j+1);
            for(i=first[j]; i<first[j+1]; i++) {
                  fprintf(file,"%.19f\t%.19f\t%.19f\t%d\t%ld\n", getPosition(i, 0)/ANGSTROM, getPosition(i, 1)/ANGSTROM, getPosition(i, 2)/ANGSTROM, getFixed(i), i);
            };
      };
      return true;
}

// Load atomic coordinates from a .con file
bool Matter::con2matter(std::string filename) {
    bool state;
    FILE *file;
    // Add the .con extension to filename if it is not already there.
    int pos=filename.find_last_of('.');
    if(filename.compare(pos+1, 3, "con")){
        filename+=".con";
    };
    file=fopen(filename.c_str(), "r");
    if (!file) {
        cerr << "File " << filename << " was not found.\n";
        return(false);
    };
    
    state = con2matter(file);
    fclose(file);
    return(state);
}
    
bool Matter::con2matter(FILE *file) {
      char line[255];// Temporary string of character to read from the file.
      fgets(headerCon1_,sizeof(line),file);
      if (strchr(headerCon1_,'\r')) {
            /* Files created on Windows or on Mac with Excell have carriage returns (\r) instead of or along
            with the new line charater (\n). C recognises only the \n as the end of line. */
            cerr << "A carriage return ('\\r') has been detected. To work correctly, new lines should be indicated by the new line character (\\n).";
            return false;// return false for error
      };
      
      long int i; int j;
      
      fgets(headerCon2_,sizeof(line),file);
      
      fgets(line,sizeof(line),file);
      
      double x, y, z;
      sscanf(line,"%lf %lf %lf", &x, &y, &z);// The third line contains the length of the periodic cell
      setBoundary(0, x * ANGSTROM);
      setBoundary(1, y * ANGSTROM);
      setBoundary(2, z * ANGSTROM);
      
      fgets(headerCon4_,sizeof(line),file);
      sscanf(headerCon4_,"%lf %lf %lf", &x, &y, &z);// The fourth line contains the angles of the cell vectors
      if ( (x != 90.0) or (y != 90.0) or (z != 90.0) ) {
            /* The code only supports cubic simulation cells*/
            cerr << "This code only supports cubic cells.";
            return false;// return false for error
      };
      
      fgets(headerCon5_,sizeof(line),file);
      fgets(headerCon6_,sizeof(line),file);
      
      fgets(line,sizeof(line),file);
      int Ncomponent;// Number of components is the number of of different types of atoms. For instance H2O (water) has two component (H and O).
      if(sscanf(line,"%d",&Ncomponent)==0) {
            std::cout << "The number of components cannot be read. One component is assumed instead\n";
            Ncomponent=1;
      };
      if((Ncomponent>MAXC)||(Ncomponent<1)) {
            cerr << "con2atoms does not support more than " << MAXC << " components (or less than 1).\n";
            return false;
      };
      /* to store the position of the 
            first atom of each element 'MAXC+1': the last element is used to store the total number of atom.*/ 
      long int first[MAXC+1];
      long int Natoms=0;
      first[0]=0;
      // Now we want to know the number of atom of each type. Ex with H2O, two hydrogens and one oxygen
      for(j=0; j<Ncomponent; j++) {
            fscanf(file, "%ld", &Natoms);
            first[j+1]=Natoms+first[j];
      };    
      fgets(line, sizeof(line), file);// Discard the rest of the line
      resize(first[Ncomponent]);// Set the total number of atoms, and allocates memory
      double mass[MAXC];
      for(j=0; j<Ncomponent; j++) {// Now we want to know the number of atom of each type. Ex with H2O, two hydrogens and one oxygen
            fscanf(file, "%lf", &mass[j]);
            mass[j]*=G_PER_MOL;// conversion of g/mol to local units. (see su.h)
      };
      fgets(line,sizeof(line),file);//Discard rest of the line
      int atomicNr;
      int fixed;
      for (j=0; j<Ncomponent; j++) {
            char symbol[3];
            fgets(line,sizeof(line),file);
            sscanf(line, "%2s\n", symbol);
            atomicNr=symbol2atomicNumber(symbol);
            fgets(line,sizeof(line),file);// skip one line
            for (i=first[j]; i<first[j+1]; i++){
                  setMass(i, mass[j]);
                  setAtomicNr(i, atomicNr);
                  fgets(line, sizeof(line), file);
                  sscanf(line,"%lf %lf %lf %d\n", &x, &y, &z, &fixed);
                  setPosition(i, 0, x * ANGSTROM);
                  setPosition(i, 1, y * ANGSTROM);
                  setPosition(i, 2, z * ANGSTROM);
                  setFixed(i, static_cast<bool>(fixed));
            };
      };
      if(usePeriodicBoundaries_){ 
            applyPeriodicBoundary();// Transform the coordinate to use the minimum image convention.
      }
      //    potential_ = new Potentials(parameters_);
      return(true);
}

void Matter::computePotential() const
{
      if(computePotential_) {
            if(potential_) {
                  assert(isnormal(positions_[0]));
                  potential_->force(nAtoms_, positions_, atomicNrs_, forces_, &potentialEnergy_, cellBoundaries_);
                  forceCalls_ = forceCalls_+1;
                  computePotential_=false;
            }
            else {
                  cerr << "No potential associated with the atomic structure." << endl;
                  exit(1);
            };
      };
}

void Matter::initialiseDataMembers(Parameters *parameters)
{
      nAtoms_ = 0;
      positions_ = 0;
      positionsBefore_ = 0;
      velocities_ = 0;
      forces_ = 0;
      masses_ = 0;
      atomicNrs_ = 0;
      isFixed_ = 0;
      constraints_ = 0;
      cellBoundaries_[0] = 0.0;
      cellBoundaries_[1] = 0.0;
      cellBoundaries_[2] = 0.0; 
      usePeriodicBoundaries_ = true;
      computePotential_ = true;
      forceCalls_ = 0;

      parameters_ = parameters;
      potential_ = new Potentials(parameters_);
}

void Matter::clearMemory()
{
      // the pointer to parameters should not be deleted, it is a reference
      // to shared data
      //delete parameters_;
      
      if (positions_!=0) {
            delete [] positions_;
            positions_=0;
      };
      if (positionsBefore_!=0) {
            delete [] positionsBefore_;
            positionsBefore_=0;
      };
      if (velocities_!=0) {
            delete [] velocities_;
            velocities_=0;
      };
      if (forces_!=0) {
            delete [] forces_;
            forces_=0;
      };
      if (masses_!=0) {
            delete [] masses_;
            masses_=0;
      };
      if (atomicNrs_!=0) {
            delete [] atomicNrs_;
            atomicNrs_=0;
      };
      if(isFixed_!=0) {
            delete [] isFixed_;
            isFixed_=0;
      };
      if (potential_!=0){
            delete potential_;
            potential_=0;
      }
}

void Matter::applyPeriodicBoundary() // Transform the coordinate to use the minimum image convention.
{
    applyPeriodicBoundary(0);
    applyPeriodicBoundary(1);
    applyPeriodicBoundary(2);
}

void Matter::applyPeriodicBoundary(int axis) {
    for(long int i=0; i<nAtoms_; i++)
        applyPeriodicBoundary(i, axis);
}

void Matter::applyPeriodicBoundary(long atom, int axis)
{
    while(positions_[atom*3+axis]>cellBoundaries_[axis]) {
        positions_[atom*3+axis]-=cellBoundaries_[axis];
    };
    while(positions_[atom*3+axis]<= 0.0) {
        positions_[atom*3+axis]+=cellBoundaries_[axis];
    };
}
