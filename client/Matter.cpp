//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include "Matter.h"
#include "Constants.h"

#include <cmath>
#include<cstdlib>
#include <cassert>


using namespace std;

namespace {

    const char *elementArray[] = {"Unknown", "H","He","Li","Be","B","C","N","O",
           "F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca","Sc",
           "Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se",
           "Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag",
           "Cd","In","Sn","Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
           "Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu","Hf","Ta",
           "W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
           "Fr","Ra","Ac","Th","Pa","U", NULL};

    // guess the atom type from the atomic mass,
    std::string mass2atom(double atomicmass) {
        return elementArray[int(atomicmass+.5)];
    }

    const int MAXC=10; // Maximum number of components for functions matter2con and con2matter

    int symbol2atomicNumber(char const * symbol)
    {
        int i=0;

        while (elementArray[i] != NULL) {
            if (strcmp(symbol, elementArray[i]) == 0) {
                return i;
            }
            i++;
        }
        // invalid symbol
        return -1;
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
    resize(nAtoms); // prepare memory for nAtoms
    initialiseDataMembers(parameters);
}

void Matter::initialiseDataMembers(Parameters *params)
{
    nAtoms = 0;
    cellBoundaries.resize(3,3);
    cellBoundaries.setZero();
    usePeriodicBoundaries = true;
    recomputePotential = true;
    forceCalls = 0;
    nsteps = 0;
    parameters = params;
    potential = new Potentials(parameters);
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
    nAtoms = matter.nAtoms;
    
    positions = matter.positions;
    forces = matter.forces;
    masses = matter.masses;
    atomicNrs = matter.atomicNrs;
    isFixed = matter.isFixed;
    cellBoundaries = matter.cellBoundaries;
    velocities = matter.velocities; 
    
    usePeriodicBoundaries = matter.usePeriodicBoundaries;
 
    potentialEnergy = matter.potentialEnergy;
    recomputePotential = matter.recomputePotential;
 
    strcpy(headerCon1,matter.headerCon1);
    strcpy(headerCon2,matter.headerCon2);
    strcpy(headerCon4,matter.headerCon4);
    strcpy(headerCon5,matter.headerCon5);
    strcpy(headerCon6,matter.headerCon6);

    // liang add here nsteps for test:
    nsteps = matter.nsteps;

    return *this;
}

// Two matter objects are considered the same if all 
// differences in positions are bellow getMaxDifferencePos.
bool Matter::operator==(const Matter& matter) {
    return (parameters->maxDifferencePos) > perAtomNorm(matter);
}

double Matter::distanceTo(const Matter& matter) 
{
    /* RT: Returns the distance to the given matter object. */

    return pbc(positions - matter.positions).norm();
}

Matrix<double, Eigen::Dynamic, 3> Matter::pbc(Matrix<double, Eigen::Dynamic, 3> diff) const
{
    Matrix<double, 3, 3> ibox = cellBoundaries.inverse();
    Matrix<double, Eigen::Dynamic, 3> ddiff = diff*ibox;
   
    int i,j;
    for(i=0; i<diff.rows(); i++)
    {
        for(j=0; j<3; j++)
        {
            ddiff(i,j) = fmod(fmod(ddiff(i,j), 1.0)  + 1.5, 1.0) -.5;
        }
    }

    return ddiff*cellBoundaries;
}

double Matter::perAtomNorm(const Matter& matter) 
{
    /* RT: Returns the maximum distance between two atoms in the Matter objects. */

    long i = 0;
    double max_distance = 0.0;

    if(matter.numberOfAtoms() == nAtoms)
    {
        Matrix<double, Eigen::Dynamic, 3> diff = pbc(positions - matter.positions);
        for(i = 0; i < nAtoms; i++)
        {
            max_distance = max(diff.row(i).norm(), max_distance);
        }
    }
    return max_distance;
}


void Matter::resize(const long int length)
{
 
    if(length>0) 
    {
        potential = new Potentials(parameters);
        
        nAtoms = length;
        positions.resize(length, 3);
        positions.setZero();
        
        velocities.resize(length ,3);
        velocities.setZero();

        forces.resize(length ,3);
        forces.setZero();
        
        masses.resize(length);
        masses.setZero();

        atomicNrs.resize(length);
        atomicNrs.setZero();
        
        isFixed.resize(length);
        isFixed.setZero();      
    }
}


long int Matter::numberOfAtoms() const {return(nAtoms);} // return the number of atoms


Vector3d Matter::getBoundary(int axis) const
{
    return(cellBoundaries.row(axis));
}

double Matter::getBoundary(int axis1, int axis2) const
{
    return cellBoundaries(axis1, axis2);
}


void Matter::setBoundary(int axis, Vector3d bound)
{
    cellBoundaries.row(axis)=bound;
    if(usePeriodicBoundaries)
    {
        applyPeriodicBoundary();
    }
    recomputePotential=true;
}

void Matter::setBoundary(int axis1, int axis2, double val)
{
    cellBoundaries(axis1,axis2)=val;
    if(usePeriodicBoundaries)
    {
        applyPeriodicBoundary();
    }
    recomputePotential=true;
}


void Matter::activatePeriodicBoundaries()
{
    usePeriodicBoundaries=true;
    applyPeriodicBoundary();
}


void Matter::deactivatePeriodicBoundaries()
{
    usePeriodicBoundaries=false;
}

double Matter::getPosition(long int indexAtom, int axis) const {
    return positions(indexAtom,axis);
}


void Matter::setPosition(long int indexAtom, int axis, double position)
{
    positions(indexAtom,axis)=position;
    if(usePeriodicBoundaries){ 
        applyPeriodicBoundary();
    }
    recomputePotential=true;
}


Matrix<double, Eigen::Dynamic, 3> Matter::getPositions() const
{//return coordinates of free atoms in array 'pos'
    return positions;
}


void Matter::setPositions(const Matrix<double, Eigen::Dynamic, 3> pos) {
// update Matter with the new positions of the free atoms given in array 'pos'
    positions = pos;
    if(usePeriodicBoundaries)
    {
        applyPeriodicBoundary();
    }
    recomputePotential=true;
}


Matrix<double, Eigen::Dynamic, 3> Matter::getForces() {
// return forces applied on all atoms in array 'force' 
    computePotential();
    Matrix<double, Eigen::Dynamic, 3> ret= forces;
    int i;
    for(i=0; i<nAtoms; i++)
    {
        if(isFixed[i])
        {
            ret.row(i).setZero();
        }
    }
    return ret;
}


double Matter::distance(long index1, long index2) const{
// return distance between the atoms with index1 and index2
    return pbc(positions.row(index1) - positions.row(index2)).norm();
}

double Matter::pdistance(long index1, long index2,int axis) const{
// return projected distance between the atoms with index1 and index2 on asix (0-x,1-y,2-z)
    Matrix<double, 1, 3> ret;
    ret.setZero();
    for(int i=0;i<3;i++){
        ret(0,i) = 0.0;
        if(i == axis){
            ret(0,i) = positions(index1,axis)-positions(index2,axis);
        }
    }
    return pbc(ret).norm();
}


double Matter::distance(const Matter& matter, long index) const {
// return the distance atom with index has moved between the current Matter object and the Matter object passed as argument
    return pbc(positions.row(index) - matter.getPositions().row(index)).norm();
}


double Matter::getMass(long int indexAtom) const
{
    return(masses[indexAtom]);
}


void Matter::setMass(long int indexAtom, double mass)
{
    masses[indexAtom]=mass;
}


long Matter::getAtomicNr(long int indexAtom) const
{
    return(atomicNrs[indexAtom]);
}


void Matter::setAtomicNr(long int indexAtom, long atomicNr)
{
    atomicNrs[indexAtom]=atomicNr;
    recomputePotential=true;
}

int Matter::getFixed(long int indexAtom) const {
    return(isFixed[indexAtom]);
}


void Matter::setFixed(long int indexAtom, int isFixed_passed)
{
    isFixed[indexAtom]=isFixed_passed;
}


double Matter::getPotentialEnergy()
{
    if(nAtoms>0) {
        computePotential();
        return potentialEnergy;
    } 
    else 
        return 0.0;
}


double Matter::getKineticEnergy() const
{
    double K=0;
    for(long int i=0; i<nAtoms; i++) {
            if(!isFixed[i]) K+=masses[i]*0.5*velocities.row(i).squaredNorm();
    };
    return K;
}


double Matter::getMechanicalEnergy()
{
    return getPotentialEnergy()+getKineticEnergy();
}


long int Matter::numberOfFreeAtoms() const {
    return nAtoms - isFixed.sum();
}

long Matter::getForceCalls() const{
    return(forceCalls);
}


void Matter::resetForceCalls(){
    forceCalls = 0;
    return;
}


// Print atomic coordinate to a .xyz file
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
        fprintf(file,"%s\t%f\t%f\t%f\n", atomicNumber2symbol(getAtomicNr(i)), getPosition(i, 0), getPosition(i, 1), getPosition(i, 2));
    };
    fclose(file);
}


// Print atomic coordinates to a .con file
bool Matter::matter2con(std::string filename) const {
    bool state;
    FILE *file;
    int pos=filename.find_last_of('.');
    if(filename.compare(pos+1, 3, "con")){
        filename+=".con";
    };
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
    int Ncomponent=0; // used to store the number of components (eg water: two components H and O)
    int first[MAXC]; // to store the position of the first atom of each component plus at the end the total number of atoms
    double mass[MAXC];
    long atomicNrs[MAXC];
    first[0]=0;
    if(numberOfAtoms()>0) {
        if(getFixed(0)) Nfix=1; // count the number of fixed atoms
        mass[0]=getMass(0);
        atomicNrs[0]=getAtomicNr(0);
    };
    j=0;
    for(i=1;i<numberOfAtoms();i++) {
        if(getFixed(i)) Nfix++; // count the number of fixed atoms
        if(getAtomicNr(i) != atomicNrs[j]) { // check if there is a second component
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
 
    fputs(headerCon1, file);
    fputs(headerCon2, file);
    double lengths[3];
    lengths[0] = cellBoundaries.row(0).norm();
    lengths[1] = cellBoundaries.row(1).norm();
    lengths[2] = cellBoundaries.row(2).norm();
    fprintf(file, "%f\t%f\t%f\n", lengths[0], lengths[1], lengths[2]);
    double angles[3];
    angles[0] = acos(cellBoundaries.row(0).dot(cellBoundaries.row(1))/lengths[0]/lengths[1])*180/M_PI;
    angles[1] = acos(cellBoundaries.row(0).dot(cellBoundaries.row(2))/lengths[0]/lengths[2])*180/M_PI;
    angles[2] = acos(cellBoundaries.row(1).dot(cellBoundaries.row(2))/lengths[1]/lengths[2])*180/M_PI;
    fprintf(file, "%f\t%f\t%f\n", angles[0], angles[1], angles[2]);
    fputs(headerCon5, file);
    fputs(headerCon6, file);
 
    fprintf(file, "%d\n", Ncomponent);
    for(j=0; j<Ncomponent; j++) {
        fprintf(file, "%d ", first[j+1]-first[j]);
    }
    fprintf(file, "\n");  
    for(j=0; j<Ncomponent; j++) {
//        mass[j]/=G_PER_MOL; // GH: I don't understand why we need to convert the mass units
        fprintf(file, "%f ", mass[j]);
    };
    fprintf(file, "\n");
    for(j=0; j<Ncomponent; j++) {
        fprintf(file, "%s\n", atomicNumber2symbol(atomicNrs[j]));
        fprintf(file, "Coordinates of Component %d\n", j+1);
        for(i=first[j]; i<first[j+1]; i++) {
            fprintf(file,"%.3f\t%.3f\t%.3f\t%d\t%ld\n", getPosition(i, 0), getPosition(i, 1), getPosition(i, 2), getFixed(i), i);
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
    file=fopen(filename.c_str(), "rb");
    if (!file) {
        cerr << "File " << filename << " was not found.\n";
        return(false);
    };
 
    state = con2matter(file);
    fclose(file);
    return(state);
}


bool Matter::con2matter(FILE *file) {
    char line[255]; // Temporary string of character to read from the file.
    fgets(headerCon1,sizeof(line),file);
    if (strchr(headerCon1,'\r')) {
        /* Files created on Windows or on Mac with Excell have carriage returns (\r) instead of or along
        with the new line charater (\n). C recognises only the \n as the end of line. */
        cerr << "A carriage return ('\\r') has been detected. To work correctly, new lines should be indicated by the new line character (\\n).";
        return false; // return false for error
    };

    long int i; int j;

    fgets(headerCon2,sizeof(line),file);


    double lengths[3];
    // The third line contains the length of the periodic cell
    fgets(line,sizeof(line),file);
    sscanf(line,"%lf %lf %lf", &lengths[0], &lengths[1], &lengths[2]);
 
    double angles[3];
    fgets(headerCon4,sizeof(line),file);
    // The fourth line contains the angles of the cell vectors
    sscanf(headerCon4,"%lf %lf %lf", &angles[0], &angles[1], &angles[2]); 

    if (angles[0] == 90.0 && angles[1] == 90.0 && angles[2] == 90.0) {
        cellBoundaries(0,0) = lengths[0];
        cellBoundaries(1,1) = lengths[1];
        cellBoundaries(2,2) = lengths[2];
    }else{
        angles[0] *= M_PI/180.0;
        angles[1] *= M_PI/180.0;
        angles[2] *= M_PI/180.0;

        cellBoundaries(0,0) = 1.0;
        cellBoundaries(1,0) = cos(angles[0]);
        cellBoundaries(1,1) = sin(angles[0]);
        cellBoundaries(2,0) = cos(angles[1]);
        cellBoundaries(2,1) = (cos(angles[2])-cellBoundaries(1,0)*cellBoundaries(2,0))/cellBoundaries(1,1);
        cellBoundaries(2,2) = sqrt(1.0-pow(cellBoundaries(2,0),2)-pow(cellBoundaries(2,1),2));

        cellBoundaries(0,0) *= lengths[0];
        cellBoundaries(1,0) *= lengths[1];
        cellBoundaries(1,1) *= lengths[1];
        cellBoundaries(2,0) *= lengths[2];
        cellBoundaries(2,1) *= lengths[2];
        cellBoundaries(2,2) *= lengths[2];
    }
 
    fgets(headerCon5,sizeof(line),file);
    fgets(headerCon6,sizeof(line),file);

    fgets(line,sizeof(line),file);
    int Ncomponent; // Number of components is the number of different types of atoms. For instance H2O (water) has two component (H and O).
    if(sscanf(line,"%d",&Ncomponent)==0) {
        std::cout << "The number of components cannot be read. One component is assumed instead\n";
        Ncomponent=1;
    }
    if((Ncomponent>MAXC)||(Ncomponent<1)) {
        cerr << "con2atoms does not support more than " << MAXC << " components (or less than 1).\n";
        return false;
    }
    /* to store the position of the 
        first atom of each element 'MAXC+1': the last element is used to store the total number of atom.*/ 
    long int first[MAXC+1];
    long int Natoms=0;
    first[0]=0;
    // Now we want to know the number of atom of each type. Ex with H2O, two hydrogens and one oxygen
    for(j=0; j<Ncomponent; j++) {
        fscanf(file, "%ld", &Natoms);
        first[j+1]=Natoms+first[j];
    }    
    fgets(line, sizeof(line), file); // Discard the rest of the line
    resize(first[Ncomponent]); // Set the total number of atoms, and allocates memory
    double mass[MAXC];
    for(j=0; j<Ncomponent; j++) { // Now we want to know the number of atom of each type. Ex with H2O, two hydrogens and one oxygen
        fscanf(file, "%lf", &mass[j]);
//        mass[j]*=G_PER_MOL; // conversion of g/mol to local units. (see su.h)
    };
    fgets(line,sizeof(line),file); // discard rest of the line
    int atomicNr;
    int fixed;
    double x,y,z;
    for (j=0; j<Ncomponent; j++) {
        char symbol[3];
        fgets(line,sizeof(line),file);
        sscanf(line, "%2s\n", symbol);
        atomicNr=symbol2atomicNumber(symbol);
        fgets(line,sizeof(line),file); // skip one line
        for (i=first[j]; i<first[j+1]; i++){
            setMass(i, mass[j]);
            setAtomicNr(i, atomicNr);
            fgets(line, sizeof(line), file);
            sscanf(line,"%lf %lf %lf %d\n", &x, &y, &z, &fixed);
            setPosition(i, 0, x);
            setPosition(i, 1, y);
            setPosition(i, 2, z);
            setFixed(i, static_cast<bool>(fixed));
        }
    }
    if(usePeriodicBoundaries)
    { 
        applyPeriodicBoundary(); // Transform the coordinate to use the minimum image convention.
    }
    //    potential_ = new Potentials(parameters_);
    return(true);
}


void Matter::computePotential()
{
    if(recomputePotential) {
        if(potential) {
            forces = potential->force(nAtoms, positions, atomicNrs, &potentialEnergy, cellBoundaries);
            forceCalls = forceCalls+1;
            recomputePotential=false;
        }
        else {
            cerr << "No potential associated with the atomic structure." << endl;
            exit(1);
        };
    };
}

void Matter::clearMemory()
{
    // the pointer to parameters should not be deleted, it is a reference to shared data
    //delete parameters_;
 
    if (potential!=0){
        delete potential;
        potential=0;
    }
}


void Matter::applyPeriodicBoundary() // Transform the coordinate to use the minimum image convention.
{ 
    Matrix<double, 3, 3> ibox = cellBoundaries.inverse();
    Matrix<double, Eigen::Dynamic, 3> ddiff = positions*ibox;
   
    int i,j;
    for(i=0; i<ddiff.rows(); i++)
    {
        for(j=0; j<3; j++)
        {
            ddiff(i,j) = fmod(ddiff(i,j), 1.0);
        }
    }

    positions= ddiff*cellBoundaries;
}



double Matter::maxForce(void)
{
    //Ensures that the forces are up to date
    computePotential();
    
    //I think this can be done in one line with the rowwise method
    double maxForce = 0.0;
    for(int i = 0; i < nAtoms; i++)
    {
        if(getFixed(i))
        {
            continue;
        }
        maxForce = max(forces.row(i).norm(), maxForce);
    }
    return maxForce;
}


bool Matter::isItConverged(double convergeCriterion)
{
    double diff=0;

    for(int i=0;i<nAtoms;i++)
    {
        if(getFixed(i))
        {
            continue;
        }
        diff = forces.row(i).norm();
        

        if(convergeCriterion < diff)
        { 
            break;
        }
    }
    return(diff < convergeCriterion);
}

Matrix<double, Eigen::Dynamic, 3> Matter::getFree() const
{
    Matrix<double, Eigen::Dynamic, 3> ret(nAtoms,3);
    int i,j;
    for(i=0;i<nAtoms;i++)
    {
        for(j=0;j<3;j++)
        {
            ret(i,j) = double(!bool(isFixed(i)));
        }
    }
    return ret;
}

Matrix<double, Eigen::Dynamic, 3> Matter::getVelocities() const
{
    return velocities.cwise() * getFree(); 
}

void Matter::setVelocities(const Matrix<double, Eigen::Dynamic, 3> v)
{
    velocities = v.cwise()*getFree(); 
}

void Matter::setForces(const Matrix<double, Eigen::Dynamic, 3> f)
{
    forces = f.cwise()*getFree();
}


Matrix<double, Eigen::Dynamic, 3> Matter::getAccelerations()
{
    Matrix<double, Eigen::Dynamic, 3> ret =  getForces().cwise() * getFree(); 
    ret.col(0).cwise()/=masses;
    ret.col(1).cwise()/=masses;
    ret.col(2).cwise()/=masses;
    return ret;
}

Matrix<double, Eigen::Dynamic, 1> Matter::getMasses() const
{
    return masses;
}
