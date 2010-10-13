/*
 *===============================================
 *  EON Matter.cpp
 *===============================================
 */

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
        //Invalid symbol.
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
    
    constraints = matter.constraints;
    usePeriodicBoundaries = matter.usePeriodicBoundaries;
 
    potentialEnergy = matter.potentialEnergy;
    recomputePotential = matter.recomputePotential;
 
    strcpy(headerCon1,matter.headerCon1);
    strcpy(headerCon2,matter.headerCon2);
    strcpy(headerCon4,matter.headerCon4);
    strcpy(headerCon5,matter.headerCon5);
    strcpy(headerCon6,matter.headerCon6);

	//liang add here nsteps for test:
	nsteps = matter.nsteps;

    return *this;
}


bool Matter::operator==(const Matter& matter) {//To compare two matter objects, if all differences in positions are bellow getMaxDifferencePos, they are considered equal
    return (parameters->maxDifferencePos) > perAtomNorm(matter);
}


double Matter::distanceTo(const Matter& matter) 
{
    /* RT: Returns the distance to the given matter object. */

    return sqrt(pbc(positions - matter.positions).cwise().square().sum());
}

Matrix<double, Eigen::Dynamic, 3> Matter::pbc(Matrix<double, Eigen::Dynamic, 3> diff) const
{
    Matrix<double, 3, 3> ibox = cellBoundaries.inverse();
    Matrix<double, Eigen::Dynamic, 3> ddiff = diff*ibox;
    
    int i,j;
    for(i=0; i<nAtoms; i++)
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
        Matrix<double, Eigen::Dynamic, 3> diff = positions - matter.positions;
        for(i = 0; i < nAtoms; i++)
        {
            max_distance = max(diff.row(i).norm(), max_distance);
        }
    }
    return max_distance;
}


void Matter::resize(const long int length)
{
 
    if(nAtoms>0) 
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


long int Matter::numberOfAtoms() const {return(nAtoms);}// return the number of atoms


Vector3d Matter::getBoundary(int axis) const
{
    return(cellBoundaries.row(axis));
}


void Matter::setBoundary(int axis, Vector3d bound)
{
    cellBoundaries.row(axis)=bound;
    if(usePeriodicBoundaries)
    {
        applyPeriodicBoundary(axis);
    }
    recomputePotential=true;
}


void Matter::activatePeriodicBoundaries()
{
    usePeriodicBoundaries=true;
    applyPeriodicBoundary(0);
    applyPeriodicBoundary(1);
    applyPeriodicBoundary(2);
}


void Matter::deactivatePeriodicBoundaries()
{
    usePeriodicBoundaries=false;
}


void Matter::setConstraints(Constraints constraints_passed)
{
    constraints=constraints_passed;
}


double Matter::getPosition(long int indexAtom, int axis) const {
    return positions(indexAtom,axis);
}


void Matter::setPosition(long int indexAtom, int axis, double position)
{
    positions(indexAtom,axis)=position;
    if(usePeriodicBoundaries){ 
        applyPeriodicBoundary(indexAtom, axis);
    }
    recomputePotential=true;
}


Matrix<double, Eigen::Dynamic, 3> Matter::getPositions() const
{//return coordinates of free atoms in array 'pos'
    return positions;
}


void Matter::setPositions(const Matrix<double, Eigen::Dynamic, 3> pos) {//Update Matter with the new positions of the free atoms given in array 'pos'
    positions = pos;
    if(usePeriodicBoundaries)
    {
        applyPeriodicBoundary();
    }
    recomputePotential=true;
}


Matrix<double, Eigen::Dynamic, 3> Matter::getForces() {// return forces applied on all atoms in array 'force' 
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


double Matter::distance(long index1, long index2) const{// return distance between the atoms with index1 and index2
    return pbc(positions.row(index1) - positions.row(index2)).norm();
}


double Matter::distance(const Matter& matter, long index) const {// return the distance atom with index has moved between the current Matter object and the Matter object passed as argument
    return pbc(positions.row(index) - matter.positions.row(index)).norm();
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

//liang add here
void Matter::setNsteps(long int Nsteps){
	 nsteps=Nsteps;
}

long Matter::getNsteps() const{
	return(nsteps);
}


int Matter::getFixed(long int indexAtom) const {
    return(isFixed[indexAtom]);
}


void Matter::setFixed(long int indexAtom, int isFixed_passed)
{
    isFixed[indexAtom]=isFixed_passed;
}


double Matter::getPotentialEnergy() const
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


double Matter::getMechanicalEnergy() const
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
        if(getFixed(i)) Nfix++; //count the number of fixed atoms
        if(getAtomicNr(i) != atomicNrs[j]) { // Check if there is a second component
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
    fprintf(file, "%f\t%f\t%f\n", cellBoundaries(0,0), cellBoundaries(1,1), cellBoundaries(2,2)); //XXX: Orthoganal boxes only!!
    fputs(headerCon4, file);
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
 
    fgets(line,sizeof(line),file);
 
    double x, y, z;
    sscanf(line,"%lf %lf %lf", &x, &y, &z); // The third line contains the length of the periodic cell
    cellBoundaries(0,0)= x;
    cellBoundaries(1,1)= y;
    cellBoundaries(2,2)= z;
 
    fgets(headerCon4,sizeof(line),file);
    sscanf(headerCon4,"%lf %lf %lf", &x, &y, &z); // The fourth line contains the angles of the cell vectors
    if ( (x != 90.0) or (y != 90.0) or (z != 90.0) ) {
        /* The code only supports cubic simulation cells*/
        cerr << "This code only supports rectangular cells.";
        return false;// return false for error
    };
 
    fgets(headerCon5,sizeof(line),file);
    fgets(headerCon6,sizeof(line),file);
 
    fgets(line,sizeof(line),file);
    int Ncomponent; // Number of components is the number of different types of atoms. For instance H2O (water) has two component (H and O).
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
    fgets(line, sizeof(line), file); // Discard the rest of the line
    resize(first[Ncomponent]); // Set the total number of atoms, and allocates memory
    double mass[MAXC];
    for(j=0; j<Ncomponent; j++) { // Now we want to know the number of atom of each type. Ex with H2O, two hydrogens and one oxygen
        fscanf(file, "%lf", &mass[j]);
//        mass[j]*=G_PER_MOL; // conversion of g/mol to local units. (see su.h)
    };
    fgets(line,sizeof(line),file); //Discard rest of the line
    int atomicNr;
    int fixed;
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
        };
    };
    if(usePeriodicBoundaries)
    { 
        applyPeriodicBoundary(); // Transform the coordinate to use the minimum image convention.
    }
    //    potential_ = new Potentials(parameters_);
    return(true);
}


void Matter::computePotential() const
{
    if(recomputePotential) {
        if(potential) {
            potential->force(nAtoms, positions, atomicNrs, forces, &potentialEnergy, cellBoundaries);
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
    applyPeriodicBoundary(0);
    applyPeriodicBoundary(1);
    applyPeriodicBoundary(2);
}


void Matter::applyPeriodicBoundary(int axis) 
{
    for(long int i=0; i<nAtoms; i++)
    {
        applyPeriodicBoundary(i, axis);
    }
}


void Matter::applyPeriodicBoundary(long atom, int axis)
{
    //TODO: implement this correctly (do we even need it??)
    return;
    //while(positions_[atom*3+axis]>cellBoundaries_[axis]) {
    //    positions_[atom*3+axis]-=cellBoundaries_[axis];
    //};
    //while(positions_[atom*3+axis]<= 0.0) {
    //    positions_[atom*3+axis]+=cellBoundaries_[axis];
    //};
}

double Matter::maxForce(void)
{
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

    for(int i=0;i<nAtoms*3;i++)
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
};

Matrix<double, Eigen::Dynamic, 3> Matter::getFree() const
{
    Matrix<double, Eigen::Dynamic, 3> ret;
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
