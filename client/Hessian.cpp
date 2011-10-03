//-----------------------------------------------------------------------------------
// eOn is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// A copy of the GNU General Public License is available at
// http://www.gnu.org/licenses/
//-----------------------------------------------------------------------------------

#include <math.h>
#include "Hessian.h"

const char Hessian::REACTANT[] =    "reactant";
const char Hessian::SADDLE[] =      "saddle";
const char Hessian::PRODUCT[] =     "product";


Hessian::Hessian(Matter *react, Matter *sad, Matter *prod, Parameters* params)
{
    reactant = react;
    saddle = sad;
    product = prod;
    parameters = params;

    for(int i=0; i<3; i++)
    {
        modes[i].resize(0);
        hessians[i].resize(0,0);
    }
}

Hessian::~Hessian()
{
}

int Hessian::whichNum(string which)
{
    if(which == REACTANT)   {return 0;}
    if(which == SADDLE)     {return 1;}
    if(which == PRODUCT)    {return 2;}
    return -1;
}    

Matrix<double, Eigen::Dynamic, Eigen::Dynamic> Hessian::getHessian(string which)
{
    if(hessians[whichNum(which)].rows() == 0)
    {
        if (!calculate(which))
        {
            hessians[whichNum(which)].resize(0,0);
        }
    }
    return hessians[whichNum(which)];
}

VectorXd Hessian::getModes(string which)
{
    if(modes[whichNum(which)].rows() == 0)
    {
        if (!calculate(which))
        {
            modes[whichNum(which)].resize(0);
        }
    }
    return modes[whichNum(which)];
}

bool Hessian::calculate(string which)
{
    Matter *curr;
    assert(saddle->numberOfAtoms() == reactant->numberOfAtoms());
    assert(saddle->numberOfAtoms() == product->numberOfAtoms());

    if(which == REACTANT)       {curr = reactant;}
    else if(which == SADDLE)    {curr = saddle;}
    else if(which == PRODUCT)   {curr = product;}
    else
    {
        cerr<<"Hessian can't deterimine which structure to use"<<endl;
        return false;
    }

    int nAtoms = curr->numberOfAtoms();

    //Determine which atoms moved in the process
    int size = 0;
    VectorXi atoms;
    atoms = movedAtoms(parameters->hessianMinDisplacement);
    size = atoms.rows()*3;
    cout<<"Hessian size: "<<size<<endl;
    assert(size > 0);

    //Build the hessian 
    Matter matterTemp(parameters);
    matterTemp = *curr;

//    double dr = 1e-6; // value used by graeme 1e-4;
    double dr = parameters->hessianFiniteDist;

    AtomMatrix pos = curr->getPositions();
    AtomMatrix posDisplace(nAtoms, 3);
    AtomMatrix posTemp(nAtoms, 3);
    AtomMatrix force1(nAtoms, 3);
    AtomMatrix force2(nAtoms, 3);

    Matrix <double, Eigen::Dynamic, Eigen::Dynamic> hessian(size, size);

    int i, j; 
    force1 = matterTemp.getForces();
    for(i = 0; i<size; i++)
    {
        posDisplace.setZero(); 

        // Displacing one coordinate
        posDisplace(atoms(i/3), i%3) = dr;

        posTemp = pos + posDisplace; 
        matterTemp.setPositions(posTemp);
        force2 = matterTemp.getForces();

        //To use central difference estimate of the hessian uncomment following (and divide by 2*dr) in the following.
        //This does use an additional 'size' forcecalls and will generally not lead to very different results.
        //In most cases, the additional accuracy is not worth the computation time.

        /*
        posTemp = pos - posDisplace; 
        matterTemp.setPositions(posTemp);
        force1 = matterTemp.getForces();
        */

        // GH: debug
        //cout <<"force"<<i<<endl<<force2<<endl;

        for(j=0; j<size; j++)
        {
            hessian(i,j) = -(force2(atoms(j/3), j%3)-force1(atoms(j/3), j%3))/dr;
            //get the effective mass of the moving atoms 
            //double effMass = sqrt(saddle->getMass(j/3)*saddle->getMass(i/3));
            double effMass = sqrt(saddle->getMass(atoms(j/3))*saddle->getMass(atoms(i/3)));
            hessian(i,j) /= effMass;
        }
    }

    //Force hessian to be symmetric
    
    //hessian = (hessian + hessian.transpose())/2;  
    //cannot be used, messes up the lower trianguler 
    //transpose does not seem to be a hardcopy, rather just an index manipulation 

    for(i=0; i<size; i++) {
        for(j=0; j<i; j++) {
            hessian(i,j) = ( hessian(i,j) + hessian(j,i) ) / 2;
            hessian(j,i) = hessian(i,j);    
        }
    }

    // GH: debug
    if(!parameters->quiet)
    {
        cout <<"writing hessian"<<endl;
        ofstream hessfile;
        hessfile.open("hessian.dat");
        hessfile <<hessian;
    }

    Eigen::SelfAdjointEigenSolver<MatrixXd> es(hessian);
    VectorXd freqs = es.eigenvalues();

    // GH debug
//    if(!parameters->quiet)
//    {
//        cout << "freqs\n" << freqs << endl;
//    }
    
    //If we are checking for rotation, then the system has no frozen atoms and
    //can rotate and translate. This gives effectively zero eigenvalues. We
    //need to remove them from the prefactor calculation. 
    //
    //the second condition requires that every atom moves. Otherwise, we don't 
    //get the 6 rotational and translational modes.
    //XXX: what happens if the entire particle rotates about one atom or a line of atoms?
    if(parameters->checkRotation && size==3*saddle->numberOfAtoms())
    {
        VectorXd newfreqs(size);
        int nremoved = 0;
        for(i=0; i<size; i++)
        {
            if(abs(freqs(i)) > 1e-6) //XXX: Hardcoded
            {
                newfreqs(i-nremoved) = freqs(i);
            }
            else
            {
                nremoved++;
            }
        }
        if(nremoved != 6)
        {
            cout<<"Error: Found "<<nremoved<<" trivial eigenmodes instead of 6."<<endl;
            return false;
        }
        freqs = newfreqs;
    }

    int nNeg = 0;
    for(i=0; i<size; i++)
    {
        if(freqs(i) < 0)
        {
            nNeg++;
        }
    }
    if(which == SADDLE)
    {
        if(nNeg!=1)
        {
            cout<<"Error: "<<nNeg<<" negative modes at the saddle"<<endl;
            return false;
        }
    }
    else
    {
        if(nNeg!=0)
        {
            cout<<"Error: "<<nNeg<<" negative modes at the reactant/product"<<endl;
            return false;
        }
    }
    modes[whichNum(which)] = freqs;
    hessians[whichNum(which)] = hessian;
    return true;
}

VectorXi Hessian::movedAtoms(double const distance)
{
    long nAtoms = saddle->numberOfAtoms();

    VectorXi moved(nAtoms);
    moved.setConstant(-1);

    AtomMatrix diffProd = saddle->pbc(saddle->getPositions() - product->getPositions());
    AtomMatrix diffReact = saddle->pbc(saddle->getPositions() - reactant->getPositions());

    diffProd.cwise() *= saddle->getFree();
    diffReact.cwise() *= saddle->getFree();

    int nMoved = 0;
    for(int i=0; i<nAtoms; i++)
    {
        if( (diffProd.row(i).norm() > distance) || (diffReact.row(i).norm() > distance))
        {
            if(!(moved.cwise() == i).any())
            {
                moved[nMoved] = i;
                nMoved++;
            }
            for(int j=0; j<nAtoms; j++)
            {
                double diffRSaddle = saddle->distance(i,j);
                
                if(diffRSaddle<parameters->hessianWithinRadius
                   && (!saddle->getFixed(j)))
                {
                    if(!(moved.cwise() == j).any())
                    {
                        moved[nMoved] = j;
                        nMoved++;
                    }
                }
            }
        }
    }
    return (VectorXi)moved.block(0,0,nMoved,1);
}

