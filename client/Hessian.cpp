#include <math.h>

#include "Hessian.h"

Hessian::Hessian(Matter *react, Matter *sad, Matter *prod, Parameters* params)
{
    reactant = react;
    saddle = sad;
    product = prod;
    parameters = params;
}

Hessian::~Hessian()
{
}

double Hessian::getModeProduct(int which)
{
    Matter *curr;
    assert(saddle->numberOfAtoms() == reactant->numberOfAtoms());
    assert(saddle->numberOfAtoms() == product->numberOfAtoms());
    switch(which)
    {
        case REACTANT:
            curr = reactant;
            break;
        case SADDLE:
            curr = saddle;
            break;
        case PRODUCT:
            curr = product;
            break;
        default:
            cerr<<"Hessian can't deterimine which structure to use"<<endl;
            return -1;
            break;
    }

    int nAtoms = curr->numberOfAtoms();

    //Determine which atoms moved in the process
    int size=0;
    VectorXi atoms;
    if(parameters->hessianMaxSize == 0)
    {
        atoms = movedAtoms(parameters->hessianMinDisplacement);
    }
    else
    {
		int loop = 0;
		do {
			double minDisp = parameters->hessianMinDisplacement + loop * 0.1;		
			atoms = movedAtoms(minDisp);
			loop = loop + 1;
		} while (parameters->hessianMaxSize < atoms.rows()*3);
    }
    size = atoms.rows()*3;
    cout<<"Hessian size: "<<size<<endl;
    assert(size > 0);
    
    //Build the hessian 
    Matter matterTemp(parameters);
    matterTemp = *curr;
    
    double dr = 1e-6; // value used by graeme 1e-4;
    double tdr = 2*dr;
    
    Matrix<double, Eigen::Dynamic, 3> pos = curr->getPositions();
    Matrix<double, Eigen::Dynamic, 3> posDisplace(nAtoms, 3);
    Matrix<double, Eigen::Dynamic, 3> posTemp(nAtoms, 3);
    Matrix<double, Eigen::Dynamic, 3> force1(nAtoms, 3);
    Matrix<double, Eigen::Dynamic, 3> force2(nAtoms, 3);

    Matrix <double, Eigen::Dynamic, Eigen::Dynamic> hessian(size, size);

    int i, j; 
    for(i = 0; i<size; i++)
    { 
        posDisplace.setZero(); 
        
        // Displacing one coordinate
        posDisplace(atoms(i/3), i%3)  = dr;
        
        posTemp = pos - posDisplace;
        matterTemp.setPositions(posTemp);
        force1 = matterTemp.getForces();
        
        posTemp = pos + posDisplace; 
        matterTemp.setPositions(posTemp);
        force2 = matterTemp.getForces();
        
        for(j=0; j<size; j++)
        {     
            hessian(i,j) = -(force2(atoms(j/3), j%3)-force1(atoms(j/3), j%3))/tdr;
            double effMass=sqrt(saddle->getMass(j/3)*saddle->getMass(i/3));
            hessian(i,j)/=effMass;
        }
    }

    //Force hessian to be symmetric
    hessian = (hessian + hessian.transpose())/2;

    
    Eigen::SelfAdjointEigenSolver<MatrixXd> es(hessian);
    VectorXd freqs = es.eigenvalues();

    int nNeg = 0;
    double prod = 1;
    for(i=0; i<size; i++)
    {
        if(freqs(i) < 0)
        {
            nNeg++;
        }
        else
        {
            prod *= freqs(i);
        }
    }
    if(which == SADDLE)
    {
        if(nNeg!=1)
        {
            return -1;
        }
    }
    else
    {
        if(nNeg!=0)
        {
            return -1;
        }
    }

    return prod;    
}

VectorXi Hessian::movedAtoms(double const distance)
{
    long nAtoms = saddle->numberOfAtoms();
    
    VectorXi moved(nAtoms);
    moved.setConstant(-1);

    Matrix<double, Eigen::Dynamic, 3> diffProd = saddle->pbc(saddle->getPositions() - product->getPositions());
    Matrix<double, Eigen::Dynamic, 3> diffReact = saddle->pbc(saddle->getPositions() - reactant->getPositions());

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
                
                if(diffRSaddle<parameters->hessianWithinRadiusDisplaced 
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


    
    
