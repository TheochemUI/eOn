#include "Hessian.h"

Hessian::Hessian(const Matter *react, const Matter *sad, const Matter *prod, Parameters* params)
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
		} while (parameters->hessianMaxSize < atoms.cols()*3);
    }
    size = atoms.cols()*3;
    
    //Build the hessian 
    int iAtom, jAtom;
    
    Matter matterTemp(parameters);
    matterTemp = *matter;
    
    double dr = 1e-6; // value used by graeme 1e-4;
    double tdr = 2*dr;
    
    Matrix<double, Eigen::Dynamic, 3> pos = matter->getPositions();
    Matrix<double, Eigen::Dynamic, 3> posDisplace(nAtoms, 3);
    Matrix<double, Eigen::Dynamic, 3> posTemp(nAtoms, 3);
    Matrix<double, Eigen::Dynamic, 3> force1(nAtoms, 3);
    Matrix<double, Eigen::Dynamic, 3> force2(nAtoms, 3);

    Matrix <double, Eigen::Dynamic, Eigen::Dynamic> hessian((int)sizeHessian, (int)sizeHessian);


    //----- Initialize end -----
    //std::cout<<"determineHessian\n";
    int i, j; 
    for(i = 0; i<size; i++)
    { 
        posDisplace.setZero(); 
        
        // Displacing one coordinate
        posDisplace(iCoord, i%3)  = dr;
        
        posTemp = pos - posDisplace;
        matterTemp.setPositions(posTemp);
        force1 = matterTemp.getForces();
        
        posTemp = pos + posDisplace; 
        matterTemp.setPositions(posTemp);
        force2 = matterTemp.getForces();
        
        for(j=0; j<size; j++)
        {     
            jCoord = coords[j];
            hessian(i,j) = -(force2(jCoord/3, jCoord%3)-force1(jCoord/3, jCoord%3))/tdr;
        }
    }
    
    forceCallsAtomsTemp = matterTemp.getForceCalls()-forceCallsAtomsTemp;
    
    totalForceCalls += forceCallsAtomsTemp;
        

}

VectorXi Hessian:movedAtoms(double const distance) const
{
    long nAtoms = saddle->numberOfAtoms();
    
    VectorXi moved(nAtoms);

    Matrix<double, Eigen::Dynamic, 3> diffProd = saddle->pbc(saddle.getPositions() - product->getPositions());
    Matrix<double, Eigen::Dynamic, 3> diffReact = saddle->pbc(saddle.getPositions() - reactant->getPositions());

    diffProd.cwise() *= saddle->getFree();
    diffReact.cwise() *= saddle->getFree();
    
    int nMoved = 0;
    for(int i=0; i<nAtoms; i++)
    {
        if(diffProd.row(i).norm() > distance || diffReact.row(i).norm() > distance)
        {
            moved[nMoved] = i;
            nMoved++;
        }
    }
    return (VectorXi)moved.block(0,0,nMoved,1);
}


    
    
