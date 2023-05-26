// EMT.cpp  --  Implements the EMT potential for multiple elements.
// When done, this will make EMTPotential.{h,cpp} obsolete.

#include "EMT.h"
#include "Atoms.h"
#include "EMTDefaultParameterProvider.h"
#include "Exception.h"
#include "GhostAtoms.h"
#include "GhostPotential.h"
#include "NeighborList.h"
#include "SuperCell.h"
#include "Vec.h"
#include "vectools.h"
// #include "Timing.h"
#include "mass.h"
#include <assert.h>
#include <math.h>
#include <set>
#include <stdlib.h>
#include <string.h>
#include <vector>
using std::cerr;
using std::endl;
using std::flush;
using std::set;
using std::vector;

int verbose = 0;

#if verbose
#define DEBUGPRINT                                                             \
  cerr << "Now in " << __FILE__ << " line " << __LINE__ << endl << flush;
#else
#define DEBUGPRINT
#endif

#if verbose
#define VERB(x)                                                                \
  if (verbose == 1)                                                            \
  cerr << x
#else
#define VERB(x)
#endif

using namespace std;

// Maximal number of elements in a single simulation
#define NMAXELEMENTS 10
// The batch--buffer size.  Avoid powers of 2
#define BUFLEN 10000

// Standard mapping of the six independent parts of the stress tensor to
// vector notation
const static int stresscomp[3][3] = {{0, 5, 4}, {5, 1, 3}, {4, 3, 2}};

EMT::EMT(EMTParameterProvider *prov) {
  DEBUGPRINT;
  if (prov == 0) {
    ownProvider = 1;
    provider = new EMTDefaultParameterProvider();
  } else {
    ownProvider = 0;
    provider = prov;
  }
  nblist = 0;
  atoms = 0;
  nAtoms = 0;
  nSize = 0;
  subtractE0 = true;
  counters.ids = counters.nblist = counters.sigma1 = counters.sigma2 =
      counters.energies = counters.forces = counters.stresses =
          counters.fullstresses = counters.beforeforces = 0;
  coef = 0;
  DEBUGPRINT;
}

void EMT::SetAtoms(Atoms *listofatoms) {
  DEBUGPRINT;
  if (atoms != 0)
    throw Exception("EMT::SetAtoms: Atoms already assigned to this potential.");
  atoms = listofatoms;
  //  try{
  //
  //	 ghostatoms = dynamic_cast<GhostAtoms *>(atoms);
  //  }catch(...)
  //  {
  ghostatoms = NULL;
  //  }
  nAtoms = atoms->GetNumberOfRealAtoms();
  nSize = 0; // Will be set when the NB list is created.  The value 0
  // should create enough trouble to detect if it is used before then.

  InitParameters();
  if (nelements == 1) {
    // Make the single element the continuum element.
    SetContinuumElement(parameters[0]->Z);
  } else {
    // Otherwise no meaningful default is set, the user must choose.
    SetContinuumElement(0);
  }
  DEBUGPRINT;
}

/// If an EMTParameterProvider was given when constructing the EMT
/// object, it will NOT be deleted, but any default parameter provider
/// and the automatically generated neigbor list will be deleted.
EMT::~EMT() {
  if (ownProvider)
    delete provider;
  if (nblist)
    delete nblist;
  if (coef)
    delete[] coef;
}

void EMT::Allocate() {
  int i;
  // USETIMER("EMT::Allocate");
  DEBUGPRINT;
  VERB(" Allocate[" << nAtoms << "," << nSize << "]" << flush);
  // WARNING: Resizing the vector may allocate way too much memory.  It
  // appears that calling reserve solves this problem.  For efficiency,
  // reserve is called with 5% extra space.  This is only necessary if the
  // atoms have ghosts, otherwise no reallocation will happen.

  // First, check if reallocation is necessary.
  if (nSize != int(Ec.size()) || nAtoms != int(Eas.size())) {
    DEBUGPRINT;
    /* Resize/intialize the internal variables. */
    sigma1.resize(nelements);
    sigma2.resize(nelements);
    // Do the reserve trick if the atoms have ghosts.
    if (ghostatoms) {
      if (int(Ec.capacity()) < nSize) {
        nSizeRes = nSize + nSize / 20;
        for (i = 0; i < nelements; i++) {
          sigma1[i].reserve(nSizeRes);
          sigma2[i].reserve(nSizeRes);
        }
        Ec.reserve(nSizeRes);
        dEds.reserve(nSizeRes);
        radius.reserve(nSizeRes);
        id.reserve(nSizeRes);
      }
      if (int(Eas.capacity()) < nAtoms) {
        nAtomsRes = nAtoms + nAtoms / 20;
        Eas.reserve(nAtomsRes);
        force.reserve(nAtomsRes);
        potentialenergy.reserve(nAtomsRes);
      }
    }

    for (i = 0; i < nelements; i++) {
      sigma1[i].resize(nSize);
      sigma2[i].resize(nSize);
    }
    DEBUGPRINT;
    Ec.resize(nSize);
    Eas.resize(nAtoms);
    potentialenergy.resize(nAtoms);
    radius.resize(nSize);
    dEds.resize(nSize);
    id.resize(nSize);
    force.resize(nAtoms);
    if (stress.size())
      AllocateStress();
    // The IDs are set to zero, then they do not have to be calculated if
    // there is only one element.
    if (nelements == 1)
      for (i = 0; i < nSize; i++)
        id[i] = 0;
  }
  DEBUGPRINT;
}

// Stresses are reallocated if they exist when Allocate is called, or if
// they do not exist, and GetStress is called.
void EMT::AllocateStress() {
  DEBUGPRINT;
  if (ghostatoms && (int(stress.capacity()) < 6 * nAtoms))
    stress.reserve(6 * nAtomsRes);
  stress.resize(6 * nAtoms);
  DEBUGPRINT;
}

void EMT::InitParameters() {
  DEBUGPRINT;
  // Extract the elements from the list of atoms, and set up an array
  // of EMT parameters.
  set<int> elements;

  // Extract the elements from the list of atoms.
  // elements.clear();
  atoms->GetListOfElements(elements);
  nelements = elements.size();

  // Get the EMT parameters
  parameters.clear();
  for (set<int>::iterator i = elements.begin(); i != elements.end(); ++i)
    parameters.push_back(provider->GetParameters(*i));

  // assert(nelements == provider->GetNumberOfElements());

  // Calculate the quantities that can only be calculated, once the
  // elements in the simulations are known.
  provider->CalcGammaEtc();
  rFermi = provider->GetCutoffDistance();
  // rNbCut = 1.04500185048 * rFermi;
  rNbCut = provider->GetListCutoffDistance();
  cutoffslope = provider->GetCutoffSlope();
  chi = provider->GetChi();
  if (verbose)
    cerr << "EMT::InitParameters:  rFermi = " << rFermi
         << "  rNbCut = " << rNbCut << "  cutoffslope = " << cutoffslope
         << endl;
  DEBUGPRINT;
}

const Vec *EMT::GetCartesianForces() {
  // USETIMER("EMT::GetCartesianForces")
  DEBUGPRINT;
  if (counters.forces != atoms->GetChangeCounter()) {
    DEBUGPRINT;
    VERB(" Force[");
    CheckNeighborLists();
    CalculateIDs();
    CalculateSigmas(0);              /* Skip sigma2 */
    CalculateEnergiesAfterSigmas(0); /* Skip actual energy calculation */
    double(*str)[6] = 0;
    if (stress.size())
      str = (double(*)[6]) & (stress[0]);
    if (nelements > 1)
      CalculateForcesAfterEnergies(&(force[0]), str);
    else
      CalculateForcesAfterEnergiesSingle(&(force[0]), str);
    VERB("]" << flush);
  } else {
    VERB(" Force[-]");
  }
  DEBUGPRINT;
  return &(force[0]);
}

const symTensor *EMT::GetStresses(const Vec *momenta) {
  int i;
  // USETIMER("EMT::GetStresses")
  DEBUGPRINT;
  VERB(" Stresses[");

  // Workaround for lack of vector<double[6]>
  double(*stress)[6];

  if (counters.stresses != atoms->GetChangeCounter()) {
    DEBUGPRINT;
    CheckNeighborLists();
    if (this->stress.size() == 0)
      AllocateStress();
    CalculateIDs();
    CalculateSigmas(0);              /* Skip sigma2 */
    CalculateEnergiesAfterSigmas(0); /* Skip actual energy calculation */

    stress = (double(*)[6]) & (this->stress[0]);
    if (nelements > 1)
      CalculateForcesAfterEnergies(&(force[0]), stress);
    else
      CalculateForcesAfterEnergiesSingle(&(force[0]), stress);
  }
  stress = (double(*)[6]) & (this->stress[0]);
  if (counters.fullstresses != atoms->GetChangeCounter()) {
    VERB("S");
    DEBUGPRINT;
    counters.fullstresses = atoms->GetChangeCounter();
    // Now calculate the dynamic part of the stress, and normalize with the
    // volume.
    DEBUGPRINT;
    for (i = 0; i < 6; i++)
      unnormalizedstress[i] = 0.0;
    DEBUGPRINT;
    totalvolume = 0.0;
    int *id = &(this->id)[0]; // This is stepped up in the loop.
    DEBUGPRINT;
    assert(this->stress.size() == 6 * nAtoms);
    for (i = 0; i < nAtoms; i++, id++) {
      const double fourpithird = 4.1887902048; // 4 * PI / 3
      double s0 = parameters[*id]->seq;
      double mass = parameters[*id]->mass; /* Get it from Python ? */
      double vol = fourpithird * s0 * s0 * s0;
      double invvol = 1.0 / vol;
      for (int alpha = 0; alpha < 3; alpha++)
        for (int beta = alpha; beta < 3; beta++) {
          int j = stresscomp[alpha][beta];
          if (momenta)
            stress[i][j] -= momenta[i][alpha] * momenta[i][beta] / mass;
          unnormalizedstress[j] += stress[i][j];
          stress[i][j] *= invvol;
        }
      totalvolume += vol;
    }
  } else {
    VERB("-");
  }
  VERB("]" << flush);
  DEBUGPRINT;
  return stress;
}

void EMT::GetStress(double totalstress[6], const Vec *momenta) {
  int i;
  // USETIMER("EMT::GetStress");
  const double(*s)[6];
  double vol;
  DEBUGPRINT;
  VERB(" Stress[");

  // Calculate the atomic stresses. Their sum is stored in unnormalizedstress.
  s = GetStresses(momenta);

  // Always use the volume of the supercell to normalize the total stress.
  // The alternative, using the sum of atomic volumes when the supercell
  // volume is ill-defined due to open boundary conditions, will fail for
  // parallel simulations.
  vol = atoms->GetSuperCell()->GetVolume();

  // Normalize the total stress with the volume.
  // Why add to old value ?  Because QCPotential need this.
  for (i = 0; i < 6; i++)
    totalstress[i] = (totalstress[i] + unnormalizedstress[i]) / vol;
  DEBUGPRINT;
  VERB("]" << flush);
}

const double *EMT::GetPotentialEnergies() {
  // USETIMER("EMT::GetPotentialEnergies");
  DEBUGPRINT;
  VERB(" Energies[");
  if (counters.energies != atoms->GetChangeCounter()) {
    DEBUGPRINT;
    CheckNeighborLists();
    CalculateIDs();
    CalculateSigmas();
  } else {
    VERB("-");
  }
  assert(potentialenergy.size() == nAtoms);
  CalculateEnergiesAfterSigmas(&(potentialenergy[0]));
  // Does its own counter checking.
  DEBUGPRINT;
  VERB("]" << flush);
  return &(potentialenergy[0]);
}

double EMT::GetPotentialEnergy() {
  int i;
  // USETIMER("EMT::GetPotentialEnergy");
  VERB(" Energy[");
  DEBUGPRINT;
  double etot = 0.0;
  const double *e;

  e = GetPotentialEnergies();
  for (i = 0; i < nAtoms; i++)
    etot += e[i];
  DEBUGPRINT;
  VERB("]" << flush);
  return etot;
}

void EMT::CheckNeighborLists() {
  // USETIMER("EMT::CheckNeighborLists");
  DEBUGPRINT;
  if (counters.nblist == atoms->GetChangeCounter())
    return;
  VERB("n");
  counters.nblist = atoms->GetChangeCounter();
  DEBUGPRINT;
  if (nblist) {
    DEBUGPRINT;
    // Check if the neighborlist is up to date.
    bool updated = nblist->CheckAndUpdateNeighborList();
    if (updated) {
      nAtoms = atoms->GetNumberOfRealAtoms();
      nSize = nAtoms;
      if (ghostatoms)
        nSize += ghostatoms->GetNumberOfGhosts();
      Allocate();
    }
  } else {
    DEBUGPRINT;
    // First call, create the neighbor list.
    nblist = new NeighborList(atoms, rNbCut);
    nAtoms = atoms->GetNumberOfRealAtoms();
    nSize = nAtoms;
    if (ghostatoms)
      nSize += ghostatoms->GetNumberOfGhosts();
    Allocate();
  }
  DEBUGPRINT;
}

/// IDs are a recoding of the atomic numbers used as indices into
/// various arrays.  Atomic numbers are not practical for this as they
/// are not contiguous.
void EMT::CalculateIDs() {
  int i;
  // USETIMER("EMT::CalculateIDs");
  DEBUGPRINT;
  if (counters.ids == atoms->GetChangeCounter())
    return;
  VERB("i");
  counters.ids = atoms->GetChangeCounter();

  DEBUGPRINT;
  // If there is only one element, all IDs are 0 and do not have to be
  // calculated.  Otherwise the ID is the offset into the list of elements.
  if (nelements > 1) {
    const int *z = atoms->GetAtomicNumbers();
    int found = 0; // REMOVE WHEN TESTED !   XXXXX
    for (i = 0; i < nelements; i++) {
      int zcand = parameters[i]->Z;
      for (int j = 0; j < nSize; j++)
        if (z[j] == zcand) {
          id[j] = i;
          found++;
        }
    }
    assert(found == nSize);
    // If we have ghost atoms, their IDs should be updated.
    //// // NO: This is now done when the ghost positions are updated.
    ////if (ghostatoms)
    ////  ghostatoms->GetGhostPotential()->CommunicateData(&id[0]);
  }
  DEBUGPRINT;
}

void EMT::CalculateSigmas(int calculatesigma2) {
  int i;
  // USETIMER("EMT::CalculateSigmas")
  DEBUGPRINT;
  int ctr = atoms->GetChangeCounter();
  if ((counters.sigma1 == ctr) &&
      (!calculatesigma2 || (counters.sigma2 == ctr)))
    return;
  counters.sigma1 = ctr;
  if (calculatesigma2) {
    counters.sigma2 = ctr;
    VERB("2");
  } else {
    VERB("1");
  }

  DEBUGPRINT;
  // const Vec *pos = atoms->GetCartesianPositions();
  int maxnblen = nblist->MaxNeighborListLength();
  // Buffer data:
  TinyMatrix<int> nbatch(nelements, nelements);
  TinyMatrix<int[BUFLEN]> self(nelements, nelements);
  TinyMatrix<int[BUFLEN]> other(nelements, nelements);
  TinyMatrix<Vec[BUFLEN]> rnb(nelements, nelements);
  TinyMatrix<double[BUFLEN]> sqdist(nelements, nelements);
  int other_buf[BUFLEN];
  Vec rnb_buf[BUFLEN];
  double sqdist_buf[BUFLEN];

  /* Set sigmas to zero */
  for (i = 0; i < nelements; i++) {
    memset(&sigma1[i][0], 0, nSize * sizeof(double));
    memset(&sigma2[i][0], 0, nSize * sizeof(double));
  }

  /* No atoms in batch pools */
  for (i = 0; i < nelements; i++)
    for (int j = 0; j < nelements; j++)
      nbatch[i][j] = 0;

  // Loop over atoms
  for (int atom = 0; atom < nAtoms; atom++) {
    int zself = id[atom];
    // Get neighbors and loop over them.  Simplest if only one element
    if (nelements == 1) {
      int nbat = nbatch[0][0]; // only element
      int size = BUFLEN - nbat;
      // WARNING: Assume we get all ghost neighbors to the real atoms!
      int n = nblist->GetNeighbors(atom, other[0][0] + nbat, rnb[0][0] + nbat,
                                   sqdist[0][0] + nbat, size);
      assert(size >= 0); // REMOVE LATER !!!
      for (i = nbat; i < nbat + n; i++)
        self[0][0][i] = atom;
      nbatch[0][0] += n;
    } else {
      int size = BUFLEN;
      int n = nblist->GetNeighbors(atom, other_buf, rnb_buf, sqdist_buf, size);
      assert(size >= 0); // REMOVE LATER !!!
      for (i = 0; i < n; i++) {
        int zother = id[other_buf[i]];
        int nbat = nbatch[zself][zother]++; // Count this atom
        self[zself][zother][nbat] = atom;
        other[zself][zother][nbat] = other_buf[i];
        rnb[zself][zother][nbat][0] = rnb_buf[i][0];
        rnb[zself][zother][nbat][1] = rnb_buf[i][1];
        rnb[zself][zother][nbat][2] = rnb_buf[i][2];
        sqdist[zself][zother][nbat] = sqdist_buf[i];
      }
    }
    // Now process any full batch
    for (int zo = 0; zo < nelements; zo++)
      if (nbatch[zself][zo] >= BUFLEN - maxnblen) {
        sigma_batch(self[zself][zo], other[zself][zo], rnb[zself][zo],
                    sqdist[zself][zo], zself, zo, nbatch[zself][zo],
                    calculatesigma2);
        nbatch[zself][zo] = 0;
      }
  } // Loop over atoms
  /* Process the remaining incomplete batches */
  for (int zs = 0; zs < nelements; zs++)
    for (int zo = 0; zo < nelements; zo++)
      if (nbatch[zs][zo])
        sigma_batch(self[zs][zo], other[zs][zo], rnb[zs][zo], sqdist[zs][zo],
                    zs, zo, nbatch[zs][zo], calculatesigma2);

  sigma2isvalid = calculatesigma2; /* Remember if it was calculated. */
  DEBUGPRINT;
}

void EMT::sigma_batch(int *self, int *other, Vec rnb[], double *sq_dist, int zs,
                      int zo, int n, int calculatesigma2) {
  int i;
  //  double dist[BUFLEN];
  //  double wght[BUFLEN];
  //  double arg[BUFLEN], ex1[BUFLEN];
  double dsigma1s[BUFLEN], dsigma2s[BUFLEN];
  double ds1o[BUFLEN], ds2o[BUFLEN];
  double *dsigma1o = 0, *dsigma2o = 0;
  double cutslopecutdist;
  double other_eta2betaseq, self_eta2betaseq;
  double other_kappaoverbeta, self_kappaoverbeta;
  double other_kappaseq, self_kappaseq;
  double *s1s, *s1o, *s2s, *s2o;
  const emt_parameters *emtself, *emtother;

  /* Get EMT parameters   REWRITE !!! XXXX */
  emtself = parameters[zs];
  emtother = parameters[zo];
  cutslopecutdist = cutoffslope * rFermi;
  other_eta2betaseq = emtother->eta2 * Beta * emtother->seq;
  self_eta2betaseq = emtself->eta2 * Beta * emtself->seq;
  other_kappaoverbeta = emtother->kappa / Beta;
  self_kappaoverbeta = emtself->kappa / Beta;
  other_kappaseq = emtother->kappa * emtother->seq;
  self_kappaseq = emtself->kappa * emtself->seq;
  double other_eta2 = emtother->eta2;
  double self_eta2 = emtself->eta2;

  assert(n <= BUFLEN);
  if ((zs == zo) && !calculatesigma2) {
    for (i = 0; i < n; i++) {
      /* dist = sqrt(sq_dist),  distances between atoms */
      double dist = sqrt(sq_dist[i]);
      /* Calculate weight factor (cutoff function) */
      double wght = 1.0 / (1.0 + exp(cutoffslope * dist - cutslopecutdist));
      /* Contribution to sigma1 */
      dsigma1s[i] = wght * exp(-other_eta2 * dist + other_eta2betaseq);
    }
    dsigma1o = dsigma1s;
  } else if ((zs != zo) && !calculatesigma2) {
    for (i = 0; i < n; i++) {
      /* dist = sqrt(sq_dist),  distances between atoms */
      double dist = sqrt(sq_dist[i]);
      /* Calculate weight factor (cutoff function) */
      double wght = 1.0 / (1.0 + exp(cutoffslope * dist - cutslopecutdist));
      /* Contributions to sigma1 */
      dsigma1s[i] = wght * exp(-other_eta2 * dist + other_eta2betaseq);
      ds1o[i] = wght * exp(-self_eta2 * dist + self_eta2betaseq);
    }
    dsigma1o = ds1o;
  } else if ((zs == zo) && calculatesigma2) {
    for (i = 0; i < n; i++) {
      /* dist = sqrt(sq_dist),  distances between atoms */
      double dist = sqrt(sq_dist[i]);
      /* Calculate weight factor (cutoff function) */
      double wght = 1.0 / (1.0 + exp(cutoffslope * dist - cutslopecutdist));
      /* Contribution to sigma1 */
      dsigma1s[i] = wght * exp(-other_eta2 * dist + other_eta2betaseq);
      dsigma2s[i] = wght * exp(-other_kappaoverbeta * dist + other_kappaseq);
    }
    dsigma1o = dsigma1s;
    dsigma2o = dsigma2s;
  } else {
    for (i = 0; i < n; i++) {
      /* dist = sqrt(sq_dist),  distances between atoms */
      double dist = sqrt(sq_dist[i]);
      /* Calculate weight factor (cutoff function) */
      double wght = 1.0 / (1.0 + exp(cutoffslope * dist - cutslopecutdist));
      /* Contributions to sigma1 */
      dsigma1s[i] = wght * exp(-other_eta2 * dist + other_eta2betaseq);
      ds1o[i] = wght * exp(-self_eta2 * dist + self_eta2betaseq);
      dsigma2s[i] = wght * exp(-other_kappaoverbeta * dist + other_kappaseq);
      ds2o[i] = wght * exp(-self_kappaoverbeta * dist + self_kappaseq);
    }
    dsigma1o = ds1o;
    dsigma2o = ds2o;
  }

  if (calculatesigma2) {
    /* Distribute contributions to sigma1 and sigma2 */
    s1o = &sigma1[zo][0];
    s1s = &sigma1[zs][0];
    s2o = &sigma2[zo][0];
    s2s = &sigma2[zs][0];

    for (i = 0; i < n; i++) {
      int s = self[i];
      int o = other[i];
      s1o[s] += dsigma1s[i];
      s2o[s] += dsigma2s[i];
      if (o < nAtoms) // Dont add to ghost atoms
      {
        s1s[o] += dsigma1o[i];
        s2s[o] += dsigma2o[i];
      }
    }
  } else {
    /* Distribute contributions to sigma1. */
    s1o = &sigma1[zo][0];
    s1s = &sigma1[zs][0];

    for (i = 0; i < n; i++) {
      int s = self[i];
      int o = other[i];
      s1o[s] += dsigma1s[i];
      if (o < nAtoms)
        s1s[o] += dsigma1o[i];
    }
  }
}

// Calculate the energies if Epot != 0, otherwise just calculate the
// derivatives needed for the forces.
void EMT::CalculateEnergiesAfterSigmas(double *Epot) {
  // USETIMER("EMT::CalculateEnergiesAfterSigmas");
  DEBUGPRINT;

  // double *sigma;
  // double *tmp1, *ex1, *ex2;
  int i;
  int zs, zo;
  double s;
  // Better performance if static ???
  assert(nelements < NMAXELEMENTS);
  double inv12gamma1[NMAXELEMENTS];
  double neginvbetaeta2[NMAXELEMENTS];
  double neglambda[NMAXELEMENTS];
  double lambdaseq[NMAXELEMENTS];
  double negkappa[NMAXELEMENTS];
  double kappaseq[NMAXELEMENTS];
  double nege0lambdalambda[NMAXELEMENTS];
  double e0lambdalambdaseq[NMAXELEMENTS];
  double neg6v0kappa[NMAXELEMENTS];
  double invbetaeta2[NMAXELEMENTS];
  double e0lambda[NMAXELEMENTS];
  double eccnst[NMAXELEMENTS];
  double sixv0[NMAXELEMENTS];
  double neghalfv0overgamma2[NMAXELEMENTS];
  double seq[NMAXELEMENTS];
  double *temporarydata, *sigma, *ex1, *ex2, *tmp1;
  int *id = &(this->id)[0];

  /* Calculate conbinations of EMT parameters */
  for (i = 0; i < nelements; i++) {
    inv12gamma1[i] = 1.0 / (12.0 * parameters[i]->gamma1);
    neginvbetaeta2[i] = -1.0 / (Beta * parameters[i]->eta2);
    neglambda[i] = -parameters[i]->lambda;
    lambdaseq[i] = parameters[i]->lambda * parameters[i]->seq;
    negkappa[i] = -parameters[i]->kappa;
    kappaseq[i] = parameters[i]->kappa * parameters[i]->seq;
    nege0lambdalambda[i] =
        -parameters[i]->e0 * parameters[i]->lambda * parameters[i]->lambda;
    e0lambdalambdaseq[i] = parameters[i]->e0 * parameters[i]->lambda *
                           parameters[i]->lambda * parameters[i]->seq;
    neg6v0kappa[i] = -6.0 * parameters[i]->V0 * parameters[i]->kappa;
    invbetaeta2[i] = 1.0 / (Beta * parameters[i]->eta2);
    e0lambda[i] = parameters[i]->e0 * parameters[i]->lambda;
    eccnst[i] =
        parameters[i]->e0 * (1.0 - parameters[i]->lambda * parameters[i]->seq);
    sixv0[i] = 6.0 * parameters[i]->V0;
    neghalfv0overgamma2[i] = -0.5 * parameters[i]->V0 / parameters[i]->gamma2;
    seq[i] = parameters[i]->seq;
  }

  /* Allocate temporary storage */
  //    sigma = dTmp[0];
  //    tmp1 = dTmp[1];
  //    ex1 = dTmp[2];
  //    ex2 = dTmp[3];
  temporarydata = new double[4 * nSize];
  sigma = temporarydata;
  tmp1 = sigma + nSize;
  ex1 = tmp1 + nSize;
  ex2 = ex1 + nSize;

  if (counters.beforeforces != atoms->GetChangeCounter() ||
      counters.energies != atoms->GetChangeCounter()) {
    counters.beforeforces = atoms->GetChangeCounter();
    VERB("b");
    DEBUGPRINT;
    /* Calculate total sigma1 */
    for (i = 0; i < nAtoms; i++) {
      s = 0.0;
      zs = id[i];
      for (zo = 0; zo < nelements; zo++)
        s += (*chi)[zs][zo] * sigma1[zo][i];
      if (s < 1.0e-9)
        s = 1.0e-9;
      sigma[i] = s;
    }
    if (ghostatoms)
      ghostatoms->GetGhostPotential()->CommunicateData(sigma);
    assert(nSize == radius.size() && nSize == Ec.size() &&
           nSize == dEds.size());

    /* radius = seq[z] - invbetaeta2[z] * log(sigma * inv12gamma1[z]) */
    /* Use ex1 as temporary variable */
    vec_mul_indir(tmp1, sigma, inv12gamma1, id, nSize);
    vlog(ex1, tmp1, &nSize);
    vec_mul_indir_add_indir(&radius[0], ex1, neginvbetaeta2, seq, id, nSize);

    /* dEds */
    /* ex1 = exp(neglambda[z] * radius + lambdaseq[z]) */
    vec_mul_indir_add_indir(tmp1, &radius[0], neglambda, lambdaseq, id, nSize);
    vexp(ex1, tmp1, &nSize);
    /* ex2 = exp(negkappa[z] * radius + kappaseq) */
    vec_mul_indir_add_indir(tmp1, &radius[0], negkappa, kappaseq, id, nSize);
    vexp(ex2, tmp1, &nSize);
    /* dEds = (nege0lambdalambda[z] * radius + e0labmdalabmdaseq[z]) * ex1
       + neg6v0kappa[z] * ex2 */

    vec_dEds(&dEds[0], &radius[0], ex1, ex2, nege0lambdalambda,
             e0lambdalambdaseq, neg6v0kappa, id, nSize);

    /* dEds = - dEds * invbetaeta2[z] / sigma */
    vrec(tmp1, sigma, &nSize);
    vec_self_mul_mul_indir(&dEds[0], tmp1, neginvbetaeta2, id, nSize);

    /* Cohesive energy */
    /* Ec = (e0lambda[z] * radius + eccnst) * ex1 */
    vec_mul_indir_add_indir_mul(&Ec[0], &radius[0], e0lambda, eccnst, ex1, id,
                                nSize);
  }
  if (Epot) {
    DEBUGPRINT;
    if (counters.energies != atoms->GetChangeCounter()) {
      counters.energies = atoms->GetChangeCounter();
      VERB("e");
      DEBUGPRINT;
      /* We also need Eas, but only for the real atoms */
      assert(sigma2isvalid);
      assert(counters.sigma2 == atoms->GetChangeCounter());
      /* Calculate total sigma2 */
      for (i = 0; i < nAtoms; i++) {
        s = 0.0;
        zs = id[i];
        for (zo = 0; zo < nelements; zo++)
          s += (*chi)[zs][zo] * sigma2[zo][i];
        sigma[i] = s;
      }

      /* Atomic-sphere energy */
      /* Eas = sixv0[z] * ex2 + neghalfv0overgamma2 * sigma */
      vec_dbl_mul_indir_add(&Eas[0], ex2, sixv0, sigma, neghalfv0overgamma2, id,
                            nAtoms);
    }
    /* Add the energies */
    // Not defined: vec_add(Epot, Ec, Eas, nAtoms);
    DEBUGPRINT;

    if (subtractE0)
      for (i = 0; i < nAtoms; i++)
        Epot[i] = Ec[i] + Eas[i] - parameters[id[i]]->e0;
    else
      for (int i = 0; i < nAtoms; i++)
        Epot[i] = Ec[i] + Eas[i];

  } // if (Epot)

  delete[] temporarydata;

  DEBUGPRINT;
}

void EMT::CalculateForcesAfterEnergies(Vec forces[], double stresses[][6]) {
  int i;
  // USETIMER("EMT::CalculateForcesAfterEnergies");
  DEBUGPRINT;
  if ((!forces || (counters.forces == atoms->GetChangeCounter())) &&
      (!stresses || (counters.stresses == atoms->GetChangeCounter())))
    return;
  if (forces) {
    counters.forces = atoms->GetChangeCounter();
    VERB("f");
  }
  if (stresses) {
    counters.stresses = atoms->GetChangeCounter();
    VERB("s");
  }

  DEBUGPRINT;
  int maxnblen = nblist->MaxNeighborListLength();
  // Buffer data:
  TinyMatrix<int> nbatch(nelements, nelements);
  TinyMatrix<int[BUFLEN]> self(nelements, nelements);
  TinyMatrix<int[BUFLEN]> other(nelements, nelements);
  TinyMatrix<Vec[BUFLEN]> rnb(nelements, nelements);
  TinyMatrix<double[BUFLEN]> sqdist(nelements, nelements);
  TinyMatrix<double[BUFLEN]> dEdss(nelements, nelements);
  TinyMatrix<double[BUFLEN]> dEdso(nelements, nelements);
  int other_buf[BUFLEN];
  Vec rnb_buf[BUFLEN];
  double sqdist_buf[BUFLEN];

  // If there is only one element, CalculateForcesAfterEnergiesSingle should
  // be used.
  assert(nelements > 1);

  /* Set forces and stresses to zero */
  if (forces)
    memset(&forces[0][0], 0, nAtoms * sizeof(Vec));
  if (stresses)
    memset(stresses, 0, nAtoms * 6 * sizeof(double));

  /* Calculate forces */

  /* No atoms in batch pools */
  for (i = 0; i < nelements; i++)
    for (int j = 0; j < nelements; j++)
      nbatch[i][j] = 0;

  // Loop over atoms
  for (int atom = 0; atom < nAtoms; atom++) {
    int zself = id[atom];
    // Get neighbors and loop over them.
    int size = BUFLEN;
    int n = nblist->GetNeighbors(atom, other_buf, rnb_buf, sqdist_buf, size);
    assert(size >= 0); // REMOVE LATER
    for (i = 0; i < n; i++) {
      int zother = id[other_buf[i]];
      int nbat = nbatch[zself][zother]++; // Count this atom
      self[zself][zother][nbat] = atom;
      other[zself][zother][nbat] = other_buf[i];
      rnb[zself][zother][nbat][0] = rnb_buf[i][0];
      rnb[zself][zother][nbat][1] = rnb_buf[i][1];
      rnb[zself][zother][nbat][2] = rnb_buf[i][2];
      sqdist[zself][zother][nbat] = sqdist_buf[i];
      dEdss[zself][zother][nbat] = dEds[atom];
      dEdso[zself][zother][nbat] = dEds[other_buf[i]];
    }
    // Now process any full batch
    for (int zo = 0; zo < nelements; zo++)
      if (nbatch[zself][zo] >= BUFLEN - maxnblen) {
        force_batch(self[zself][zo], other[zself][zo], rnb[zself][zo],
                    sqdist[zself][zo], dEdss[zself][zo], dEdso[zself][zo],
                    zself, zo, nbatch[zself][zo], forces, stresses);
        nbatch[zself][zo] = 0;
      }
  } // Loop over atoms
  /* Process the remaining incomplete batches */
  for (int zs = 0; zs < nelements; zs++)
    for (int zo = 0; zo < nelements; zo++)
      if (nbatch[zs][zo])
        force_batch(self[zs][zo], other[zs][zo], rnb[zs][zo], sqdist[zs][zo],
                    dEdss[zs][zo], dEdso[zs][zo], zs, zo, nbatch[zs][zo],
                    forces, stresses);
  DEBUGPRINT;
}

void EMT::CalculateForcesAfterEnergiesSingle(Vec forces[],
                                             double stresses[][6]) {
  int i;
  // USETIMER("EMT::CalculateForcesAfterEnergiesSingle");
  DEBUGPRINT;
  if ((!forces || (counters.forces == atoms->GetChangeCounter())) &&
      (!stresses || (counters.stresses == atoms->GetChangeCounter())))
    return;
  if (forces) {
    counters.forces = atoms->GetChangeCounter();
    VERB("f");
  }
  if (stresses) {
    counters.stresses = atoms->GetChangeCounter();
    VERB("s");
  }

  int maxnblen = nblist->MaxNeighborListLength();
  DEBUGPRINT;
  // Buffer data:
  int nbatch;
  int self[BUFLEN];
  int other[BUFLEN];
  Vec rnb[BUFLEN];
  double sqdist[BUFLEN];
  double dEdss[BUFLEN];
  double dEdso[BUFLEN];

  DEBUGPRINT;
  // If there is more than one element, CalculateForcesAfterEnergies should
  // be used.
  assert(nelements == 1);

  DEBUGPRINT;
  /* Set forces and stresses to zero */
  if (forces)
    memset(&forces[0][0], 0, nAtoms * sizeof(Vec));
  if (stresses)
    memset(stresses, 0, nAtoms * 6 * sizeof(double));

  /* Calculate forces */

  /* No atoms in batch pool */
  DEBUGPRINT;
  nbatch = 0;

  // Loop over atoms
  DEBUGPRINT;
  for (int atom = 0; atom < nAtoms; atom++) {
    DEBUGPRINT;

    // Get neighbors and loop over them.
    int size = BUFLEN - nbatch;
    int n = nblist->GetNeighbors(atom, other + nbatch, rnb + nbatch,
                                 sqdist + nbatch, size);
    DEBUGPRINT;
    double dEdsatom = dEds[atom];
    for (i = nbatch; i < nbatch + n; i++) {
      self[i] = atom;
      dEdss[i] = dEdsatom;
      dEdso[i] = dEds[other[i]];
    }
    nbatch += n;
    // Now process any full batch
    if (nbatch >= BUFLEN - maxnblen) {
      DEBUGPRINT;
      force_batch(self, other, rnb, sqdist, dEdss, dEdso, 0, 0, nbatch, forces,
                  stresses);
      DEBUGPRINT;
      nbatch = 0;
    }
  } // Loop over atoms
  /* Process the remaining incomplete batches */
  DEBUGPRINT;
  if (nbatch)
    force_batch(self, other, rnb, sqdist, dEdss, dEdso, 0, 0, nbatch, forces,
                stresses);
  DEBUGPRINT;
}

void EMT::force_batch(int *self, int *other, Vec rnb[BUFLEN], double *sq_dist,
                      double dEdss[BUFLEN], double dEdso[BUFLEN], int zs,
                      int zo, int n, Vec *force, double (*stresses)[6]) {
  int i;
  double df[BUFLEN];
  double cutslopecutdist, other_eta2betaseq, other_kappaoverbeta;
  double other_kappaseq, self_eta2betaseq, self_kappaoverbeta;
  double self_kappaseq, cnst_s, cnst_o;
  const emt_parameters *emtself, *emtother;
  // double pairA, exprcut, pairD;
  // double *tmp;

  /* Get EMT parameters */
  emtself = parameters[zs];
  emtother = parameters[zo];
  cutslopecutdist = cutoffslope * rFermi;
  other_eta2betaseq = emtother->eta2 * Beta * emtother->seq;
  other_kappaoverbeta = emtother->kappa / Beta;
  other_kappaseq = emtother->kappa * emtother->seq;
  self_eta2betaseq = emtself->eta2 * Beta * emtself->seq;
  self_kappaoverbeta = emtself->kappa / Beta;
  self_kappaseq = emtself->kappa * emtself->seq;
  double other_eta2 = emtother->eta2;
  double self_eta2 = emtself->eta2;

  cnst_s = -0.5 * emtself->V0 * (*chi)[zs][zo] / emtself->gamma2;
  cnst_o = -0.5 * emtother->V0 * (*chi)[zo][zs] / emtother->gamma2;
  double chi_zs_zo = (*chi)[zs][zo];
  double chi_zo_zs = (*chi)[zo][zs];
  if (zs == zo) {
    for (i = 0; i < n; i++) {
      /* Get the distances from their squares */
      double dist = sqrt(sq_dist[i]);
      double inv_dist = 1.0 / dist;
      /* Calculate cutoff function (weight factor) */
      double wght = 1.0 / (1.0 + exp(cutoffslope * dist - cutslopecutdist));
      /* Calculate derivative of the cutoff function. */
      double dwdr = -cutoffslope * wght * (1 - wght);
      double dsigma1dr = (dwdr - wght * other_eta2) *
                         exp(-other_eta2 * dist + other_eta2betaseq);
      double dsigma2dr = (dwdr + wght * -other_kappaoverbeta) *
                         exp(-other_kappaoverbeta * dist + other_kappaseq);
      df[i] =
          inv_dist * (dsigma1dr * dEdss[i] * chi_zs_zo + cnst_s * dsigma2dr +
                      dsigma1dr * dEdso[i] * chi_zo_zs + cnst_o * dsigma2dr);
    }
  } else {
    for (i = 0; i < n; i++) {
      /* Get the distances from their squares */
      double dist = sqrt(sq_dist[i]);
      double inv_dist = 1.0 / dist;
      /* Calculate cutoff function (weight factor) */
      double wght = 1.0 / (1.0 + exp(cutoffslope * dist - cutslopecutdist));
      /* Calculate derivative of the cutoff function. */
      double dwdr = -cutoffslope * wght * (1 - wght);
      double dsigma1dr_o = (dwdr - wght * other_eta2) *
                           exp(-other_eta2 * dist + other_eta2betaseq);
      double dsigma2dr_o = (dwdr + wght * -other_kappaoverbeta) *
                           exp(-other_kappaoverbeta * dist + other_kappaseq);
      double dsigma1dr_s =
          (dwdr - wght * self_eta2) * exp(-self_eta2 * dist + self_eta2betaseq);
      double dsigma2dr_s = (dwdr + wght * -self_kappaoverbeta) *
                           exp(-self_kappaoverbeta * dist + self_kappaseq);
      df[i] = inv_dist *
              (dsigma1dr_o * dEdss[i] * chi_zs_zo + cnst_s * dsigma2dr_o +
               dsigma1dr_s * dEdso[i] * chi_zo_zs + cnst_o * dsigma2dr_s);
    }
  }

  /* Distribute force contributions */
  if (force) {
    for (i = 0; i < n; i++)
      for (int j = 0; j < 3; j++) {
        int o = other[i];
        double dfx = df[i] * rnb[i][j];
        force[self[i]][j] += dfx;
        if (o < nAtoms)
          force[o][j] -= dfx;
      }
  }
  /* Distribute stress contributions */
  if (stresses) {
    for (i = 0; i < n; i++)
      for (int alpha = 0; alpha < 3; alpha++)
        for (int beta = alpha; beta < 3; beta++) {
          int o = other[i];
          double dsx = 0.5 * df[i] * rnb[i][alpha] * rnb[i][beta];
          int ii = stresscomp[alpha][beta];
          stresses[self[i]][ii] += dsx;
          if (o < nAtoms)
            stresses[o][ii] += dsx;
        }
  }
}

void EMT::UpdateSuperCell(const SuperCell *newSuperCell) {
  DEBUGPRINT;
  if (nblist != 0)
    nblist->UpdateSuperCell(newSuperCell);
  DEBUGPRINT;
}

// Quasicontinuum stuff

void EMT::SetContinuumElement(int z) {
  int i;
  static int shellPop[6] = {1, 12, 6, 24, 12, 24};
  int nShells = 3;
  int indices[][3] = {
      {1, 0, 0},  {0, 1, 0},  {-1, 1, 0}, {0, 0, 1},  {-1, 0, 1}, {0, -1, 1},
      {1, 1, -1}, {1, -1, 1}, {-1, 1, 1}, {1, 1, 0},  {1, 0, 1},  {0, 1, 1},
      {2, -1, 0}, {0, 2, -1}, {-1, 0, 2}, {2, 0, -1}, {-1, 2, 0}, {0, -1, 2},
      {-2, 1, 1}, {1, -2, 1}, {1, 1, -2}, {2, 0, 0},  {0, 2, 0},  {0, 0, 2},
      {0, 2, -2}, {2, 0, -2}, {2, -2, 0}};

  nHalfNeighbors = 0;

  continuumelement = 0;
  if (z) {
    vector<const emt_parameters *>::iterator param;
    for (param = parameters.begin(); param != parameters.end(); ++param) {
      if ((*param)->Z == z)
        continuumelement = *param;
    }
    assert(continuumelement != 0);
  }

  for (int shell = 1; shell <= nShells; shell++)
    nHalfNeighbors += shellPop[shell] / 2;
  assert(nHalfNeighbors <= 27);
  if (coef)
    delete[] coef;
  coef = new int[nHalfNeighbors][6];
  for (int n = 0; n < nHalfNeighbors; n++) {
    int k = 0;
    for (i = 0; i < 3; i++) {
      for (int j = i; j < 3; j++) {
        coef[n][k] = indices[n][i] * indices[n][j];
        if (j > i)
          coef[n][k] *= 2;
        k++;
      }
    }
  }
}

double EMT::CalculateLatticeEnergy(const Vec a[3]) {
  // USETIMER("EMT::CalculateLatticeEnergy");
  double sigma1;
  double sigma2;
  double diffs2[50];
  double wght[50];
  double dsigma1[50];
  double dsigma2[50];
  CalculateStuff(a, sigma1, sigma2, diffs2, wght, dsigma1, dsigma2);
  double tmp = sigma1 / (12.0 * continuumelement->gamma1);
  double radiusMinusSeq = -log(tmp) / (Beta * continuumelement->eta2);
  double ex1 = exp(-continuumelement->kappa * radiusMinusSeq);
  tmp = continuumelement->lambda * radiusMinusSeq;
  double Ec = continuumelement->e0 * (tmp + 1.0) * exp(-tmp);
  double EAS = continuumelement->V0 *
               (6.0 * ex1 - 0.5 * sigma2 / continuumelement->gamma2);
  return Ec + EAS - continuumelement->e0;
}

void EMT::CalculateDerivatives(const Vec a[3], double dEdaDota[6]) {
  double sigma1;
  double sigma2;
  double diffs2[50];
  double wght[50];
  double dsigma1[50];
  double dsigma2[50];
  int k;
  // USETIMER("EMT::CalculateDerivatives");
  CalculateStuff(a, sigma1, sigma2, diffs2, wght, dsigma1, dsigma2);
  vrec(diffs2, diffs2, &nHalfNeighbors);
  vec_mul_scal_add_scal(wght, wght, cutoffslope,
                        -cutoffslope - continuumelement->eta2, nHalfNeighbors);
  vec_mul(dsigma1, dsigma1, wght, nHalfNeighbors);
  vec_mul(dsigma1, dsigma1, diffs2, nHalfNeighbors);
  vec_mul_scal_add_scal(wght, wght, 1.0,
                        continuumelement->eta2 - continuumelement->kappa / Beta,
                        nHalfNeighbors);
  vec_mul(dsigma2, dsigma2, wght, nHalfNeighbors);
  vec_mul(dsigma2, dsigma2, diffs2, nHalfNeighbors);
  double dsigma1daDota[6];
  double dsigma2daDota[6];
  for (k = 0; k < 6; k++) {
    double d1 = 0.0;
    double d2 = 0.0;
    for (int n = 0; n < nHalfNeighbors; n++) {
      d1 += coef[n][k] * dsigma1[n];
      d2 += coef[n][k] * dsigma2[n];
    }
    dsigma1daDota[k] = d1;
    dsigma2daDota[k] = d2;
  }
  double tmp = sigma1 / (12.0 * continuumelement->gamma1);
  double radiusMinusSeq = -log(tmp) / (Beta * continuumelement->eta2);
  double x = continuumelement->lambda * radiusMinusSeq;
  double dEdsigma1 =
      (6.0 * continuumelement->V0 * continuumelement->kappa *
           exp(-continuumelement->kappa * radiusMinusSeq) +
       continuumelement->e0 * continuumelement->lambda * x * exp(-x)) /
      (Beta * continuumelement->eta2 * sigma1);
  double dEdsigma2 = -0.5 * continuumelement->V0 / continuumelement->gamma2;
  for (k = 0; k < 6; k++)
    dEdaDota[k] = dEdsigma1 * dsigma1daDota[k] + dEdsigma2 * dsigma2daDota[k];
}

void EMT::CalculateStuff(const Vec a[3], double &sigma1, double &sigma2,
                         double *diffs2, double *wght, double *dsigma1,
                         double *dsigma2) {
  int i;
  // USETIMER("EMT::CalculateStuff");
  double aDota[6];
  int n, k = 0;
  if (!continuumelement) {
    cerr << "Use EMT::SetContinuumElement() to choose element in the"
         << " QC part" << endl;
    exit(1);
  }

  for (i = 0; i < 3; i++)
    for (int j = i; j < 3; j++) {
      double dot = 0.0;
      for (int coord = 0; coord < 3; coord++)
        dot += a[i][coord] * a[j][coord];
      aDota[k++] = dot;
    }
  for (n = 0; n < nHalfNeighbors; n++) {
    double r2 = 0.0;
    for (int k = 0; k < 6; k++)
      r2 += coef[n][k] * aDota[k];
    diffs2[n] = r2;
  }
  // diffs2 = sqrt(diffs2), distance between atoms
  vsqrt(diffs2, diffs2, &nHalfNeighbors);
  // dsigma1 = exp(-emeta2 * dist + emeta2BetaSeq)
  vec_mul_scal_add_scal(wght, diffs2, -continuumelement->eta2,
                        continuumelement->eta2 * Beta * continuumelement->seq,
                        nHalfNeighbors);
  vexp(dsigma1, wght, &nHalfNeighbors);
  // dsigma2 = exp(-emKappaOverBeta * dist + emKappaSeq)
  vec_mul_scal_add_scal(wght, diffs2, -continuumelement->kappa / Beta,
                        continuumelement->kappa * continuumelement->seq,
                        nHalfNeighbors);
  vexp(dsigma2, wght, &nHalfNeighbors);
  bool veryAccurate = true;
  if (veryAccurate) {
    // wght = exp(cutoffslope * dist - cutslopeCutdist)
    vec_mul_scal_add_scal(wght, diffs2, cutoffslope, -cutoffslope * rFermi,
                          nHalfNeighbors);
    vexp(wght, wght, &nHalfNeighbors);
    // wght = 1 / (wght + 1.0)
    vec_mul_scal_add_scal(wght, wght, 1.0, 1.0, nHalfNeighbors);
    vrec(wght, wght, &nHalfNeighbors);
    // dsigmai = wght * dsigmai
    vec_mul(dsigma1, wght, dsigma1, nHalfNeighbors);
    vec_mul(dsigma2, wght, dsigma2, nHalfNeighbors);
  }
  sigma1 = 0.0;
  sigma2 = 0.0;
  for (n = 0; n < nHalfNeighbors; n++) {
    sigma1 += dsigma1[n];
    sigma2 += dsigma2[n];
  }
  sigma1 *= 2;
  sigma2 *= 2;
  dataSigma1 = sigma1;
}

double EMT::GetLatticeConstant() const {
  assert(continuumelement != 0);
  return Beta * continuumelement->seq * sqrt(2.0);
}

void EMT::PrintParameters() {
  int i;
  for (i = 0; i < nelements; i++) {
    const emt_parameters *p = parameters[i];
    cerr << "Parameters for element " << i << " (" << p->name << ")" << endl;
    cerr << "E0:" << p->e0 << "  s0:" << p->seq << "  V0:" << p->V0
         << "  eta2:" << p->eta2 << "  kappa:" << p->kappa
         << "  lambda:" << p->lambda << "  rFermi:" << rFermi << "  cutSlope"
         << cutoffslope << "  gamma1:" << p->gamma1 << "  gamma2:" << p->gamma2
         << endl
         << endl;
  }
}
