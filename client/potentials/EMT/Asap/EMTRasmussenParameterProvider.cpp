// Copyright (C) 2008-2011 Jakob Schiotz and Center for Individual
// Nanoparticle Functionality, Department of Physics, Technical
// University of Denmark.  Email: schiotz@fysik.dtu.dk
//
// This file is part of Asap version 3.
//
// This program is free software: you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public License
// version 3 as published by the Free Software Foundation.  Permission
// to use other versions of the GNU Lesser General Public License may
// granted by Jakob Schiotz or the head of department of the
// Department of Physics, Technical University of Denmark, as
// described in section 14 of the GNU General Public License.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// and the GNU Lesser Public License along with this program.  If not,
// see <http://www.gnu.org/licenses/>.

#include "EMTRasmussenParameterProvider.h"
#include "Exception.h"
#include <cmath>
#include <iostream>
#include <string.h>

using namespace std;

emt_parameters *EMTRasmussenParameterProvider::GetNewParameters(int element) {
  double E0, S0, n0, V0, eta2, kappa, lambda, mass;
  double ls;
  int Z;
  string name;
  emt_parameters *p;

  if (element == 13) {
    name = "Al";
    E0 = -3.280;
    S0 = 3.000;
    V0 = 1.725;
    Z = 13;
    eta2 = 1.380;
    kappa = 2.311;
    lambda = 1.591;
    mass = 26.98;
    n0 = 0.007 /* CHECK n0 */;
    ls = 7.54871784;
  } else if (element == 29) {
    name = "Cu";
    //    E0=-3.510; S0=2.67 ; V0=3.010; Z=29; eta2=1.494; kappa=2.500;
    //    lambda=1.942; mass=63.54 ; n0=0.0091 /* CHECK n0 */ ;
    E0 = -3.510;
    S0 = 2.67;
    V0 = 2.643;
    Z = 29;
    eta2 = 1.506;
    kappa = 2.492;
    lambda = 1.942;
    mass = 63.54;
    n0 = 0.0091 /* CHECK n0 */;
    ls = 6.789382809;
  } else if (element == 47) {
    name = "Ag";
    E0 = -2.96;
    S0 = 3.01;
    V0 = 2.679;
    Z = 47;
    eta2 = 1.400;
    kappa = 2.365;
    lambda = 1.956;
    mass = 107.87;
    n0 = 0.0059;
    ls = 7.6790043;
  } else if (element == 79) {
    name = "Au";
    E0 = -3.80;
    S0 = 3.00;
    V0 = 2.703;
    Z = 79;
    eta2 = 1.310;
    kappa = 2.221;
    lambda = 2.192;
    mass = 196.97;
    n0 = 0.0064;
    ls = 7.66504117182;
  }
#if 0
  else if(element == 28){
    name = "Ni";
    E0=-4.44 ; S0=2.60 ; V0=3.673; Z=28; eta2=1.669; kappa=2.757;
    lambda=1.948; mass=58.71 ; n0=0.0103 ;
    ls = 6.598896;
  }
  else if(element == 46){
    name = "Pd";
    E0=-3.90 ; S0=2.87 ; V0=2.773; Z=46; eta2=1.818; kappa=3.107;
    lambda=2.155; mass=106.4 ; n0=0.00688;
    ls = 7.330378;
  }
  else if(element == 78){
    name = "Pt";
    E0=-5.85 ; S0=2.90 ; V0=4.067; Z=78; eta2=1.812; kappa=3.145;
    lambda=2.192; mass=195.09; n0=0.00802;
    ls = 7.41119853;
  }
#endif
  else {
    throw Exception("This element isn't defined in EMT.");
  }
  //  else
  //      throw AsapError("Unknown element Z = ") << element;

  double bohr = 0.5291772; // Angstrom

  p = new emt_parameters;
  p->e0 = E0;
  p->seq = S0 * bohr;
  p->neq = n0 / (bohr * bohr * bohr);
  p->V0 = V0;
  p->eta2 = eta2 / bohr;
  p->kappa = kappa / bohr;
  p->lambda = lambda / bohr;
  p->mass = mass;
  p->invmass = 1.0 / mass;
  p->gamma1 = 0.0; // These are calculated later!
  p->gamma2 = 0.0;
  p->Z = Z;
  //  assert(element == Z);
  p->name = name;
  p->lengthscale = ls / sqrt(2.0) * bohr;

  return p;
}
