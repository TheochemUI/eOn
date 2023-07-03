#pragma once
#ifndef FORCEFIELDS_UNIT_SYSTEM_HPP
#define FORCEFIELDS_UNIT_SYSTEM_HPP

/** @file
Physical constants, unit conversion.
@author Jean-Claude C. Berthet
@date 2006-2007
University of Iceland
@section introduction Introduction
This file contains constants to convert between different units and also to
choose a default system of units. There are two types of constants, those in @em
lowercase letters and those in @em capital letters. Constants in lowercase
letters are universal constants usually expressed in the international sytem of
unit (SI). The contents of these constants do not depend on the default system
of unit chosen. The content of constants is capital letters depends on the
default system. Those latter are further devided into two subgroups: unit name
constants and converters.

@section Default system of units:
The default system of units is selected at compilation by setting the macro
#FORCEFIELDS_UNIT_SYSTEM_HPP .

The content of the constants called @ref units "unit name constants" and @ref
converters "converters" depends on the system chosen. \n The default system like
the SI is based on the four dimensions distance, mass, time, electric current.
All other units are extracted from the four units chosen for these four
dimensions. For example if the default system is SI, the unit of energy will be
Joule because @f$ 1\ J = 1\ m^2\ kg\ s^{-2}  @f$ . As well the unit of pressure
will be Pascal ( @f$ 1\ Pa=1\ m^{-1}\ kg\ s^{-2} @f$ ). Quantities of same
dimensions are always expressed with the same units. For instance the
temperature has the dimension of the energy, therefore is default unit is the
unit of the energy (Joule if SI is the default system).

@section units Default unit names:
A unit constant is a constant string of characters named after a quantity (e.g.
ENERGY, PRESSURE, etc). It contains the name of the default unit used for that
quantity (e.g. @em kJ , @em eV , etc for an energy). These constants can be used
to display a quantity in the default system of unit. For example, let us assume
that <em> kJ/mol </em> is the default unit for the energy. We wish to display
the content of a variable @em energy which contains an energy expressed in the
default system of unit. We shall write the following code:
@code
cout << energy << ' ' << ENERGY << endl;
@endcode

@section converters Converters: converting quantities.
A converter is a double constant which has the name of a unit (e.g. PASCAL,
KELVIN, ANGSTROM...). The constant contains the value of one quantity of the
unit espressed in the default system of unit. Here is an example on how to use
the converters. Let us assume that the default system of unit uses the meter as
length unit. We wish to display the content of a variable @a x which contains a
distance in meter, but we want to display it in Angstrom. To make the conversion
we need to divide the quantity by one Angstrom:
@f[
    \frac{x\ meter}{1\ Angstrom}=\frac{x\ meter}{10^{-10}\ meter}=x \time
10^{-10}
    @f]
We shall add to the code the following line which will display the quantity in
Angstrom.
@code
cout << "x = " << x / ANGSTROM << " Angstrom." << endl;
@endcode
Reversively, we may also wish to input a value into @em x. We shall write for
example:
@code
double x_angstrom;
cout << "Enter distance in Angstrom: ";
cin >> x_angstrom;
x= x_angstrom*ANGSTROM;
cout << "The distance is " << x/ANGSTROM << " Angstrom." << endl;
@endcode
*/

/// Systeme international
#define FORCEFIELDS_UNIT_SYSTEM_SI_MKSA 0
/// System based on Metre, kg/mol, second, Ampere/mol
#define FORCEFIELDS_UNIT_SYSTEM_METRE_KILOGRAMME_PER_MOL_SECOND_AMPERE_PER_MOL 1
/// System based on angstrom, kg/mol, femtosecond, Ampere/mol
#define FORCEFIELDS_UNIT_SYSTEM_ANGSTROM_KILOGRAMME_PER_MOL_FEMTOSECOND_AMPERE_PER_MOL \
  2
/// Angstrom, gramme/mol, femtosecond, e (proton charge)
#define FORCEFIELDS_UNIT_SYSTEM_ANGSTROM_GRAMME_PER_MOL_FEMTOSECOND_ECHARGE 3
/// Electron volt, angstrom, femtosecond, e (proton charge)
#define FORCEFIELDS_UNIT_SYSTEM_ELECTRONVOLT_ANGSTROM_FEMTOSECOND_ECHARGE 4
/// Kilo calory per mole, angstrom, femtosecond, e (proton charge)
#define FORCEFIELDS_UNIT_SYSTEM_KILOCALORY_PER_MOL_ANGSTROM_FEMTOSECOND_ECHARGE \
  5

/// Default system of units.
#undef FORCEFIELDS_UNIT_SYSTEM_HPP
#define FORCEFIELDS_UNIT_SYSTEM_HPP                                            \
  FORCEFIELDS_UNIT_SYSTEM_ELECTRONVOLT_ANGSTROM_FEMTOSECOND_ECHARGE
namespace forcefields {
/// Physical constants and unit conversion
namespace unit_system {
// Some physical constants expressed in SI unit. Source for constants:
// http://physics.nist.gov/cuu/Constants/index.html
const double R_gas = 8.314472;            // J/mol/K
const double k_boltzmann = 1.3806505e-23; // J/K
const double Na_avogadro = 6.0221415e23;  // mol^-1
const double /* 1 */ calory = 4.1868;     // J
const double e_charge = 1.60217653e-19;   // C
const double /* 1 */ debye =
    3.33564e-30; // C m
                 /// This is the quantity 1/(4*pi*epsilon0) in SI units.
const double one_over_4_pi_epsilon0 =
    8.98755178495272255e+09; // SI
                             /// This is the quantity 1/sqrt(4*pi*epsilon0) in
                             /// SI units.
const double sqrt_one_over_4_pi_epsilon0 = 9.48026992492973368e+04; // SI

// System independent units
const double pi = 3.14159265358979312;
const double MOL = Na_avogadro;
const double DEGREE = pi / 180; // rad

#if FORCEFIELDS_UNIT_SYSTEM_HPP == FORCEFIELDS_UNIT_SYSTEM_SI_MKSA
const char LENGTH[] = "m";
const char MASS[] = "kg";
const char TIME[] = "s";
const char CURRENT[] = "A"; // Ampere
const char VELOCITY[] = "m/s";
const char ACCELERATION[] = "m s^-2";
const char FORCE[] = "N";
const char ENERGY[] = "J";
const char PRESSURE[] = "Pa";
const char CHARGE[] = "C";

// Main SI units Converters
const double /* 1 */ METRE = 1.0;  // m
const double /* 1 */ KG = 1.0;     // kg
const double /* 1 */ SECOND = 1.0; // s
const double /* 1 */ AMPERE = 1.0; // A
#elif FORCEFIELDS_UNIT_SYSTEM_HPP ==                                           \
    FORCEFIELDS_UNIT_SYSTEM_METRE_KILOGRAMME_PER_MOL_SECOND_AMPERE_PER_MOL
const char LENGTH[] = "m";
const char MASS[] = "kg/mol";
const char TIME[] = "s";
const char CURRENT[] = "A/mol";
const char VELOCITY[] = "m/s";
const char ACCELERATION[] = "m s^-2";
const char FORCE[] = "N/mol";
const char ENERGY[] = "J/mol";
const char PRESSURE[] = "Pa/mol";
const char CHARGE[] = "C/mol";

// Main SI units Converters
const double /* 1 */ METRE = 1.0;  // m
const double /* 1 */ KG = MOL;     // kg/mol
const double /* 1 */ SECOND = 1.0; // s
const double /* 1 */ AMPERE = MOL; // A/mol
#elif FORCEFIELDS_UNIT_SYSTEM_HPP ==                                           \
    FORCEFIELDS_UNIT_SYSTEM_ANGSTROM_KILOGRAMME_PER_MOL_FEMTOSECOND_AMPERE_PER_MOL
const char LENGTH[] = "Angstrom";
const char MASS[] = "kg/mol";
const char TIME[] = "fs";
const char CURRENT[] = "A/mol";
const char VELOCITY[] = "Angstrom/fs";
const char ACCELERATION[] = "Angstrom fs^-2";
const char FORCE[] = "1e30 N/mol";
const char ENERGY[] = "1e20 J/mol";
const char PRESSURE[] = "1e50 Pa/mol";
const char CHARGE[] = "C/mol";

// Main SI units Converters
const double /* 1 */ METRE = 1e10;  // m
const double /* 1 */ KG = MOL;      // kg/mol
const double /* 1 */ SECOND = 1e15; // s
const double /* 1 */ AMPERE = MOL;  // A/mol
#elif FORCEFIELDS_UNIT_SYSTEM_HPP ==                                           \
    FORCEFIELDS_UNIT_SYSTEM_ANGSTROM_GRAMME_PER_MOL_FEMTOSECOND_ECHARGE
const char LENGTH[] = "Angstrom";
const char MASS[] = "g/mol";
const char TIME[] = "fs";
const char CURRENT[] = "e/fs";
const char VELOCITY[] = "Angstrom/fs";
const char ACCELERATION[] = "Angstrom fs^-2";
const char FORCE[] = "1e17 J / Angstrom";
const char ENERGY[] = "1e17 J";
const char PRESSURE[] = "1e17 J Angstrom^-3";
const char CHARGE[] = "e";

// Main SI units Converters
const double /* 1 */ METRE = 1e10;    // Angstrom
const double /* 1 */ KG = 1000 * MOL; // g/mol
const double /* 1 */ SECOND =
    1e15; // fs
          // 1 e/fs= e_charge C / 1e-15 s=e_charge*1e15 A
const double /* 1 */ AMPERE = 1e-15 / e_charge; // e/fs
#elif FORCEFIELDS_UNIT_SYSTEM_HPP ==                                           \
    FORCEFIELDS_UNIT_SYSTEM_ELECTRONVOLT_ANGSTROM_FEMTOSECOND_ECHARGE
const char LENGTH[] = "A"; // Angstrom
const char TIME[] = "fs";  // second
const char MASS[] = "eV A^-2 fs^2";
const char CURRENT[] = "e/fs";
const char VELOCITY[] = "A fs^-1";
const char ACCELERATION[] = "A fs^-2";
const char FORCE[] = "eV/A";
const char ENERGY[] = "eV"; // electron volt
const char PRESSURE[] = "eV A^-3";
const char CHARGE[] = "(eV A)^(-1/2)";

// Main SI units Converters
const double /* 1 */ METRE = 1e10;         // Angstrom
                                           //  1 J = 1 kg m^2 s^-2
                                           //  e_charge C V = e_charge J = 1 eV
const double /* 1 */ KG = 1e10 / e_charge; // eV Angstrom^-2 fs^2
const double /* 1 */ SECOND = 1e15;        // fs
const double /* 1 */ AMPERE = 1e-15 / e_charge; // e/fs
#elif FORCEFIELDS_UNIT_SYSTEM_HPP ==                                           \
    FORCEFIELDS_UNIT_SYSTEM_KILOCALORY_PER_MOL_ANGSTROM_FEMTOSECOND_ECHARGE
const char LENGTH[] = "A"; // Angstrom
const char TIME[] = "fs";  // second
const char MASS[] = "kcal/mol A^-2 fs^2";
const char CURRENT[] = "e/fs";
const char VELOCITY[] = "A fs^-1";
const char ACCELERATION[] = "A fs^-2";
const char FORCE[] = "kcal mol A^-1";
const char ENERGY[] = "kcal/mol"; // electron volt
const char PRESSURE[] = "kcal/mol A^-3";
const char CHARGE[] = "(kcal/mol A)^(-1/2)";

// Main SI units Converters
const double /* 1 */ METRE = 1e10; // Angstrom
                                   //  1 J = 1 kg m^2 s^-2
                                   //  e_charge C V = e_charge J = 1 eV
const double /* 1 */ KG =
    Na_avogadro * 1e7 / calory;                 // kcal/mol Angstrom^-2 fs^2
const double /* 1 */ SECOND = 1e15;             // fs
const double /* 1 */ AMPERE = 1e-15 / e_charge; // e/fs
#endif
const char *const TEMPERATURE = ENERGY;
// SI Units
const double METRE_PER_SECOND = METRE / SECOND;
const double METRE_PER_SECOND2 = METRE_PER_SECOND / SECOND;
const double /* 1 */ NEWTON = KG * METRE_PER_SECOND2;
const double /* 1 */ JOULE = NEWTON * METRE;
const double /* 1 */ COULOMB = AMPERE * SECOND;
const double WATT = JOULE / SECOND;
const double VOLT = WATT / AMPERE;
const double /* 1 */ KELVIN = k_boltzmann * JOULE;
// SI Prefix
const double /* 1 */ KILO = 1e3;
const double NANO = 1e-9;
const double FEMTO = 1e-15;
// Derived SI units
const double NM = NANO * METRE;
const double GRAM = KG / KILO;
const double /* 1 */ GRAM_PER_MOL = GRAM / MOL;
const double /* 1 */ G_PER_MOL = GRAM_PER_MOL;
const double /* 1 */ KJ_PER_MOL = KILO * JOULE / MOL;
const double /* 1 */ KJ = KILO * JOULE;
const double FS = FEMTO * SECOND;
// Common non SI Units
const double /* 1 */ ANGSTROM = 1e-10 * METRE;
const double /* 1 */ ANGSTROM2 = ANGSTROM * ANGSTROM;
const double /* 1 */ ANGSTROM3 = ANGSTROM * ANGSTROM * ANGSTROM;
const double ECHARGE = e_charge * COULOMB;
const double EV = ECHARGE * VOLT;
// Deprecated units
const double /* 1 */ CALORY = calory * JOULE;
const double /* 1 */ KCAL_PER_MOL = KILO * CALORY / MOL;
// Atomic Units
const double BOHR = 5.291772108e-11 * METRE;
const double HARTREE = 4.35974417e-18 * JOULE;
// CGS units
const double ERGS = 1e-7 * JOULE;
// Constants
const double ONE_OVER_4_PI_EPSILON0 =
    one_over_4_pi_epsilon0 * METRE * JOULE / COULOMB / COULOMB;
} // namespace unit_system
} // namespace forcefields

#endif
