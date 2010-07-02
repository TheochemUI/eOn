#ifndef SYTEM_UNIT_H

//Jean-Claude C. Berthet, University of Iceland, 2006

/** @file
This file contains the definition of some common constants. These constants are usually expressed
in the Systeme International (SI). It also contains local unit names (e.g. E_UNIT) and unit converters.

COHERENTE LOCAL UNIT SYSTEM AND CONVETERS:
Two constraints:
First, minimizing the number of operations during heavy computation, requires the use a coherent system of unit.
Second, changing to another system must be easy.
The unit system is based on for starting dimensions: distance, mass, time. Therefore, once the 
units for those three dimensions are chosen, the units for all other quantities are set. For 
instance, if we choose the units meter, kg/mol and second, the energy will be in 
J/mol(=kg/mol m^2 s^-2), the pressure in Pa/mol (=kg/mol m^-1 s^-2), etc... Different quantities of same 
dimension will always be expressed with the same unit. The temperature will always be stored as an 
energy (J/mol in the example given before).

UNIT NAMES (E_UNIT, etc...)
The constants named after a quantity (e.g. ENERGY) contain the unit used by the program for this quantity. It is recommended to use these 
constants every time the program needs to display a quantity in the units used by the program. Then, if the system 
of unit is changed only the constant need to be updated. 
Example: the local unit system uses kJ/mol as the energy unit. We wish to print an energy returned by a 
function 'E()'. Then the expression should look like:
cout << E() << ' ' << E_UNIT << endl;
If the program is later modified to use a different energy unit, eV for instance, 
we simply need to change the value of E_UNIT by "eV".

However, this is not the best way to display results as we will certainly end up displaying quantities with 
uncommon units. For example, the temperature will always be displayed in the energy unit, most of the time kJ/mol 
or eV. I have then defined what I shall call 'converters'.

CONVERTERS:
A converter is a double constant which has the name of a unit (e.g. PASCAL, KELVIN, ANGSTROM...). The constant contains 
the value of one quantity of the unit it names expressed in the system of unit used by the program.
Example with ANGSTROM: We assume that the program uses the meter as unit of distance. We wish to display the 
content of a variable 'x' which is a distance in meter, but we want to display it in Angstrom. To make the 
conversion we simply need to divide the quantity in meter by one Angstrom:
(x meter)/(1 angstrom)=(x meter)/(10^-10 meter)=x*10^-10. So, the expression to display will be like:
cout << x/ANGSTROM << " Angstrom" << endl;
As well to input a value in angstrom the code will be similar to:
double x_ang;
cin >> x_ang;
x= x_ang*ANGSTROM;
Of course, if the unit system used by the program already uses the angstrom, the constant ANGSTROM will be equal to 
one. But it is still better to use the converter as it makes the program easier to update if things were to change.
*/

/** Units based on eV, Angstrom, and femtosecond. If the SYTEM_UNIT_H is defined as EV_A_FS, the system currently used is a coherent
system based on eV (electron Volt), Angstrom, and femtosecond. Then the energies returned are in eV, the forces in eV/A, and the distance 
in Angstrom.*/
#define EV_A_FS 0
#define SYTEM_UNIT_H EV_A_FS
namespace system_unit {
      // Some physical constants expressed in SI unit
      const double				R=8.31441;					//	J/mol/K
      const double				k=1.380662e-23;				//	J/K
      const double				Na=6.022045e23;				//	/mol
      const double      /* 1 */                        calory=4.1868;				//	J
      const double				e=1.6021892e-19;				//	C
      const double				pi=3.14159265358979312;
      const double      /* 1 */                        debye=3.33564e-30;			// C m
      /// This is the quantity 1/(4*pi*epsilon0) in SI units.
      const double				coulomb=8.98755178495272255e+09;	//SI
      /// This is the quantity 1/sqrt(4*pi*epsilon0) in SI units.
      const double				sqrtcoulomb=9.48026992492973368e+04;	//SI

      // System independent units
      const double			MOL=Na;
      const double			DEGREE=pi/180;//rad
      #if SYTEM_UNIT_H==EV_A_FS
            const char			ENERGY[]="eV";          //electron volt
            const char			LENGTH[]="A";           //Angstrom
            const char			TIME[]="s";                  //second

            const char			FORCE[]="eV/A";
            const char			MASS[]="eV A^-2 fs^2";
            const char			ACCELERATION[]="A fs^-2";
            const char			VELOCITY[]="A fs^-1";
            const char *const                 TEMPERATURE=ENERGY;
            const char			PRESSURE[]="eV A^-3";
            const char                              CHARGE[]="(eV A)^(-1/2)";

            //Converters
            const double	/* 1 */	EV=1;                                                             // eV
            const double	/* 1 */	KJ_PER_MOL=1000/(e*Na);		// eV	
            const double	/* 1 */	KELVIN=k/e;                                                 // eV

            const double	/* 1 */	ANGSTROM=1;                                             // A
            const double /* 1 */            METER=1e10;

            const double /* 1 */	FS=1;                                                             //s
            const double /* 1 */	SECOND=1e15;                                            //  s
      #endif
      const double	/* 1 */           KCAL_PER_MOL=KJ_PER_MOL*calory;
      const double	/* 1 */           ANGSTROM3=ANGSTROM*ANGSTROM*ANGSTROM;
      const double /* 1 */                 KG_PER_MOL=1e-3*KJ_PER_MOL*SECOND*SECOND/METER/METER;
      const double /* 1 */                 G_PER_MOL=1e-3*KG_PER_MOL;                        
}
#undef EV_A_FS

#endif
