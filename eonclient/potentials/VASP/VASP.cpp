#include "VASP.h"

VASP::VASP(void){
    return;
}

void VASP::cleanMemory(void){
    return;
}

// pointer to number of atoms, pointer to array of positions	
// pointer to array of forces, pointer to internal energy
// adress to supercell size
void VASP::force(long N, const double *R, const long *atomicNrs, double *F, double *U, const double *box){

    writePositionsToFile(N, R, atomicNrs, box);
    system("vasp >> vaspOutput");
    readForcesFromFile(N, F, U);
    
    return;
}

void VASP::writePositionsToFile(long N, const double *R, long const *atomicNrs, const double *box){
    // Positions are scaled 
    int i=0, i_old=0;
    FILE *file;
    file = fopen("POSCAR","w");

    // header line (treated as a comment)
    i_old=0;
    fprintf(file, "%li ",atomicNrs[0]);
    for(i=0; i<N; i++){
        if(atomicNrs[i]!=atomicNrs[i_old]){
            fprintf(file, "%li ",atomicNrs[i]);
            i_old = i;
        }
    }
    fprintf(file, ": Atomic numbers\n");
    
    // boundary box
    fprintf(file, "1.0\n");
    fprintf(file, " %.8lf\t%.8lf\t%.8lf\n", box[0]/system_unit::ANGSTROM, 0.0, 0.0);
    fprintf(file, " %.8lf\t%.8lf\t%.8lf\n", 0.0, box[1]/system_unit::ANGSTROM, 0.0);
    fprintf(file, " %.8lf\t%.8lf\t%.8lf\n", 0.0, 0.0, box[2]/system_unit::ANGSTROM);

    // the number of atoms of each of the the different atomic types
    i_old=0;
    for(i=0; i<N; i++){
        if(atomicNrs[i]!=atomicNrs[i_old]){
            fprintf(file, "%li ",i-i_old);
            i_old = i;
        }
    }
    fprintf(file, "%li\n", N-i_old);

    // coordinates for all atoms
    fprintf(file, "Cartesian\n");
    for(i=0; i<N; i++){
        fprintf(file, "%.19lf\t%.19lf\t%.19lf\t F F F\n", 
                       R[0+3*i]/system_unit::ANGSTROM,
                       R[1+3*i]/system_unit::ANGSTROM,
                       R[2+3*i]/system_unit::ANGSTROM);
    }
    fclose(file);
    return;
}

void VASP::readForcesFromFile(long N, double *F, double *U){
    int i=0, line=0, length=500;
    double garbageDouble;
    char garbageChar[length], word1[length], word2[length];
    char lineAll[length];// Temporary string of character to read from the file.
    FILE *file;
    file = fopen("OUTCAR","r");
    
    while (!feof(file)){
        fgets(lineAll, length, file);
        std::sscanf(lineAll, "%s %s", word1, word2);
        
        if((!strcmp(word1, "FREE")) && (!strcmp(word2, "ENERGIE"))){
            // remove three lines and read the fourth
            for(i=0; i<4; i++)
                fgets(lineAll, length, file);
            sscanf(lineAll,"%s %s %s %lf %s %s %lf", 
                           garbageChar, garbageChar, garbageChar,
                           U,
                           garbageChar, garbageChar, &garbageDouble);
        }
        
        if((!strcmp(word1, "POSITION")) && (!strcmp(word2, "TOTAL-FORCE"))){
            // remove one line
            fgets(lineAll, length, file);
            for(i=0; i<N; i++) {
                fgets(lineAll, length, file);
                sscanf(lineAll,"%lf %lf %lf %lf %lf %lf\n", 
                               &garbageDouble, &garbageDouble, &garbageDouble,
                               &F[0+i*3], &F[1+i*3], &F[2+i*3]);
            }
        }
    };
    fclose(file);
    return;
}
