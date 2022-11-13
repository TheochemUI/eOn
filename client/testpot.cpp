#include "Matter.h"
#include "Parameters.h"
#include <cstdlib>

int main(void){
  string confile("pos.con");
  Parameters *parameters = new Parameters;
  parameters->potential = "cuh2_pot";
  Matter *matter = new Matter(parameters);
  matter->con2matter(confile);
  std::cout<<matter->getPotentialEnergy()<<std::endl;
  delete matter;
  delete parameters;
  return EXIT_SUCCESS;
}
