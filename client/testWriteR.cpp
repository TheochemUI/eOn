#include "Matter.h"
#include "Parameters.h"
#include <cstdlib>

int main(void){
  string confile("pos.con");
  Parameters *parameters = new Parameters;
  Matter *matter = new Matter(parameters);
  matter->con2matter(confile);
  matter-> writeTibble("rSysdat.txt"s);
  delete matter;
  delete parameters;
  return EXIT_SUCCESS;
}
