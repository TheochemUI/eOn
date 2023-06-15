#include "Matter.h"
#include "Parameters.h"
#include <cstdlib>

int main(void) {
  string confile("pos.con");
  auto parameters = std::make_shared<Parameters>();
  auto pot = helper_functions::makePotential(parameters);
  Matter *matter = new Matter(pot, parameters);
  matter->con2matter(confile);
  matter->writeTibble("rSysdat.txt"s);
  delete matter;
  return EXIT_SUCCESS;
}
