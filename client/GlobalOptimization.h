#ifndef GLOBALOPTIMIZATION_H
#define GLOBALOPTIMIZATION_H
#include "Parameters.h"

class GlobalOptimization {
public:
  GlobalOptimization(Parameters *params);
  ~GlobalOptimization(void);
  void run(void);

private:
  Parameters *parameters;
};

#endif
