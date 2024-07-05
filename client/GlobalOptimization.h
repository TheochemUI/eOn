#pragma once
#include "Parameters.h"

class GlobalOptimization {
public:
  GlobalOptimization(Parameters *params);
  ~GlobalOptimization(void);
  void run(void);

private:
  Parameters *parameters;
};
