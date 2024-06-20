#ifndef FPE_HANDLER_H_
#define FPE_HANDLER_H_

// Floating Point Trapping. It is platform specific!
// This causes the program to crash on divison by zero,
// invalid operations, and overflows.
void enableFPE(void);

#endif // FPE_HANDLER_H_
