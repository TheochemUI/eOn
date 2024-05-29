#ifndef _EXCEPTIONSEON_H
#define _EXCEPTIONSEON_H

void enableFPE(void);
void disableFPE(void);
int getFPEState(void);
bool isFPEEnabled(void);

#endif  // _EXCEPTIONSEON_H