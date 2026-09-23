On Apple x86_64, the SIGFPE handler now sets MXCSR exception masks in the saved context, and disableFPE writes those masks back instead of restoring the unmasked environment from enableFPE.
