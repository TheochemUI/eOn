The SIGFPE continue path now demotes ARM FPCR trap-enable bits (Linux fpsimd_context and Apple neon state). enableFPE on Apple aarch64 sets `__fpcr`, not `__fpsr`. x86 MXCSR masking is unchanged.
