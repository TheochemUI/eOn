Suppress FPE trapping during libtorch operations in `MetatomicPotential`
constructor and `force()`, preventing SIGFPE from benign NaN/Inf produced by
SiLU (sleef) and autograd internals. Follows the existing `FPEHandler` pattern
from ASE_ORCA, ASE_NWCHEM, and AtomicGPDimer.
