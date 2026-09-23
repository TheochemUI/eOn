AtomicGPDimer restores floating-point traps if execute throws, so a saddle search that catches the failure does not leave later force calls running with traps masked.
