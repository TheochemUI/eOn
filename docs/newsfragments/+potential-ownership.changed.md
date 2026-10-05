makePotential returns exclusive ownership. Point and minimization keep that
unique_ptr until one Matter takes it. Shared ownership is explicit via
sharePotential at multi-owner boundaries.
