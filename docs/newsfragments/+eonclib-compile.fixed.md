eonclib compiles after the Parameters accessor split: pot TUs include
Parameters.h, write-holes use ParametersLoadAccess, and unit tests alias
eonc types (including EigenmodeStrategy) from the shared test header.
