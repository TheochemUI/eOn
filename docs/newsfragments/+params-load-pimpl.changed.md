``Parameters`` stores load state and option groups in a private ``Impl``.
``sizeof(Parameters)`` is that pointer. Const accessors and
``ParametersLoadAccess`` are the read and write surface. Option-group
types stay in ``ParametersOptions.h``. ``Matter`` still exposes Eigen.
