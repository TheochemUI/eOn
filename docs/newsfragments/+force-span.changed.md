``Potential::force`` has a ``std::span`` overload that checks sizes
before the raw C-array virtual. Matter, ``get_ef``, and surrogate
``get_ef_var`` use it. Fortran/FFI loaders keep the pointer API.
