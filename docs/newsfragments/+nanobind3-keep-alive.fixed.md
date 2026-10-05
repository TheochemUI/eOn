`tie_lifetime` calls `keep_alive_py` on nanobind 3 and `nb::detail::keep_alive` on nanobind 2, so a returned path keeps the Parameters owner alive.
