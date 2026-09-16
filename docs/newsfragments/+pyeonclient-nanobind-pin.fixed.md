Wheel and sdist builds pin nanobind to ``>=2.2,<3``, matching pixi.
Unconstrained ``pip install -U nanobind`` pulled 3.0, which dropped
``nb::detail::keep_alive`` used by pyeonclient ``tie_lifetime``.
