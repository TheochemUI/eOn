"""Validated constructors for chemist-facing algorithm objects.

When ``pydantic`` is installed (``pip/uv install 'pyeonclient[models]'``),
:class:`~pyeonclient.models.DimerSpec` / :class:`~pyeonclient.models.NebSpec`
provide typed validation. Without pydantic, a small pure-Python check still
enforces ``accelerant="gp"`` ⇒ ``method="improved"`` (same rule as C++).
"""

from __future__ import annotations

from typing import Any, Optional, Union


def _light_dimer_kwargs(
    method: str = "improved",
    accelerant: str | None = None,
) -> dict[str, str]:
    m = (method or "improved").strip().lower()
    if m in ("", "dimer"):
        m = "improved"
    acc = (accelerant or "").strip().lower()
    if acc in ("", "none"):
        acc = ""
    elif acc == "gp":
        if m != "improved":
            raise ValueError(
                'accelerant="gp" is only valid with method="improved" '
                f'(got method="{method}"). '
                "GP-dimer is improved + GP (AtomicGPDimer), not "
                "classic/lanczos/davidson. "
                "For typed errors install: pip/uv install 'pyeonclient[models]'"
            )
    else:
        raise ValueError(
            f'accelerant must be "" / "none" / "gp", got {accelerant!r}'
        )
    if m not in ("improved", "classic", "lanczos", "davidson"):
        raise ValueError(
            f'method must be improved|classic|lanczos|davidson, got {method!r}'
        )
    return {"method": m, "accelerant": acc}


def _light_neb_accelerant(accelerant: str | None = None) -> str:
    acc = (accelerant or "").strip().lower()
    if acc in ("", "none"):
        return ""
    if acc == "gp":
        return "gp"
    raise ValueError(
        f'NEB accelerant must be "" / "none" / "gp", got {accelerant!r}'
    )


def Dimer(
    matter: Any,
    parameters: Any,
    potential: Any,
    method: str = "improved",
    accelerant: str | None = None,
    *,
    spec: Any = None,
) -> Any:
    """Min-mode entry. Default ``method="improved"``.

    ``accelerant="gp"`` only with improved. Prefer
    ``spec=DimerSpec(...)`` when ``pyeonclient[models]`` is installed.
    """
    from pyeonclient import _core

    if spec is not None:
        try:
            from pyeonclient.models import DimerSpec
        except ImportError as e:
            raise ImportError(
                "Dimer(spec=...) requires pydantic. "
                "Install with: pip/uv install 'pyeonclient[models]'"
            ) from e
        if not isinstance(spec, DimerSpec):
            raise TypeError("spec must be a DimerSpec")
        kw = spec.core_kwargs()
    else:
        try:
            from pyeonclient.models import parse_dimer_spec

            kw = parse_dimer_spec(method=method, accelerant=accelerant).core_kwargs()
        except ImportError:
            kw = _light_dimer_kwargs(method=method, accelerant=accelerant)
    return _core.Dimer(matter, parameters, potential, **kw)


def NEB(
    *args: Any,
    accelerant: str | None = None,
    path_init: str = "idpp",
    n_images: int = 10,
    spec: Any = None,
    parameters: Any = None,
    potential: Any = None,
    **kwargs: Any,
) -> Any:
    """NEB entry. ``accelerant="gp"`` selects GP-surrogate band when built."""
    from pyeonclient import _core

    if spec is not None:
        try:
            from pyeonclient.models import NebSpec
        except ImportError as e:
            raise ImportError(
                "NEB(spec=...) requires pydantic. "
                "Install with: pip/uv install 'pyeonclient[models]'"
            ) from e
        if not isinstance(spec, NebSpec):
            raise TypeError("spec must be a NebSpec")
        acc = spec.core_accelerant()
        apply_params = spec.apply_to_parameters
    else:
        try:
            from pyeonclient.models import parse_neb_spec

            s = parse_neb_spec(
                accelerant=accelerant,
                path_init=path_init,
                n_images=n_images,
            )
            acc = s.core_accelerant()
            apply_params = s.apply_to_parameters
        except ImportError:
            acc = _light_neb_accelerant(accelerant)
            apply_params = None

    # parameters / potential are named above, so they never reach **kwargs;
    # dispatch on what is bound rather than on the kwargs dict.
    initial = kwargs.pop("initial", None)
    final = kwargs.pop("final", None)
    path = kwargs.pop("path", None)
    if len(args) == 4:
        initial, final, parameters, potential = args
    elif len(args) == 3:
        path, parameters, potential = args
    elif len(args) == 2:
        initial, final = args
    elif len(args) == 1:
        if final is not None:
            initial = args[0]
        else:
            path = args[0]
    elif args:
        raise TypeError(f"NEB takes at most 4 positional arguments, got {len(args)}")
    if kwargs:
        raise TypeError(f"NEB got unexpected keyword arguments: {sorted(kwargs)}")

    if initial is not None and final is not None:
        if apply_params is not None and parameters is not None:
            apply_params(parameters)
        return _core.NEB(initial, final, parameters, potential, acc)
    if path is not None:
        return _core.NEB(path, parameters, potential, acc)
    raise TypeError(
        "NEB expected (initial, final, parameters, potential) or "
        "(path, parameters, potential)"
    )


__all__ = ["Dimer", "NEB"]
