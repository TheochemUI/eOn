"""NEB API composition (L2)."""

from __future__ import annotations

from typing import Any, Optional, Union

from eon_schema._deps import ensure_import
from eon_schema.fields.enums import Accelerant, PathInit

_NebSpecCls: Any = None


def _neb_spec_cls() -> Any:
    global _NebSpecCls
    if _NebSpecCls is not None:
        return _NebSpecCls
    pd = ensure_import(
        "pydantic",
        install_spec="pydantic>=2",
        extra_hint="Or: pip install 'eon-schema[pydantic]'",
    )
    BaseModel = pd.BaseModel
    ConfigDict = pd.ConfigDict
    Field = pd.Field

    class NebSpec(BaseModel):
        """Chemist-facing NEB request: path init, CI, energy weighting, OCI-MMF.

        Call :meth:`apply_to_parameters` (or pass ``spec=`` into
        :func:`pyeonclient.NEB`) so one model owns the knobs instead of
        hand-assigning ``Parameters.neb_*`` fields.
        """

        model_config = ConfigDict(extra="forbid", use_enum_values=False)

        accelerant: Accelerant = Field(default=Accelerant.none)
        path_init: PathInit = Field(default=PathInit.idpp)
        n_images: int = Field(default=10, ge=1)
        # FILE init: list file (one .con path per line), e.g. idppPath.dat
        path_list: Optional[str] = Field(default=None)
        minimize_endpoints: bool = False
        # Climbing image
        climbing_image: bool = True
        climbing_converged_only: bool = True
        climbing_band_slack: float = Field(default=10.0, gt=0.0)
        ci_after: float = Field(default=0.5, ge=0.0)
        ci_after_rel: float = Field(default=0.8, ge=0.0)
        # Energy-weighted springs
        energy_weighted: bool = False
        ew_ksp_min: float = 0.972
        ew_ksp_max: float = 9.72
        # Off-path CI / hybrid MMF (OCI-NEB)
        ci_mmf: bool = False
        ci_mmf_after: float = Field(default=0.1, ge=0.0)
        ci_mmf_after_rel: float = Field(default=0.5, ge=0.0)
        ci_mmf_angle: float = Field(default=0.9, ge=0.0, le=1.0)
        ci_mmf_nsteps: int = Field(default=1000, ge=1)
        # Zoom-NEB image packing around the climbing image. Composes with ci_mmf.
        zoom: bool = False
        zoom_alpha: float = Field(default=0.5, gt=0.0, lt=1.0)
        zoom_offset: int = Field(default=1, ge=1)
        zoom_mode: str = "auto"
        zoom_after: float = Field(default=0.0, ge=0.0)
        zoom_interpolation: str = "cubic"
        zoom_ci_stability: int = Field(default=5, ge=0)
        zoom_max_iterations: int = Field(default=0, ge=0)
        # Optimizer (mapped onto neb + opt Parameters fields)
        max_iterations: int = Field(default=1000, ge=1)
        force_tolerance: float = Field(default=0.01, gt=0.0)
        max_move: float = Field(default=0.1, gt=0.0)
        write_movies: bool = False
        random_seed: Optional[int] = None

        def core_accelerant(self) -> str:
            return (
                ""
                if self.accelerant is Accelerant.none
                else self.accelerant.value
            )

        def apply_to_parameters(self, params: Any) -> Any:
            try:
                from pyeonclient._core import NEBInit  # type: ignore
            except ImportError as e:
                raise ImportError(
                    "NebSpec.apply_to_parameters requires pyeonclient"
                ) from e
            mapping = {
                PathInit.linear: NEBInit.LINEAR,
                PathInit.idpp: NEBInit.IDPP,
                PathInit.idpp_collective: NEBInit.IDPP_COLLECTIVE,
                PathInit.sidpp: NEBInit.SIDPP,
                PathInit.sidpp_zbl: NEBInit.SIDPP_ZBL,
                PathInit.file: NEBInit.FILE,
            }
            params.neb_images = int(self.n_images)
            if self.path_list:
                params.neb_init_method = NEBInit.FILE
                params.neb_initial_path = str(self.path_list)
            else:
                params.neb_init_method = mapping[self.path_init]
                params.neb_initial_path = ""
            params.neb_minimize_endpoints = bool(self.minimize_endpoints)
            params.neb_climbing_image = bool(self.climbing_image)
            params.neb_climbing_converged_only = bool(self.climbing_converged_only)
            params.neb_climbing_band_slack = float(self.climbing_band_slack)
            params.neb_ci_after = float(self.ci_after)
            params.neb_ci_after_rel = float(self.ci_after_rel)
            params.neb_energy_weighted = bool(self.energy_weighted)
            params.neb_ew_ksp_min = float(self.ew_ksp_min)
            params.neb_ew_ksp_max = float(self.ew_ksp_max)
            params.neb_ci_mmf = bool(self.ci_mmf)
            params.neb_ci_mmf_after = float(self.ci_mmf_after)
            params.neb_ci_mmf_after_rel = float(self.ci_mmf_after_rel)
            params.neb_ci_mmf_angle = float(self.ci_mmf_angle)
            params.neb_ci_mmf_nsteps = int(self.ci_mmf_nsteps)
            params.neb_zoom = bool(self.zoom)
            params.neb_zoom_alpha = float(self.zoom_alpha)
            params.neb_zoom_offset = int(self.zoom_offset)
            params.neb_zoom_mode = str(self.zoom_mode)
            params.neb_zoom_after = float(self.zoom_after)
            params.neb_zoom_interpolation = str(self.zoom_interpolation)
            params.neb_zoom_ci_stability = int(self.zoom_ci_stability)
            params.neb_zoom_max_iterations = int(self.zoom_max_iterations)
            params.neb_max_iterations = int(self.max_iterations)
            params.neb_force_tolerance = float(self.force_tolerance)
            params.opt_max_iterations = int(self.max_iterations)
            params.opt_converged_force = float(self.force_tolerance)
            params.opt_max_move = float(self.max_move)
            params.write_movies = bool(self.write_movies)
            if self.random_seed is not None:
                params.random_seed = int(self.random_seed)
            return params

    _NebSpecCls = NebSpec
    return _NebSpecCls


def __getattr__(name: str) -> Any:
    if name == "NebSpec":
        return _neb_spec_cls()
    raise AttributeError(name)


def parse_neb_spec(
    accelerant: Union[str, Accelerant, None] = None,
    *,
    path_init: Union[str, PathInit] = "idpp",
    n_images: int = 10,
    path_list: Optional[str] = None,
    energy_weighted: bool = False,
    ci_mmf: bool = False,
    max_iterations: int = 1000,
    force_tolerance: float = 0.01,
    spec: Any = None,
    **extra: Any,
) -> Any:
    if spec is not None:
        return spec
    acc: Any = Accelerant.none if accelerant in (None, "") else accelerant
    payload: dict[str, Any] = {
        "accelerant": acc,
        "path_init": path_init,
        "n_images": n_images,
        "path_list": path_list,
        "energy_weighted": energy_weighted,
        "ci_mmf": ci_mmf,
        "max_iterations": max_iterations,
        "force_tolerance": force_tolerance,
    }
    payload.update(extra)
    return _neb_spec_cls().model_validate(payload)
