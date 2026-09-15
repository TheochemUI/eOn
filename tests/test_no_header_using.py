"""Job headers in this slice must not inject using eonc:: into the global ns."""

from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
# Optimizer.h still has the factory alias.
JOB_HEADERS = sorted(
    (ROOT / "include" / "eon").glob("*Job.h")
) + [
    ROOT / "include" / "eon" / name
    for name in (
        "ConjugateGradients.h",
        "FIRE.h",
        "LBFGS.h",
        "Quickmin.h",
        "SteepestDescent.h",
        "Dimer.h",
        "ImprovedDimer.h",
        "Lanczos.h",
        "Davidson.h",
        "LORRotation.h",
        "LowestEigenmode.h",
        "MinModeSaddleSearch.h",
        "EigenmodeStrategy.h",
        "GleThermostat.h",
        "BondBoost.h",
        "Hessian.h",
        "MonteCarlo.h",
        "GlobalOptimization.h",
        "IRACompare.h",
        "IDPPObjectiveFunction.hpp",
        "ObjectiveFunction.h",
        "SaddleSearchMethod.h",
        "SurrogatePotential.h",
        "PotRegistry.h",
        "BasinHoppingSaddleSearch.h",
        "BiasedGradientSquaredDescent.h",
        "DynamicsSaddleSearch.h",
        "ARTnSaddleSearch.h",
        "AtomicGPDimer.h",
        "Dynamics.h",
        "ServeMode.h",
        "NudgedElasticBand.h",
        "Optimizer.h",
        "ServeRpcServer.h",
        "Potential.h",
    )
]


def test_job_headers_have_no_file_scope_using():
    offenders: list[str] = []
    for path in JOB_HEADERS:
        for i, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1):
            stripped = line.strip()
            if stripped.startswith("using eonc::") or stripped.startswith(
                "using namespace "
            ):
                offenders.append(f"{path.name}:{i}:{stripped}")
    assert not offenders, "file-scope using remains:\n" + "\n".join(offenders)
