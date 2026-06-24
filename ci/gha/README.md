# CI workflow generation (Nickel)

Release-related GitHub Actions workflows are **generated** from Nickel sources
(rgpot `ci/gha` tooling layer — pins/steps/drift/cigen — **not** the full
orchestrator/build matrix). eOn uses Nickel for the **release trio only**;
hand-maintained `ci_*.yml` (build/test/docs/bench) stay outside this contract.

| Nickel source | Generated YAML (do not hand-edit) |
|---------------|-----------------------------------|
| [`lib/pins.ncl`](lib/pins.ncl) | Shared action/runner pins (imported only) |
| [`lib/steps.ncl`](lib/steps.ncl) | Step builders (`checkout_*`, `setup_python_*`, `run_script`, …) |
| [`pypi.ncl`](pypi.ncl) | PyPI distribution config (`eon-akmc`; never `eon`) |
| [`release.ncl`](release.ncl) | [`.github/workflows/release.yml`](../../.github/workflows/release.yml) |
| [`release_prepare.ncl`](release_prepare.ncl) | [`.github/workflows/release-prepare.yml`](../../.github/workflows/release-prepare.yml) |
| [`towncrier_check.ncl`](towncrier_check.ncl) | [`.github/workflows/towncrier.yml`](../../.github/workflows/towncrier.yml) |
| [`common.ncl`](common.ncl) | Thin shim re-exporting `lib/steps.ncl` (compat only) |

## Regenerate (autonomous)

```bash
# preferred (hermetic nickel via pixi feature cigen)
pixi r -e cigen gen-gha

# or nickel on PATH
./ci/gha/gen.sh
```

Drift gate (fails if committed YAML ≠ export):

```bash
pixi r -e cigen gha-drift
# or
bash ci/gha/check-drift.sh
```

Optional format:

```bash
pixi r -e cigen nickel-fmt
```

**CI enforcement:** `release-prepare` runs `check-drift.sh` when `nickel` or `pixi`
is available; `ci_precommit.yml` runs the same best-effort step so hand-edited
release/prepare/towncrier YAML cannot merge unnoticed on paths that install tools.

## Design notes (rgpot-aligned, eOn-native)

- **Modular pins/steps**: workflow modules import `lib/steps.ncl` only — no raw
  `actions/checkout@…` strings in `release.ncl` / `release_prepare.ncl` /
  `towncrier_check.ncl`. Bump versions in `lib/pins.ncl` once.
- **Caching**: `setup_python_default` sets `cache: pip` so pytest/towncrier/build
  pip installs on prepare/release are faster without new workflow families.
- **Staged `release.ncl` jobs**: `gate` → `tarball` → `gh-release` → `pypi` |
  `pypi-skip-rc` (RC tags skip PyPI; stable uses environment `release` + OIDC/token).
- **Gate uses** `scripts/release_assert.py` (lockstep + CHANGELOG), not rgpot `potctl`.
- **`release_prepare.ncl`**: PR + `workflow_dispatch` dry-run (drift, assert, pytest,
  towncrier draft, optional cog); path filters include `ci/gha/**`.
- **`towncrier_check.ncl`**: `towncrier check` on shipped-path PRs.
- **Contract**: edit `.ncl` → `pixi r -e cigen gen-gha` → commit **both** source and
  generated YAML. Drift between only-YAML edits and Nickel is a bug (`gha-drift`).
- **Non-goals here**: porting `ci_build_akmc.yml` / `ci_xtb.yml` / etc. into Nickel.

Hand `ci_*.yml` may still pin `checkout@v4` while this tree centralizes pins for the
release trio only; dual pins are intentional until those workflows migrate.

## Maintainer docs

Full cut checklist (cog, PyPI, feedstock, incomplete releases, Doxygen deferral):
[`docs/source/devdocs/release.md`](../../docs/source/devdocs/release.md).

## PyPI (`pypi.ncl`)

[`pypi.ncl`](pypi.ncl) holds the **distribution** name (`eon-akmc`), blocked name
(`eon` — epidemics package on PyPI), and trusted-publisher / environment fields.
Imported by `release.ncl` for the `pypi` job `environment.url`. Keep
`pyproject.toml` `[project].name` equal to `distribution_name`.
