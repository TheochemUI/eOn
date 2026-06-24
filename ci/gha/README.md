# CI workflow generation (Nickel)

Release-related GitHub Actions workflows are **generated** from Nickel sources
(rgpot `ci/gha/workflow.ncl` pattern; see internal note
`Software/rgpot/release_process_v1.2.0_lockstep_potctl_cog_towncrier.org` § Nickel
build matrix — eOn extends the idea to **release** workflows because we do not
have a combinatorial build matrix in Nickel yet).

| Nickel source | Generated YAML (do not hand-edit) |
|---------------|-----------------------------------|
| [`common.ncl`](common.ncl) | Shared pins / step builders (imported) |
| [`release.ncl`](release.ncl) | [`.github/workflows/release.yml`](../../.github/workflows/release.yml) |
| [`release_prepare.ncl`](release_prepare.ncl) | [`.github/workflows/release-prepare.yml`](../../.github/workflows/release-prepare.yml) |
| [`towncrier_check.ncl`](towncrier_check.ncl) | [`.github/workflows/towncrier.yml`](../../.github/workflows/towncrier.yml) |

Other workflows under `.github/workflows/` (`ci_*.yml`) remain hand-maintained
YAML unless/until ported.

## Regenerate

Requires [Nickel](https://nickel-lang.org/) ≥ 1.x (`nickel` on PATH, or pixi
`cigen` environment when configured):

```bash
# all three release-related workflows
./ci/gha/gen.sh
# or one at a time:
nickel export --format yaml ci/gha/release.ncl -o .github/workflows/release.yml
nickel export --format yaml ci/gha/release_prepare.ncl -o .github/workflows/release-prepare.yml
nickel export --format yaml ci/gha/towncrier_check.ncl -o .github/workflows/towncrier.yml
```

With pixi (after `cigen` feature is in the lockfile on your machine):

```bash
# pixi cigen env optional; use ./ci/gha/gen.sh with nickel on PATH
pixi r -e cigen nickel format ci/gha/*.ncl   # optional format
```

## Design notes (rgpot-aligned, eOn-native)

- **Modular steps**: `common.ncl` holds action pins (`checkout@v4`, `setup-python@v5`,
  artifact actions, `pypa/gh-action-pypi-publish`) and `run_script` / `setup_python`
  builders — same idea as rgpot's `run_pixi` / `checkout_full` / `ensure_potctl`.
- **Staged `release.ncl` jobs**: `gate` → `tarball` → `gh-release` → `pypi` |
  `pypi-skip-rc` (RC tags skip PyPI; stable uses environment `release` + OIDC/token).
- **Gate uses** `scripts/release_assert.py` (lockstep + CHANGELOG), not rgpot `potctl`.
- **`release_prepare.ncl`**: PR + `workflow_dispatch` dry-run (assert, pytest,
  towncrier draft, optional cog); path filters include `ci/gha/**` so Nickel edits
  re-run prepare.
- **`towncrier_check.ncl`**: `towncrier check` on shipped-path PRs.
- **Contract**: edit `.ncl` → run `gen.sh` → commit **both** source and generated
  YAML in the same change (like rgpot `chore/modernize-nickel-gha`). Drift between
  only-YAML edits and Nickel source is a bug.

## Maintainer docs

Full cut checklist (cog, PyPI, feedstock, incomplete releases, Doxygen deferral):
[`docs/source/devdocs/release.md`](../../docs/source/devdocs/release.md).

## PyPI (`pypi.ncl`)

[`pypi.ncl`](pypi.ncl) holds the **distribution** name (`eon-akmc`), blocked name
(`eon` — epidemics package on PyPI), and trusted-publisher / environment fields.
Imported by `release.ncl` for the `pypi` job `environment.url`. Keep
`pyproject.toml` `[project].name` equal to `distribution_name`.
