# PyPI distribution config for eOn (TheochemUI/eOn).
#
# IMPORTANT: pypi.org/project/eon/ is already taken by an unrelated package
# ("EoN" — Epidemics on Networks, springer-math). Do NOT publish as `eon`.
# Importable top-level package in this repo remains `eon/` (Python layout);
# the *distribution* / pip install name is below.
#
# First-time maintainer: create empty project on pypi.org (or TestPyPI) with
# this exact name, then attach Trusted Publishing:
#   owner/repo: TheochemUI/eOn (or your fork), workflow: release.yml,
#   environment: release, project: <distribution_name>
# Or set repository secret PYPI_API_TOKEN scoped to that project.
#
# Wire-up: imported by release.ncl; keep pyproject.toml [project].name in sync.
# Regenerate workflows after edits: ./ci/gha/gen.sh
{
  # Distribution / PyPI project name (pip install <this>)
  distribution_name = "eon-akmc",

  # Human label in logs / docs
  product_label = "eOn (kinetic Monte Carlo / AKMC server+client)",

  # GitHub Environment name used by release.yml pypi job
  gha_environment = "release",

  # Project page URLs (must match the distribution_name path segment)
  pypi_project_url = "https://pypi.org/p/eon-akmc",
  pypi_simple_check = "https://pypi.org/simple/eon-akmc/",
  testpypi_project_url = "https://test.pypi.org/p/eon-akmc",

  # Trusted publisher / token notes (rendered in docs, not in YAML)
  trusted_publisher_workflow = "release.yml",
  trusted_publisher_environment = "release",
  token_secret_name = "PYPI_API_TOKEN",

  # Occupied name we must not target
  blocked_pypi_name = "eon",
  blocked_pypi_reason = "Epidemics on Networks (EoN) — https://pypi.org/project/eon/",
}
