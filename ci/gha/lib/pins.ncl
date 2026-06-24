# Single place to bump action / runner / tool pins for generated release workflows.
# Workflow modules import via steps.ncl (`S.pins`) — do not hardcode action refs in
# release.ncl / release_prepare.ncl / towncrier_check.ncl.
# Hand-maintained ci_*.yml may lag on @v4; dual pins are intentional until those migrate.
{
  # --- GitHub Actions ---
  checkout = "actions/checkout@v4",
  setup_python = "actions/setup-python@v5",
  upload_artifact = "actions/upload-artifact@v4",
  download_artifact = "actions/download-artifact@v4",
  pypi_publish = "pypa/gh-action-pypi-publish@release/v1",

  # --- defaults for generated jobs ---
  python_version = "3.12",
  ubuntu = "ubuntu-latest",
}
