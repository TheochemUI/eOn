# Single place to bump action / runner / tool pins for generated release workflows.
# Workflow modules import via steps.ncl (`S.pins`) — do not hardcode action refs in
# release.ncl / release_prepare.ncl / towncrier_check.ncl.
# Hash pins (zizmor unpinned-uses). Comments name the tag the hash peeled from.
{
  # --- GitHub Actions ---
  # actions/checkout@v4.2.2
  checkout = "actions/checkout@11bd71901bbe5b1630ceea73d27597364c9af683",
  # actions/setup-python@v5.6.0
  setup_python = "actions/setup-python@a26af69be951a213d495a4c3e4e4022e16d87065",
  # actions/upload-artifact@v4.6.2
  upload_artifact = "actions/upload-artifact@ea165f8d65b6e75b540449e92b4886f43607fa02",
  # actions/download-artifact@v4.3.0
  download_artifact = "actions/download-artifact@d3f86a106a0bac45b974a628896c90dbdf5c8093",
  # pypa/gh-action-pypi-publish@v1.9.0
  pypi_publish = "pypa/gh-action-pypi-publish@ec4db0b4ddc65acdf4bff5fa45ac92d78b56bdf0",

  # --- defaults for generated jobs ---
  python_version = "3.12",
  ubuntu = "ubuntu-latest",
}
