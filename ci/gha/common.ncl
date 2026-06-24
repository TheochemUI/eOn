# Shared GHA pins and step builders (rgpot ci/gha/workflow.ncl style).
# Used by release.ncl, release_prepare.ncl, towncrier_check.ncl.
{
  checkout_v4 = "actions/checkout@v4",
  setup_python_v5 = "actions/setup-python@v5",
  upload_artifact_v4 = "actions/upload-artifact@v4",
  download_artifact_v4 = "actions/download-artifact@v4",
  pypi_publish = "pypa/gh-action-pypi-publish@release/v1",
  python_version = "3.12",

  checkout_full = {
    uses = "actions/checkout@v4",
    with = { "fetch-depth" = 0 },
  },

  checkout_default = { uses = "actions/checkout@v4" },

  setup_python = fun ver =>
    {
      uses = "actions/setup-python@v5",
      with = { "python-version" = ver },
    },

  # Multistring runner script (avoid `name = name` field/param shadowing).
  run_script = fun step_name script_body =>
    {
      name = step_name,
      run = script_body,
    },
}
