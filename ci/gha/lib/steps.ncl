# Reusable GHA step builders (rgpot ci/gha/lib/steps.ncl style, eOn-native).
# Import from workflow modules: let S = import "./lib/steps.ncl" in …
let P = import "pins.ncl" in
{
  pins = P,

  checkout_minimal = { uses = P.checkout },

  checkout_full = {
    uses = P.checkout,
    with = { "fetch-depth" = 0 },
  },

  checkout_default = { uses = P.checkout },

  setup_python = fun ver =>
    {
      uses = P.setup_python,
      with = {
        "python-version" = ver,
        # Cache pip downloads for pytest/towncrier/build installs on prepare/release.
        cache = "pip",
      },
    },

  setup_python_default = {
    uses = P.setup_python,
    with = {
      "python-version" = P.python_version,
      cache = "pip",
    },
  },

  upload_artifact = fun artifact_name artifact_path retention_days =>
    {
      name = "Upload artifact",
      uses = P.upload_artifact,
      with = {
        name = artifact_name,
        path = artifact_path,
        "retention-days" = retention_days,
      },
    },

  download_artifact = fun artifact_name dest_path =>
    {
      name = "Download artifact",
      uses = P.download_artifact,
      with = {
        name = artifact_name,
        path = dest_path,
      },
    },

  # Multistring runner script (avoid `name = name` field/param shadowing).
  run_script = fun step_name script_body =>
    {
      name = step_name,
      run = script_body,
    },

  pypi_publish_step = fun token_expr =>
    {
      name = "Publish to PyPI",
      uses = P.pypi_publish,
      with = {
        password = token_expr,
        "skip-existing" = true,
        verbose = true,
      },
    },
}
