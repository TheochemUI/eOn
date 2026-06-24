# towncrier check when shipped surfaces change (rgpot towncrier.yml intent).
# Regenerate: nickel export --format yaml ci/gha/towncrier_check.ncl -o .github/workflows/towncrier.yml
let C = import "./common.ncl" in

let shipped_paths = [
  "eon/**",
  "client/**",
  "meson.build",
  "meson_options.txt",
  "CMakeLists.txt",
  "pyproject.toml",
  "pixi.toml",
  "scripts/**",
  "docs/newsfragments/**",
  "towncrier.toml",
  "CHANGELOG.md",
]
in

{
  name = "Towncrier check",
  on = { pull_request = { paths = shipped_paths } },
  permissions = {
    contents = "read",
    "pull-requests" = "read",
  },
  jobs = {
    towncrier = {
      "runs-on" = "ubuntu-latest",
      steps = [
        C.checkout_full,
        C.setup_python C.python_version,
        C.run_script "Install towncrier" "python3 -m pip install --quiet towncrier",
        {
          name = "towncrier check",
          env = { BASE_REF = "${{ github.event.pull_request.base.ref }}" },
          run = m%"
            git fetch origin "$BASE_REF"
            if towncrier check --compare-with "origin/${BASE_REF}"; then
              echo "towncrier check ok"
            else
              echo "::error::Missing or invalid newsfragment for shipped-path changes."
              echo "Add docs/newsfragments/<issue>.<type>.md (see docs/newsfragments/README.md)"
              echo "Types: security removed deprecated added dev changed fixed"
              exit 1
            fi
          "%,
        },
      ],
    },
  },
}
