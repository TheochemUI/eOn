# Advisory release preparation (does not publish).
# Regenerate: nickel export --format yaml ci/gha/release_prepare.ncl -o .github/workflows/release-prepare.yml
let C = import "./common.ncl" in

let release_paths = [
  "cog.toml",
  "towncrier.toml",
  "scripts/release_assert.py",
  "tests/release/**",
  "pyproject.toml",
  "pixi.toml",
  "docs/newsfragments/**",
  "CHANGELOG.md",
  ".github/workflows/release.yml",
  ".github/workflows/release-prepare.yml",
  ".github/workflows/towncrier.yml",
  "ci/gha/**",
]
in

let job_prepare = {
  "runs-on" = "ubuntu-latest",
  steps = [
    C.checkout_full,
    C.setup_python C.python_version,
    C.run_script "Assert lockstep surfaces agree" m%"
      python3 scripts/release_assert.py
      python3 -c "import tomllib, pathlib; tomllib.loads(pathlib.Path('cog.toml').read_text()); print('cog.toml parse ok')"
      python3 -c "import tomllib, pathlib; tomllib.loads(pathlib.Path('pyproject.toml').read_text()); print('pyproject.toml parse ok')"
    "%,
    C.run_script "Unit tests for release_assert" m%"
      python3 -m pip install --quiet pytest
      python3 -m pytest tests/release/ -q
    "%,
    C.run_script "Towncrier draft (fragments must build cleanly)" m%"
      python3 -m pip install --quiet towncrier
      if compgen -G "docs/newsfragments/*.md" > /dev/null 2>&1 || \
         compgen -G "docs/newsfragments/+*.md" > /dev/null 2>&1; then
        # towncrier.toml carries name/version so draft does not import eon
        towncrier build --draft
      else
        echo "No newsfragments; skip draft (ok if this PR only touches tooling)"
      fi
    "%,
    C.run_script "Optional cog dry-run" m%"
      if command -v cog >/dev/null 2>&1; then
        VER=$(python3 scripts/release_assert.py --print-version)
        cog bump --version "$VER" --dry-run 2>&1 | head -40 || \
          echo "note: cog dry-run non-zero (often dirty/unpushed PR tip); lockstep ok above"
      else
        echo "cog not installed on runner; skip cog dry-run (run locally with cocogitto 7+)"
      fi
    "%,
    {
      name = "Cog bump dry-run (workflow_dispatch)",
      "if" = "github.event_name == 'workflow_dispatch'",
      run = m%"
        if ! command -v cog >/dev/null 2>&1; then
          echo "Install cocogitto 7+ locally for cog dry-run; runner has no cog."
          exit 0
        fi
        BUMP="${{ inputs.bump }}"
        case "$BUMP" in
          auto)  cog bump --auto --dry-run ;;
          patch) cog bump --patch --dry-run ;;
          minor) cog bump --minor --dry-run ;;
          major) cog bump --major --dry-run ;;
        esac
      "%,
    },
  ],
}
in

{
  name = "Release prepare (dry-run)",
  on = {
    workflow_dispatch = {
      inputs = {
        bump = {
          description = "SemVer component for cog dry-run (if cog installed on runner)",
          required = true,
          default = "auto",
          type = "choice",
          options = ["auto", "patch", "minor", "major"],
        },
      },
    },
    pull_request = { paths = release_paths },
  },
  concurrency = {
    group = "release-prepare-${{ github.ref }}",
    "cancel-in-progress" = true,
  },
  permissions = { contents = "read" },
  jobs = { prepare = job_prepare },
}
