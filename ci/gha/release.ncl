# Tag-triggered release: gate → tarball → GitHub Release → PyPI (stable only).
# Source of truth; regenerate: nickel export --format yaml ci/gha/release.ncl -o .github/workflows/release.yml
# See docs/source/devdocs/release.md and ci/gha/README.md.
# Mirrors rgpot staged release.yml; eOn uses scripts/release_assert.py (not potctl).
let C = import "./common.ncl" in
let Pypi = import "./pypi.ncl" in

let job_gate = {
  "runs-on" = "ubuntu-latest",
  outputs = {
    is_prerelease = "${{ steps.prerelease.outputs.is_prerelease }}",
    version = "${{ steps.meta.outputs.version }}",
    tag = "${{ steps.meta.outputs.tag }}",
  },
  steps = [
    C.checkout_full,
    C.run_script "Tag / version meta" m%"
      TAG="${GITHUB_REF#refs/tags/}"
      VERSION="${TAG#v}"
      echo "tag=${TAG}" >> "$GITHUB_OUTPUT"
      echo "version=${VERSION}" >> "$GITHUB_OUTPUT"
      echo "Release gate for ${TAG} (${VERSION})"
    "%
    & { id = "meta" },
    C.run_script "Assert lockstep + CHANGELOG section" m%"
      VERSION="${{ steps.meta.outputs.version }}"
      python3 scripts/release_assert.py "$VERSION" --require-changelog
      NOTES="${RUNNER_TEMP}/release_notes.md"
      if [[ -f "docs/source/releases/v${VERSION}/release-notes.md" ]]; then
        cp "docs/source/releases/v${VERSION}/release-notes.md" "$NOTES"
      else
        python3 scripts/release_assert.py "$VERSION" --notes "$NOTES"
      fi
      echo "--- release notes (head) ---"
      head -40 "$NOTES" || true
    "%,
    {
      name = "Upload release notes",
      uses = C.upload_artifact_v4,
      with = {
        name = "release-notes-${{ github.sha }}",
        path = "${{ runner.temp }}/release_notes.md",
        "retention-days" = 5,
      },
    },
    C.run_script "Detect pre-release (RC) tag" m%"
      VERSION="${{ steps.meta.outputs.version }}"
      if [[ "$VERSION" == *-* ]]; then
        echo "is_prerelease=true" >> "$GITHUB_OUTPUT"
        echo "Pre-release: will skip PyPI (stable index only)."
      else
        echo "is_prerelease=false" >> "$GITHUB_OUTPUT"
      fi
    "%
    & { id = "prerelease" },
  ],
}
in

let job_tarball = {
  needs = ["gate"],
  "runs-on" = "ubuntu-latest",
  steps = [
    C.checkout_full,
    C.run_script "Build source tarball (git archive, not GitHub auto-archive)" m%"
      TAG="${{ needs.gate.outputs.tag }}"
      VERSION="${{ needs.gate.outputs.version }}"
      OUT="eon-v${VERSION}.tar.xz"
      git archive --format=tar "$TAG" | xz -9 > "$OUT"
      ls -la "$OUT"
      sha256sum "$OUT" | tee tarball.sha256
    "%,
    {
      name = "Upload tarball artifact",
      uses = C.upload_artifact_v4,
      with = {
        name = "eon-tarball-${{ needs.gate.outputs.version }}",
        path = m%"
eon-v${{ needs.gate.outputs.version }}.tar.xz
tarball.sha256
"%,
        "retention-days" = 30,
      },
    },
  ],
}
in

let job_gh_release = {
  needs = ["gate", "tarball"],
  "runs-on" = "ubuntu-latest",
  steps = [
    C.checkout_default,
    {
      name = "Download release notes",
      uses = C.download_artifact_v4,
      with = {
        name = "release-notes-${{ github.sha }}",
        path = "${{ runner.temp }}/release-notes",
      },
    },
    {
      name = "Download tarball",
      uses = C.download_artifact_v4,
      with = {
        name = "eon-tarball-${{ needs.gate.outputs.version }}",
        path = "${{ runner.temp }}/tarball",
      },
    },
    {
      name = "Create or update GitHub Release",
      env = {
        GH_TOKEN = "${{ github.token }}",
        TAG = "${{ needs.gate.outputs.tag }}",
        VERSION = "${{ needs.gate.outputs.version }}",
        IS_PRE = "${{ needs.gate.outputs.is_prerelease }}",
      },
      run = m%"
        set -euo pipefail
        ASSET="${RUNNER_TEMP}/tarball/eon-v${VERSION}.tar.xz"
        NOTES="${RUNNER_TEMP}/release-notes/release_notes.md"
        PRE_FLAG=()
        if [[ "$IS_PRE" == "true" ]]; then
          PRE_FLAG=(--prerelease)
        fi
        if gh release view "$TAG" --repo "${{ github.repository }}" >/dev/null 2>&1; then
          echo "Release $TAG exists; ensuring notes/asset..."
          gh release edit "$TAG" --repo "${{ github.repository }}" --notes-file "$NOTES" "${PRE_FLAG[@]}" || true
          gh release upload "$TAG" "$ASSET" --repo "${{ github.repository }}" --clobber
        else
          gh release create "$TAG" "$ASSET" \
            --repo "${{ github.repository }}" \
            --title "$TAG" \
            --notes-file "$NOTES" \
            "${PRE_FLAG[@]}"
        fi
      "%,
    },
  ],
}
in

let job_pypi = {
  needs = ["gate", "tarball", "gh-release"],
  "if" = "needs.gate.outputs.is_prerelease != 'true'",
  "runs-on" = "ubuntu-latest",
  environment = {
    name = "release",
    url = Pypi.pypi_project_url,
  },
  permissions = {
    contents = "read",
    "id-token" = "write",
  },
  steps = [
    C.checkout_full,
    C.setup_python C.python_version,
    C.run_script "Install build frontend" "python3 -m pip install --upgrade pip build",
    C.run_script "Build sdist/wheel (meson-python)" m%"
      python3 scripts/release_assert.py "${{ needs.gate.outputs.version }}" --require-changelog
      if python3 -m build 2>&1 | tee "${RUNNER_TEMP}/pypi-build.log"; then
        echo "build ok"
      else
        echo "::warning::python -m build failed (native deps missing on runner). Trying sdist via pip wheel isolation off."
        exit 1
      fi
      ls -la dist/
    "%,
    {
      name = "Publish to PyPI",
      uses = C.pypi_publish,
      with = {
        password = "${{ secrets.PYPI_API_TOKEN }}",
        "skip-existing" = true,
        verbose = true,
      },
    },
  ],
}
in

let job_pypi_skip_rc = {
  needs = ["gate", "tarball", "gh-release"],
  "if" = "needs.gate.outputs.is_prerelease == 'true'",
  "runs-on" = "ubuntu-latest",
  steps = [
    C.run_script "Skip PyPI for RC" m%"
      echo "Skipped PyPI for pre-release ${{ needs.gate.outputs.tag }}."
      echo "GitHub Release created as prerelease with tarball asset."
    "%,
  ],
}
in

{
  name = "Release",
  on = { push = { tags = ["v*"] } },
  permissions = {
    contents = "write",
    "id-token" = "write",
  },
  concurrency = {
    group = "release-${{ github.ref }}",
    "cancel-in-progress" = false,
  },
  jobs = {
    gate = job_gate,
    tarball = job_tarball,
    "gh-release" = job_gh_release,
    pypi = job_pypi,
    "pypi-skip-rc" = job_pypi_skip_rc,
  },
}
