---
myst:
  html_meta:
    "description": "The eOn release workflow: cog-driven version bump, towncrier changelog, downstream cookbook gate, GitHub release, and conda-forge feedstock update."
    "keywords": "eOn release, cocogitto, towncrier, conventional commits, conda-forge feedstock, atomistic-cookbook"
---

# Release workflow

eOn releases use [cocogitto (`cog`)](https://docs.cocogitto.io/) to derive the
next semver tag from the conventional-commit log, run the mechanical bump
steps, create an atomic signed tag, and emit the source tarball. The
`cog.toml` at the repo root wires this up; the hand checklist below documents
what it automates and the manual steps that precede and follow it.

## 1. Pre-release checklist

Before invoking `cog bump`, verify:

- [ ] Every PR intended for the release is merged into `main`.
- [ ] `CHANGELOG.md` is up to date *through the previous release only* --
  towncrier will prepend the new section from `docs/newsfragments/*`. If
  fragments are missing for user-facing changes, add them first.
- [ ] `docs/source/releases/v<X.Y.Z>/release-notes.md` and `index.md` exist
  and are linked from `docs/source/releases/index.md` (cog cannot generate
  these -- they are curated).
- [ ] CI is green on `main`.
- [ ] `pixi install --locked` succeeds cleanly (so the lockfile regeneration
  hook has something stable to work from).
- [ ] **Downstream integration test passes.** The
  [lab-cosmo/atomistic-cookbook](https://github.com/lab-cosmo/atomistic-cookbook)
  example `eon-pet-neb` exercises eOn + PET-MAD on the oxadiazole NEB and
  is the canonical smoke gate for metatomic consumers. Point a cookbook
  checkout at the prospective release commit of eOn (edit its
  `environment.yml` to override the `eon` dependency, or
  `pip install -e` the eOn source tree into the nox-managed env) and run:

  ``` bash
  nox -e eon-pet-neb
  ```

  A pass means the nox session exits 0, the generated gallery artefacts
  regenerate, and `client_traceback.log` inside the example directory is
  empty. A regression here blocks the release.

## 2. Running `cog bump`

The one-command happy path on `main`:

``` bash
cog bump --auto --annotated 'Release v<X.Y.Z>

<one-paragraph body listing headline features, e.g.:
ARTn saddle search (pARTn), OCINEB hybrid CI-NEB + min-mode,
IRA structure comparison, ...>
See CHANGELOG.md.'
```

Override cog's commit-log inference with `--minor`, `--major`, or
`--version X.Y.Z` when needed (e.g. the only `feat:` commits since the
last tag are `feat(docs)` and a patch bump is incorrect).

### What `cog bump` does

Defined in `cog.toml`:

1. Picks the next version (semver increment inferred from conventional
   commit types since the last tag, or taken from the flag).
2. Runs `pre_bump_hooks`:
   - `sed` the new version into `pyproject.toml` and `pixi.toml`. All
     other consumers (`docs/source/conf.py`, `client/get_version.py`,
     the meson and CMake builds) read from `pyproject.toml`.
   - `uvx pixi-to-conda-lock pixi.lock --output condaEnvs` regenerates
     every `condaEnvs/*.conda-lock.yml`.
   - `pixi r -e dev towncrier build --version {{version}} --date
     $(date +%Y-%m-%d) --yes` consumes every file in
     `docs/newsfragments/` and prepends the new section to
     `CHANGELOG.md`.
3. Stages the resulting diff and writes the atomic release commit
   (`chore(version): <X.Y.Z>` by default; override via the
   `--annotated` body if you want custom wording).
4. Creates the signed annotated tag `v<X.Y.Z>`.
5. Runs `post_bump_hooks`:
   - `git archive --format=tar v<X.Y.Z> | xz -9 > eon-v<X.Y.Z>.tar.xz`
     produces the source tarball. Automatic GitHub tag archives are
     avoided because their contents have a
     [documented stability caveat](https://github.blog/open-source/git/update-on-the-future-stability-of-source-code-archives-and-hashes/).
   - Prints a reminder for the remaining manual steps.

`cog.toml` sets `disable_changelog = true` so cog stays out of
`CHANGELOG.md` (towncrier owns it) and `branch_whitelist = ["main"]`
so `cog bump` refuses to run on a PR branch.

## 3. Post-bump actions

1. Inspect the release commit:

   ``` bash
   git show HEAD
   ```

   The diff should look like `1edaf4d8` from v2.12.0 -- version bumps
   in `pyproject.toml`/`pixi.toml`, regenerated `condaEnvs/*.conda-lock.yml`,
   an appended block in `CHANGELOG.md`, every
   `docs/newsfragments/*` file deleted, and the new release-notes page.

2. Push the commit and tag together:

   ``` bash
   git push origin main --follow-tags
   ```

3. Create the GitHub release with the tarball attached:

   ``` bash
   gh release create v<X.Y.Z> eon-v<X.Y.Z>.tar.xz \
     --title v<X.Y.Z> \
     --notes-file docs/source/releases/v<X.Y.Z>/release-notes.md
   ```

   The feedstock URL assumes this exact tarball filename; do not rename
   the asset.

## 4. conda-forge feedstock bump

Once the GitHub release is published (so the asset URL resolves), update
the feedstock at
[conda-forge/eon-feedstock](https://github.com/conda-forge/eon-feedstock):

1. Branch `release-<X.Y.Z>` off `main`.
2. In `recipe/recipe.yaml`, bump `context.version` to `"<X.Y.Z>"`.
3. Compute the tarball sha256 and replace the first `source[].sha256`:

   ``` bash
   curl -sSL "https://github.com/TheochemUI/eOn/releases/download/v<X.Y.Z>/eon-v<X.Y.Z>.tar.xz" \
     | sha256sum
   ```

4. Reset `build.number: 0`.
5. Verify `fix_vesin_const.patch` and `fix_capnpc.patch` still apply
   cleanly; if not, regenerate them against the new source tree.
6. Commit, push, `gh pr create` against `conda-forge/eon-feedstock`.
7. Watch the Azure Pipelines matrix (linux-64, osx-64, osx-arm64, win-64)
   via `cf-ci watch conda-forge/eon-feedstock#<N>`.

## 5. Manual fallback

When `cog bump` fails mid-run (hook script error, unexpected prompt, etc.)
the state is recoverable: no commit has been written yet. Fix the root
cause and re-run. If you must cut by hand, this is the equivalent
sequence; `1edaf4d8 chore: release v2.12.0` is the canonical reference.

``` bash
# Version
sed -i 's/^version = "[^"]*"/version = "<X.Y.Z>"/' pyproject.toml pixi.toml

# Lockfiles
uvx pixi-to-conda-lock pixi.lock --output condaEnvs

# Changelog
pixi r -e dev towncrier build --version <X.Y.Z> \
  --date "$(date +%Y-%m-%d)" --yes

# Release notes page (manual)
# Create docs/source/releases/v<X.Y.Z>/{index.md,release-notes.md} and
# add to docs/source/releases/index.md toctree.

# Atomic commit + signed tag
git add -A
git commit -m "chore: release v<X.Y.Z>"
git tag -s -a v<X.Y.Z> -m "Release v<X.Y.Z>

<body>"

# Tarball
git archive --format=tar v<X.Y.Z> | xz -9 > eon-v<X.Y.Z>.tar.xz
```

Then continue from [§3 post-bump actions](#post-bump-actions).

## 6. Paper tags

Preprints and other paper-related tags do not need the release
machinery; just a lightweight annotated tag:

``` bash
git tag -a arXiv_2510.06030v3 45c077cb
# feat(arxiv): descriptive one-line summary
#
# body, including the arXiv link
```
