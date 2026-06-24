# Security Policy

## Supported Versions

Security fixes are applied on a best-effort basis to the latest released semver
on the `main` branch and the corresponding GitHub Release / conda-forge build.

| Version | Supported |
|---------|-----------|
| Latest `2.x` release on `main` / PyPI / conda-forge | Yes |
| Older minor/patch lines | Best effort only; please upgrade |
| Pre-release / RC tags (`X.Y.Z-rc.N`) | No stability or security SLA |

See [docs/source/devdocs/release.md](docs/source/devdocs/release.md) for how
releases are cut (cog, towncrier, GitHub, PyPI, conda-forge/eon-feedstock).

## Reporting a Vulnerability

**Do not** open a public GitHub issue for unreleased security problems.

Preferred:

1. Email the maintainers listed in `pyproject.toml` / feedstock maintainers with
   a clear description, impact, and reproduction if possible.
2. Or use GitHub **Private vulnerability reporting** on
   [TheochemUI/eOn](https://github.com/TheochemUI/eOn) if enabled for the repo.

Please allow reasonable time for a fix and coordinated disclosure before public
discussion. Credit is given in release notes / `CHANGELOG.md` (`security`
towncrier fragments) when reporters wish to be named.

## Supply chain notes

- Source of truth for a cut is the annotated tag `vX.Y.Z` and the
  `eon-vX.Y.Z.tar.xz` asset (git archive), not ad-hoc tarballs from forks.
- PyPI `eon` publishes only from `release.yml` on stable tags (trusted publisher
  or scoped token); RC tags skip PyPI.
- conda-forge builds are reviewed via the feedstock PR process separately.
