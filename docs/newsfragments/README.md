# Changelog fragments (towncrier)

User-facing changes intended for the **next** semver should add a fragment here.
On `cog bump`, towncrier consumes these files into `CHANGELOG.md` and deletes them.

## Create a fragment

```bash
uvx towncrier create --content "Short user-facing note." 42.added.md
# or interactive:
uvx towncrier create
```

Filename pattern: `<issue_or_pr_number>.<type>.md` (or `+short_slug.<type>.md`).

## Types (`towncrier.toml`)

| Directory / type | Section in CHANGELOG |
|------------------|----------------------|
| `security` | Security |
| `removed` | Removed |
| `deprecated` | Deprecated |
| `added` | Added |
| `dev` | Developer |
| `changed` | Changed |
| `fixed` | Fixed |

## Preview

```bash
uvx towncrier build --draft
# or: python3 -m towncrier build --draft
```

## CI

PRs that touch shipped surfaces (`eon/`, `client/`, build files, `pyproject.toml`,
etc.) run `.github/workflows/towncrier.yml` (`towncrier check`). Pure documentation
or CI-only changes outside those paths may omit fragments.

See [docs/source/devdocs/release.md](../source/devdocs/release.md) for the full
release checklist (GitHub, PyPI, conda-forge/eon-feedstock).
