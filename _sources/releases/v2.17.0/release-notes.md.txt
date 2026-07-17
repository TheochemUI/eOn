# eOn 2.17.0

Class-first **pyeonclient** Matter/NEB surface, shared **eon-schema** L1/L2
(including expanded NebSpec), metatomic backends (fat / RGPOT / ASE), and GIL
safety for torch autograd from C++.

See monorepo `CHANGELOG.md` section 2.17.0 for the full fragment list.

## Companion split packages (same train)

| Package | Version | Tag |
|---------|---------|-----|
| eon-schema | 0.2.1 | `eon-schema-v0.2.1` |
| pyeonclient | 0.3.2 | `pyeonclient-v0.3.2` |
| eon / eon-akmc | 2.17.0 | `v2.17.0` |

Publish order: eon-schema → pyeonclient → monorepo fat tag.
