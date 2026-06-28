rgpot_pot = library(
    'rgpot_pot',
    sources: ['RgpotPot.cpp', 'RgpotClientRpc.cpp'],
    include_directories: _incdirs,
    dependencies: _deps + [ptlrpc_dep, dependency('capnp-rpc')],
    link_with: _linkto,
    install: true,
)
