py.install_sources(
    [
        'akmc.py',
        'akmcstatelist.py',
        'akmcstate.py',
        'askmc.py',
        'atoms.py',
        'basinhopping.py',
        'communicator.py',
        'config.py',
        'config.yaml',
        'displace.py',
        'eon_kdb.py',
        'escaperate.py',
        'explorer.py',
        'fileio.py',
        '__init__.py',
        'locking.py',
        'server.py',
        'schema.py',
        # TODO(rg): Spin mcamc out and use as a dependency
        # This will only use the python variant for now
        'mcamc/__init__.py',
        'mcamc/mcamc.py',
        'movie.py',
        'mpiwait.py',
        'parallelreplica.py',
        'prstatelist.py',
        'prstate.py',
        'recycling.py',
        'statelist.py',
        'state.py',
        'superbasin.py',
        'superbasinscheme.py',
    ],
    pure: false,  # install next to compiled extension
    subdir: 'eon',
)


eon_dir = py.get_install_dir() / 'eon'

# Generate version.py for sdist
meson.add_dist_script(
    ['../tools/gitversion.py', '--meson-dist', '--write', 'eon/version.py'],
)
fs = import('fs')
if not fs.exists('version.py')
    generate_version = custom_target(
        'generate-version',
        install: true,
        build_always_stale: true,
        build_by_default: true,
        output: 'version.py',
        input: '../tools/gitversion.py',
        command: [
            '../tools/gitversion.py',
            '--meson-dist',
            '--write',
            'eon/version.py',
        ],
        install_dir: eon_dir,
        install_tag: 'python-runtime',
    )
else
    py.install_sources(['version.py'], subdir: 'eon')
endif
