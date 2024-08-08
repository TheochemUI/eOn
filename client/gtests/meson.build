test_args = _args
test_deps = _deps
test_array = [  #
    ['test_element', 'ElementTest.cc', ''],
    # ['test_matter_cache', 'MatterCache.cc', 'sulfolene'],
    # ['test_enum', 'EnumTest.cc', ''],
    # ['test_iohelp', 'IOHelpersTest.cc', ''],
    # ['test_potcount', 'PotCountTest.cc', '../saddle_search'],
    # # ['test_impldim', 'ImpDimerTest.cpp', '/gtests/data/saddle_search'],
    # # ['strparse_run', 'StringHelpersTest.cpp', ''],
    # # ['test_matter', 'MatterTest.cpp', 'sulfolene'],
    # # ['test_pot', 'PotTest.cpp', 'sulfolene'],
    # # ['approve_params', 'ApproveParams.cpp', ''],
    # ['approve_matter', 'ApproveMatter.cpp', 'sulfolene'],
    # ['approve_matter_writers', 'ApproveMatterWriters.cc', 'sulfolene'],
    # # ['neb_test', 'NEBTest.cpp', '/gtests/data/saddle_search'],
]
foreach test : test_array
    test(
        test.get(0),
        executable(
            test.get(0),
            sources: [test.get(1), catch_cpp],
            dependencies: test_deps,
            include_directories: _incdirs,
            cpp_args: test_args,
            link_with: _linkto,
        ),
        workdir: meson.project_source_root() + '/client/gtests/data/systems/' + test.get(
            2,
        ),
    )
endforeach
