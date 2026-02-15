add_library(ExtPot SHARED ExtPot.cpp)
target_link_libraries(ExtPot PRIVATE eoncbase ${EON_DEPS})
target_include_directories(ExtPot PRIVATE ${EON_INCLUDE_DIRS})
install(TARGETS ExtPot)
