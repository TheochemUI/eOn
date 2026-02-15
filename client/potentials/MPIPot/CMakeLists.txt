add_library(MPIPot SHARED MPIPot.cpp)
target_link_libraries(MPIPot PRIVATE eoncbase ${EON_DEPS})
target_include_directories(MPIPot PRIVATE ${EON_INCLUDE_DIRS})
install(TARGETS MPIPot)
