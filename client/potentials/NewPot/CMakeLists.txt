add_library(NewPot SHARED NewPot.cpp)
target_link_libraries(NewPot PRIVATE eoncbase ${EON_DEPS})
target_include_directories(NewPot PRIVATE ${EON_INCLUDE_DIRS})
install(TARGETS NewPot)
