add_library(Morse SHARED Morse.cpp)
target_link_libraries(Morse PRIVATE eoncbase ${EON_DEPS})
target_include_directories(Morse PRIVATE ${EON_INCLUDE_DIRS})
install(TARGETS Morse)
