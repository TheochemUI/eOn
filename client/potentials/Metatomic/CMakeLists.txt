add_library(Metatomic SHARED MetatomicPotential.cpp)
target_link_libraries(Metatomic PRIVATE eoncbase ${EON_DEPS})
target_include_directories(Metatomic PRIVATE ${EON_INCLUDE_DIRS})
install(TARGETS Metatomic)
