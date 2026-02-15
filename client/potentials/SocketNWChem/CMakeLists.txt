add_library(SocketNWChem SHARED SocketNWChemPot.cpp)
target_link_libraries(SocketNWChem PRIVATE eoncbase ${EON_DEPS})
target_include_directories(SocketNWChem PRIVATE ${EON_INCLUDE_DIRS})
install(TARGETS SocketNWChem)
