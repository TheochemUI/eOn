add_library(ZBL SHARED ZBLPot.cpp)
target_link_libraries(ZBL PRIVATE eoncbase ${EON_DEPS})
target_include_directories(ZBL PRIVATE ${EON_INCLUDE_DIRS})
install(TARGETS ZBL)
