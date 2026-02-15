add_library(CatLearnPot SHARED CatLearnPot.cpp)
target_link_libraries(CatLearnPot PRIVATE eoncbase ${EON_DEPS})
target_include_directories(CatLearnPot PRIVATE ${EON_INCLUDE_DIRS})
install(TARGETS CatLearnPot)
