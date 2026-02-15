add_library(QSC SHARED QSC.cpp)
target_link_libraries(QSC PRIVATE eoncbase ${EON_DEPS})
target_include_directories(QSC PRIVATE ${EON_INCLUDE_DIRS})
install(TARGETS QSC)
