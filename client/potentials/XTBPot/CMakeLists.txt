add_library(XTBPot SHARED XTBPot.cpp)
target_link_libraries(XTBPot PRIVATE eoncbase ${EON_DEPS})
target_include_directories(XTBPot PRIVATE ${EON_INCLUDE_DIRS})
install(TARGETS XTBPot)
