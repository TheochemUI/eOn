add_library(AMS_IO SHARED AMS_IO.cpp)
target_link_libraries(AMS_IO PRIVATE eoncbase ${EON_DEPS})
target_include_directories(AMS_IO PRIVATE ${EON_INCLUDE_DIRS})
install(TARGETS AMS_IO)
