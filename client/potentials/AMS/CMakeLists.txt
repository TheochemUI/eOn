find_package(absl CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)
find_package(Boost CONFIG REQUIRED COMPONENTS iostreams)

add_library(AMS SHARED AMS.cpp)
target_link_libraries(AMS PRIVATE eoncbase ${EON_DEPS} fmt::fmt absl::strings Boost::headers Boost::iostreams)
target_include_directories(AMS PRIVATE ${EON_INCLUDE_DIRS})
install(TARGETS AMS)
