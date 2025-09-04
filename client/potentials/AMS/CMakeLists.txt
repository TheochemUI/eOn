set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)

# For the AMS interface
find_package(absl CONFIG REQUIRED)
find_package(fmt CONFIG REQUIRED)
find_package(Boost CONFIG REQUIRED COMPONENTS iostreams)

set(DEPS fmt::fmt absl::strings Boost::headers Boost::iostreams)

add_library(
  AMS SHARED
  AMS.cpp
  )

target_link_libraries(AMS INTERFACE ${DEPS})

# Install
install(TARGETS AMS)
