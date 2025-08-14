# UF Python Bindings

# This will put the dependencies in a top-level _deps directory
# so if we delete the build directory, we don't need to re-download.
set(FETCHCONTENT_BASE_DIR "${CMAKE_SOURCE_DIR}/_deps")

include(FetchContent)
FetchContent_Declare(
  pybind11
  GIT_REPOSITORY https://github.com/pybind/pybind11.git
  GIT_TAG v2.12.0
)

FetchContent_MakeAvailable(pybind11)

find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)
message(STATUS "Using Python: ${Python3_EXECUTABLE}")
set(Python3_FIND_VIRTUALENV FIRST)

file(GLOB SRCS "*.cpp")

pybind11_add_module(uf_python MODULE ${SRCS})

# copy .pyi file alongside the built .so file
set(PYI_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/uf_python.pyi)
set(PYI_TARGET ${CMAKE_CURRENT_BINARY_DIR}/uf_python.pyi)

add_custom_command(
  OUTPUT ${PYI_TARGET}
  DEPENDS ${PYI_SOURCE}
  COMMAND ${CMAKE_COMMAND} -E copy_if_different ${PYI_SOURCE} ${PYI_TARGET}
  COMMENT "Copying uf_python.pyi to output directory"
)

add_custom_target(uf_python_copy_pyi ALL DEPENDS ${PYI_TARGET})
add_dependencies(uf_python uf_python_copy_pyi)
target_link_libraries(uf_python PRIVATE uf svg)