file(GLOB_RECURSE SOURCES "SVG/*.cpp")
add_library(svg STATIC ${SOURCES})

file(GLOB_RECURSE SOURCES "UF/*.cpp")
message(STATUS "UF sources: ${SOURCES}")
add_library(uf STATIC ${SOURCES})
