file(GLOB_RECURSE SOURCES "SVG/*.cpp")
add_library(svg ${SOURCES})

file(GLOB_RECURSE SOURCES "UF/*.cpp")
add_library(uf ${SOURCES})
