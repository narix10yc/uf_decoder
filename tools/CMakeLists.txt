add_executable(scratch scratch.cpp)
target_link_libraries(scratch PRIVATE svg uf)

add_executable(scratch3d scratch3d.cpp)
target_link_libraries(scratch3d PRIVATE svg uf)