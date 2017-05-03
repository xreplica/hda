include_directories("${PROJECT_SOURCE_DIR}" "${PROJECT_SOURCE_DIR}/src/core")

add_library(cluster OBJECT cluster.cpp)
add_library(distance OBJECT distance.cpp)
add_library(matrix OBJECT matrix.cpp)
add_library(params OBJECT params.cpp)
add_library(point OBJECT point.cpp)
add_library(validation OBJECT validation.cpp)
