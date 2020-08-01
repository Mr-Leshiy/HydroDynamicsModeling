cmake_minimum_required(VERSION 3.5)


set(PROJECT_NAME HydroDynamcis)
set(VERSION 0.1)


project(${PROJECT_NAME} VERSION ${VERSION} LANGUAGES C CXX)

include_directories(include)
add_subdirectory(src)