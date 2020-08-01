include_directories(${CMAKE_CURRENT_LIST_DIR})

add_library(source OBJECT 
        calculation.cpp
        triangulation.cpp
        )

add_executable(start
            $<TARGET_OBJECTS:source>
            main.cpp
            )