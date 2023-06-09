project(bezmesh)


OPTION (ALL_EXACT "Exact Algorithm" OFF)
IF(ALL_EXACT)
    add_definitions( -DALL_EXACT )
    message("Using Exact Version: needs CGAL")
    include(cgal)
ELSE()
    include(gmp)
ENDIF()

### Register all sources
file(GLOB_RECURSE SOURCE_LIST "*.cpp")
add_library(${PROJECT_NAME} ${SOURCE_LIST})

### Link with dependencies
target_link_libraries(${PROJECT_NAME}
    triangle
    aabbcc
    gmp::gmp
)

## EXACT Version: Link with CGAL
IF(ALL_EXACT)
    target_link_libraries( ${PROJECT_NAME} CGAL::CGAL CGAL::CGAL_Core)
ENDIF()
