project(bezmesh)

### Find dependencies
find_package(GMP REQUIRED)
find_package(GMPXX REQUIRED)

OPTION (ALL_EXACT "Exact Algorithm" OFF)
IF(ALL_EXACT)
    add_definitions( -DALL_EXACT )
    message("Using Exact Version: needs CGAL")
    find_package(CGAL COMPONENTS)
    if(NOT CGAL_FOUND)
        message(FATAL_ERROR "CGAL NOT found.")
    ENDIF()
    include( ${CGAL_USE_FILE} )
    set(CMAKE_BUILD_TYPE "Release")
ENDIF()

### Register all sources
file(GLOB_RECURSE SOURCE_LIST "*.cpp")
add_library(${PROJECT_NAME} ${SOURCE_LIST})

### Link with dependencies
target_link_libraries(${PROJECT_NAME}
    triangle
    aabbcc
    ${GMP_LIBRARIES}
    ${GMPXX_LIBRARIES}
)

## EXACT Version: Link with CGAL
IF(ALL_EXACT)
    target_link_libraries( ${PROJECT_NAME} ${CGAL_LIBRARIES} ${CGAL_3RD_PARTY_LIBRARIES})
ENDIF()
