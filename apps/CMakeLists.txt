project(apps)

### Define local paths
set(RESOURCE_PATH ${CMAKE_SOURCE_DIR}/resource/ CACHE PATH "location of local resources")
add_library(localpaths INTERFACE)
target_compile_definitions(localpaths INTERFACE
    "-DRESOURCE_PATH=\"${RESOURCE_PATH}\"")
target_include_directories(localpaths INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})

### Add all applications
add_subdirectory(bezmeshCLI)
