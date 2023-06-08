project(bezmeshCLI)

add_executable(${PROJECT_NAME} main.cpp)
target_link_libraries(${PROJECT_NAME} localpaths bezmesh triangle)
