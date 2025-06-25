include_directories(../src)

file(GLOB src ${CMAKE_CURRENT_LIST_DIR}/*.cpp)
foreach(s IN LISTS src)
    get_filename_component(test_WLE ${s} NAME_WLE)
    add_executable(${test_WLE} ${s})
    target_link_libraries(${test_WLE} PUBLIC meep)
endforeach()