file(GLOB src ${CMAKE_CURRENT_LIST_DIR}/*.cpp)
add_library(meep SHARED ${src})
target_link_libraries(meep PUBLIC MKL::MKL ctl ctlgeom)