project(garment_transfer)
add_executable(${PROJECT_NAME} main.cpp)
target_compile_definitions(${PROJECT_NAME} PUBLIC "-DDATA_PATH=\"${DATA_PATH}\"")
target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen pmp_jie optimizer CGAL::CGAL)