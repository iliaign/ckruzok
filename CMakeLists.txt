cmake_minimum_required(VERSION 3.10)
project(kep_dec CXX)

find_package(Eigen3 REQUIRED NO_MODULE)


add_library(kep_dec_lib STATIC
    kep_dec/dec_to_kep.cpp
    kep_dec/dec_to_kep.hpp
    kep_dec/kep_to_dec.cpp
    kep_dec/kep_to_dec.hpp
)

target_include_directories(kep_dec_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/kep_dec)
target_link_libraries(kep_dec_lib PUBLIC Eigen3::Eigen)
target_compile_features(kep_dec_lib PRIVATE cxx_std_17)

add_executable(kep_dec main.cpp)
target_link_libraries(kep_dec PRIVATE kep_dec_lib)
target_compile_features(kep_dec PRIVATE cxx_std_17)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)