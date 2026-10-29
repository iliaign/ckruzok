cmake_minimum_required(VERSION 3.10)
project(kep_dec CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Eigen3 REQUIRED NO_MODULE)

add_executable(kep_dec 
    main.cpp 
    kep_dec/dec_to_kep.cpp 
    kep_dec/kep_to_dec.cpp
    kep_dec/dec_to_kep.hpp
    kep_dec/kep_to_dec.hpp
)

target_include_directories(kep_dec PRIVATE kep_dec)

target_link_libraries(kep_dec Eigen3::Eigen)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)