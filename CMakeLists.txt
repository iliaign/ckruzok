cmake_minimum_required(VERSION 3.10)
project(kep_dec CXX)

find_package(Eigen3 REQUIRED NO_MODULE)


add_library(kep_dec_lib STATIC
    kep_dec/dec_to_kep.cpp
    kep_dec/dec_to_kep.hpp
    kep_dec/kep_to_dec.cpp
    kep_dec/kep_to_dec.hpp
    kep_dec/kepler.hpp
    kep_dec/kep_dec.hpp
)

add_library(predict_lib STATIC
    predict/ean_mean.cpp
    predict/mean_motion.cpp
    predict/solve_kepler.cpp
    predict/ecsrt_v_ist.cpp
    predict/ist_v_ecrt.cpp
    predict/predict.cpp
    predict/predict.hpp
)
target_include_directories(predict_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/predict)
target_include_directories(kep_dec_lib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/kep_dec)

target_link_libraries(kep_dec_lib PUBLIC Eigen3::Eigen)
target_link_libraries(predict_lib PUBLIC kep_dec_lib)

add_executable(kep_dec main.cpp)
target_link_libraries(kep_dec PRIVATE kep_dec_lib predict_lib)
target_compile_features(kep_dec PRIVATE cxx_std_17)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)