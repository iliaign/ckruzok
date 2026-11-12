cmake_minimum_required(VERSION 3.15...3.27)
project(pycosmos LANGUAGES CXX)


if (CMAKE_VERSION VERSION_LESS 3.18)
  set(DEV_MODULE Development)
else()
  set(DEV_MODULE Development.Module)
endif()
find_package(Python 3.9 COMPONENTS Interpreter ${DEV_MODULE} REQUIRED)

execute_process(
  COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
  OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_ROOT)
find_package(nanobind CONFIG REQUIRED)

# 3. Ваши существующие статические библиотеки
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


nanobind_add_module(pycosmos bindings/module.cpp)

target_link_libraries(pycosmos PRIVATE kep_dec_lib predict_lib)
target_include_directories(pycosmos PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

install(TARGETS pycosmos DESTINATION .)