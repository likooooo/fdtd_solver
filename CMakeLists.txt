cmake_minimum_required(VERSION 3.10) # 根据你的需求调整最低CMake版本
project(Meep VERSION 1.27.0 LANGUAGES C CXX Fortran) # 假设Meep也用Fortran

# 共享库 (AC_ENABLE_SHARED)
option(BUILD_SHARED_LIBS "Build shared libraries" ON) # CMake默认就是ON，但明确设置更好

# C++11 支持
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# fftw
find_package(PkgConfig REQUIRED)
pkg_search_module(FFTW REQUIRED fftw3 IMPORTED_TARGET)
include_directories(PkgConfig::FFTW )
link_libraries     (PkgConfig::FFTW)
pkg_search_module(FFTWF REQUIRED fftw3f IMPORTED_TARGET)
include_directories(PkgConfig::FFTWF )
link_libraries     (PkgConfig::FFTWF)
add_definitions(-DHAVE_LIBFFTW3=1)
# mkl
set(MKL_THREADING gnu_thread)
set(MKL_MPI openmpi)
set(MKL_LINK static)
if(EXISTS "/opt/intel/oneapi/mkl/latest/lib/cmake/mkl/")
    set(MKL_DIR "/opt/intel/oneapi/mkl/latest/lib/cmake/mkl/")
    set(MKL_ARCH intel64)
    set(MKL_INTERFACE lp64)
    find_package(MKL CONFIG REQUIRED)
else()
    find_package(MKL CONFIG REQUIRED)
endif()
# hdf5
pkg_search_module(HDF5 REQUIRED hdf5 IMPORTED_TARGET)
include_directories(PkgConfig::HDF5 )
link_libraries     (PkgConfig::HDF5)

add_subdirectory(src)
add_subdirectory(tests)
