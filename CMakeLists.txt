cmake_minimum_required(VERSION 2.8.4)
project(inclined)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_LIST_DIR}/cmake")

set(WARNING_FLAGS "-Wall -Wextra -Wsign-conversion -Wconversion -Wundef -Wcast-align -Wwrite-strings -Wlogical-op -Wmissing-declarations -Wredundant-decls -Wshadow -Woverloaded-virtual")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=c++11 ${WARNING_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

# Show make output (default is not)
option(VERBOSE_MAKE "VERBOSE_MAKE" OFF)
IF(VERBOSE_MAKE)
    set(CMAKE_VERBOSE_MAKEFILE ON)
ENDIF()

# Additional unsafe optimization (default is not)
option(ADDITIONAL_UNSAFE_OPTIMIZE "ADDITIONAL_UNSAFE_OPTIMIZE" OFF)
set(ADDITIONAL_OPTIMIZATION_FLAGS "-march=native -Ofast -msse3 -ffast-math -funroll-loops")
IF(ADDITIONAL_UNSAFE_OPTIMIZE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${ADDITIONAL_OPTIMIZATION_FLAGS}")
ENDIF()

IF(CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0")
ENDIF()

# MPI library
find_package(MPI REQUIRED)
include_directories(${MPI_CXX_INCLUDE_PATH})
# GSL library
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
# VTK library
find_package(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})
if(VTK_MAJOR_VERSION EQUAL 5)
    set(CONFIG_VTK_5 1)
elseif (VTK_MAJOR_VERSION EQUAL 6)
    set(CONFIG_VTK_6 1)
elseif (VTK_MAJOR_VERSION EQUAL 7)
    set(CONFIG_VTK_7 1)
else()
    message(FATAL_ERROR "Unsupported VTK version found: ${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION}")
endif()

# cmake generated config
set(CMAKE_CONFIG_SHARE_GCM "${CMAKE_INSTAL_PREFIX}/share/gcm")

set(CMAKE_CONFIG_ENABLE_LOGGING       1)
set(CMAKE_CONFIG_ENABLE_LOGGING_FATAL 1)
set(CMAKE_CONFIG_ENABLE_LOGGING_ERROR 1)
set(CMAKE_CONFIG_ENABLE_LOGGING_WARN  1)
set(CMAKE_CONFIG_ENABLE_LOGGING_INFO  1)
set(CMAKE_CONFIG_ENABLE_LOGGING_DEBUG 0)
set(CMAKE_CONFIG_ENABLE_LOGGING_TRACE 0)

IF(CMAKE_BUILD_TYPE MATCHES RelWithDebInfo)
    set(CMAKE_CONFIG_ENABLE_LOGGING_INFO  1)
ELSEIF(CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_CONFIG_ENABLE_LOGGING_INFO  1)
    set(CMAKE_CONFIG_ENABLE_LOGGING_DEBUG 1)
    set(CMAKE_CONFIG_ENABLE_LOGGING_TRACE 1)
ENDIF()

# source files
include_directories("${CMAKE_CURRENT_LIST_DIR}/")
FILE(GLOB_RECURSE HEADERS "${CMAKE_CURRENT_LIST_DIR}/src/*.hpp")
file(GLOB_RECURSE INCLINED_SOURCES "${CMAKE_CURRENT_LIST_DIR}/main.cpp")
file(GLOB_RECURSE TEST_SOURCES "${CMAKE_CURRENT_LIST_DIR}/test/*.cpp")

# inclined
add_executable(inclined ${INCLINED_SOURCES} ${HEADERS})
target_link_libraries(inclined ${MPI_CXX_LIBRARIES} ${GSL_LIBRARIES} ${VTK_LIBRARIES})

# tests
#add_executable(gcm_tests ${SEQUENCE_TEST_SOURCES})
#target_link_libraries(gcm_tests gcm cgal3dmesher ${GTEST_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${CGAL_LIBRARIES} ${GSL_LIBRARIES} ${LOG4CXX_LIBRARIES} ${MPI_CXX_LIBRARIES} ${VTK_LIBRARIES})

install(
    TARGETS inclined
    RUNTIME DESTINATION bin
)
