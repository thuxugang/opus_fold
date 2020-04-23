cmake_minimum_required(VERSION 3.5)
project(opus_fold)

set(CMAKE_CXX_STANDARD 11)

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Ofast")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Ofast")

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

find_package(OpenMP)
if (OPENMP_FOUND OR OpenMP_CXX_FOUND) 
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    message("USE OpenMP...")
endif()

message("=================================")
message(${CMAKE_C_FLAGS})
message(${CMAKE_CXX_FLAGS})

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
message("=================================")
message(${EIGEN3_INCLUDE_DIR})

find_package(libmongocxx REQUIRED)
find_package(libbsoncxx REQUIRED)
include_directories(${LIBMONGOCXX_INCLUDE_DIR})
include_directories(${LIBBSONCXX_INCLUDE_DIR})
include_directories("/usr/local/include/mongocxx/v_noabi")
include_directories("/usr/local/include/bsoncxx/v_noabi")
include_directories("/usr/local/include/libmongoc-1.0")
include_directories("/usr/local/include/libbson-1.0")
include_directories("/usr/local/lib")

message("=================================")
message(${LIBMONGOCXX_LIBRARIES})
message(${LIBBSONCXX_LIBRARIES})

link_directories(${PROJECT_SOURCE_DIR}/lib)
message("=================================")
message(${PROJECT_SOURCE_DIR}/lib)

aux_source_directory(./build DIR_SRCS)
add_executable(opus_fold ${DIR_SRCS})

target_link_libraries(opus_fold -lmkdssp)
target_link_libraries(opus_fold -lboost_date_time)
target_link_libraries(opus_fold -lboost_filesystem)
target_link_libraries(opus_fold -lboost_iostreams)
target_link_libraries(opus_fold -lboost_program_options)
target_link_libraries(opus_fold -lboost_regex)
target_link_libraries(opus_fold -lboost_system)
target_link_libraries(opus_fold -lboost_thread)
target_link_libraries(opus_fold -lpthread)

target_link_libraries(opus_fold ${LIBMONGOCXX_LIBRARIES})
target_link_libraries(opus_fold ${LIBBSONCXX_LIBRARIES})
