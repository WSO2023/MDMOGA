cmake_minimum_required(VERSION 3.14)
project(Frame)

include_directories(${PROJECT_SOURCE_DIR}/include)
set(CMAKE_CXX_STANDARD 14)
SET(CMAKE_CXX_FLAGS "-fopenmp")
add_subdirectory(${PROJECT_SOURCE_DIR}/lib)

add_executable(Frame main.cpp lib/pugixml.cpp include/pugixml.hpp  include/pugiconfig.hpp lib/classAndVarDefine.cpp include/classAndVarDefine.h include/common.h lib/config.cpp include/config.h lib/tools.cpp include/tools.h GRPHEFT.cpp GRPHEFT.h MOELS.cpp MOELS.h MDMOGA.cpp MDMOGA.h lib/GenOper.cpp include/GenOper.h lib/ChrOper.cpp include/ChrOper.h   EMS.cpp EMS.h  MOHEFT.cpp MOHEFT.h HMCMB.cpp HMCMB.h)