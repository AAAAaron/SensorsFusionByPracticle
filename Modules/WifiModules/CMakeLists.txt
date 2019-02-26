cmake_minimum_required( VERSION 2.8 )
project( WifiModule )

# 添加Eigen头文件
include_directories( "/usr/include/eigen3" )
set(SRC WifiModule.cpp WifiModule.h)                                                      
add_library(WifiModule SHARED ${SRC})  

#add_executable( eigenGeometry eigenGeometry.cpp )
# add_library(CoordinateTR SHARED CoordinateTR.cpp)
