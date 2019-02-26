cmake_minimum_required( VERSION 2.8 )
project( Practicle_Filter )

# 添加Eigen头文件
include_directories( "/usr/include/eigen3" )
set(SRC Practicle_Filter.cpp  Practicle.cpp WifiModule.cpp LeastSquare.cpp Practicle_Filter.h)                                                      
add_library(Practicle_Filter SHARED ${SRC})  

#add_executable( eigenGeometry eigenGeometry.cpp )
# add_library(CoordinateTR SHARED CoordinateTR.cpp)
