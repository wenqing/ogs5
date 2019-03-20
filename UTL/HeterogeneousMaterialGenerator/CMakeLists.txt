##-------------------------------------------------------
##
##    Tool to generate element wise heterogeneous material parameters 
##                By Wenqing Wang
##
##                     03.2019
##-------------------------------------------------------

if(COMPILER_IS_GCC)
	add_definitions(--std=c++14)
endif()

include_directories(
	${CMAKE_SOURCE_DIR}
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/FEM
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/FileIO/MeshIO
	${CMAKE_SOURCE_DIR}/GEO
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MSH
	${CMAKE_SOURCE_DIR}/MSHGEOTOOLS
	${CMAKE_SOURCE_DIR}/UTL/MSHGEOTOOLS
)

FILE(GLOB hete_p_H RELATIVE  ${CMAKE_CURRENT_SOURCE_DIR}   *.h)
FILE(GLOB hete_p_CPP RELATIVE  ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

add_executable(hete_mat ${hete_p_H} ${hete_p_CPP})

target_link_libraries( hete_mat
	MSH
	MSHGEOTOOLS
	FileIO
)
