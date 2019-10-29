
include_directories(
	${CMAKE_SOURCE_DIR}
	${CMAKE_SOURCE_DIR}/Base
	${CMAKE_SOURCE_DIR}/FileIO
	${CMAKE_SOURCE_DIR}/FileIO/MeshIO
	${CMAKE_SOURCE_DIR}/FEM
        ${CMAKE_SOURCE_DIR}/GEO
	${CMAKE_SOURCE_DIR}/MathLib
	${CMAKE_SOURCE_DIR}/MSH
	${CMAKE_SOURCE_DIR}/MSHGEOTOOLS
	${CMAKE_SOURCE_DIR}/UTL/VariableInterpolation
)

FILE(GLOB VariabbleInterpolation_H RELATIVE  ${CMAKE_CURRENT_SOURCE_DIR}   *.h)
FILE(GLOB VariabbleInterpolation_CPP RELATIVE  ${CMAKE_CURRENT_SOURCE_DIR} *.cpp)

SET (VariabbleInterpolation_H, ${VariabbleInterpolation_H},
              ${CMAKE_SOURCE_DIR}/FEM/fem_ele.h
              ${CMAKE_SOURCE_DIR}/MSH/MSHEnums.h)
SET (VariabbleInterpolation_CPP, ${VariabbleInterpolation_CPP},
           ${CMAKE_SOURCE_DIR}/FEM/fem_ele.cpp,
           ${CMAKE_SOURCE_DIR}/MSH/MSHEnums.cpp)

add_executable(extract_point_value ${VariabbleInterpolation_CPP}
                                   ${VariabbleInterpolation_H})

target_link_libraries( extract_point_value MSH)

if (WIN32)
	if (MSVC)
		target_link_libraries( extract_point_value
			winmm.lib # for timer
		)
	endif()
endif()
