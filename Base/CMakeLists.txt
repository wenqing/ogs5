# Source files
set( HEADERS
	binarySearch.h
	BuildInfo.h
	Configure.h.in
        display.h
	DateTools.h
	FileFinder.h
	FileTools.h
	MemWatch.h
	quicksort.h
	printList.h
	uniqueListInsert.h
	StringTools.h
	swap.h
	wait.h	
)

set( SOURCES
	DateTools.cpp
	MemWatch.cpp
        ${CMAKE_CURRENT_BINARY_DIR}/BuildInfo.cpp
	StringTools.cpp
	binarySearch.cpp
        display.cpp
)

# Create the library
add_library( Base STATIC ${HEADERS} ${SOURCES} )

set_target_properties(Base PROPERTIES LINKER_LANGUAGE CXX)

include_directories(
	${CMAKE_SOURCE_DIR}/Base
)

if(OGS_BUILD_CORRECT_SCIENTIFIC)
	# Create utility executable
	add_executable( correctScientificNotation CorrectScientificNotationMain.cpp )
endif()
