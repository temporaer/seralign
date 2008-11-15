SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/../CMakeModules )

# const-char* -> char* conversion, needed for Nana
IF( CMAKE_COMPILER_IS_GNUCXX )
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-write-strings -Wall" )
ENDIF()

include_directories( 
	. 
	${CMAKE_SOURCE_DIR}/configuration
	${CMAKE_SOURCE_DIR}/sdp_wrappers
	)

add_library( prob_dt SHARED 
	ProbData.cpp  ProbData.hpp 
	AdjMat.cpp    AdjMat.hpp 
	SimMat.cpp    SimMat.hpp 
	PermMat.cpp   PermMat.hpp 
	Laplacian.cpp Laplacian.hpp)
