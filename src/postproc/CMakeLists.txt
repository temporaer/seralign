SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/../CMakeModules )

# const-char* -> char* conversion, needed for Nana
IF( CMAKE_COMPILER_IS_GNUCXX )
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-write-strings -Wall" )
ENDIF()

include_directories( 
	. 
	${CMAKE_SOURCE_DIR}/configuration
	${CMAKE_SOURCE_DIR}/graph_gen
	${CMAKE_SOURCE_DIR}/seriation_gen
	${CMAKE_SOURCE_DIR}/prob_dt
	${CMAKE_SOURCE_DIR}/tools
	${CMAKE_SOURCE_DIR}/third_party/seqan
	)

SET(SRC 
	postproc.cpp     postproc.hpp
	RealXMLPrint.cpp RealXMLPrint.hpp
	GraphVizPrint.cpp GraphVizPrint.hpp
	PlainPrinter.cpp   PlainPrinter.hpp
	localign_distmat.cpp localign_distmat.hpp
		)

add_library( postproc SHARED ${SRC} )
