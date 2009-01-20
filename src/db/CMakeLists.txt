SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/../CMakeModules )

# const-char* -> char* conversion, needed for Nana
IF( CMAKE_COMPILER_IS_GNUCXX )
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-write-strings -Wall" )
ENDIF()

SET(CMAKE_INCLUDE_PATH ${CMAKE_INCLUDE_PATH} ${CMAKE_SOURCE_DIR}/third_party/ann_1.1.1/include )
SET(CMAKE_LIBRARY_PATH ${CMAKE_LIBRARY_PATH} ${CMAKE_SOURCE_DIR}/third_party/ann_1.1.1/lib )
FIND_PACKAGE(ANN)

include_directories ( 
	. 
	${CMAKE_SOURCE_DIR}/tools
	${CMAKE_SOURCE_DIR}/sdp_wrappers 
	${CMAKE_SOURCE_DIR}/seriation_gen
	${CMAKE_SOURCE_DIR}/prob_dt
	${CMAKE_SOURCE_DIR}/graph_gen
	${ANN_INCLUDE_DIR}
	)


add_library( hrdb SHARED 
	icpdb.cpp icpdb.hpp 
	anndb.cpp anndb.hpp 
	db_cfg.cpp
	)
target_link_libraries(hrdb  ${Boost_LIBRARIES} ${NANA_LIBRARIES} ${ANN_LIBRARIES} )
