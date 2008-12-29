SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/../CMakeModules )

# const-char* -> char* conversion, needed for Nana
IF( CMAKE_COMPILER_IS_GNUCXX )
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-write-strings -Wall" )
ENDIF()

include_directories ( 
	. 
	${CMAKE_SOURCE_DIR}/tools
	${CMAKE_SOURCE_DIR}/sdp_wrappers 
	${CMAKE_SOURCE_DIR}/seriation_gen
	${CMAKE_SOURCE_DIR}/prob_dt
	${CMAKE_SOURCE_DIR}/graph_gen
	)

add_library( hrdb SHARED 
	gdist_projected_db.cpp gdist_projected_db.hpp 
	db_cfg.cpp
	)
target_link_libraries(hrdb  ${Boost_LIBRARIES} ${NANA_LIBRARIES} )
