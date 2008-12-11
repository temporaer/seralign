# vim:ft=cmake
ENABLE_TESTING()

SET(Boost_USE_STATIC_LIBS ON)                                                                                                                    
find_package( Boost 1.35 COMPONENTS unit_test_framework program_options filesystem ) 
IF(DSDP_FOUND)
  MESSAGE(STATUS "Found DSDP, wrapping it..." )
  SET(SRC ${SRC} dsdp_wrapper.cpp dsdp_wrapper.hpp )
  ADD_DEFINITIONS(-DDSDP_BINARY="${DSDP_BINARY}")
  ADD_DEFINITIONS(-DHAVE_DSDP)
ENDIF(DSDP_FOUND)
IF(CSDP_FOUND)
  MESSAGE(STATUS "Found CSDP, wrapping it..." )
  SET(SRC ${SRC} csdp_wrapper.cpp csdp_wrapper.hpp)
  ADD_DEFINITIONS(-DHAVE_CSDP)
  ADD_DEFINITIONS(-DCSDP_BINARY="${CSDP_BINARY}")
ENDIF(CSDP_FOUND)
IF(SDPA_FOUND)
  MESSAGE(STATUS "Found SDPA, wrapping it..." )
  SET(SRC ${SRC} sdpa_wrapper.cpp sdpa_wrapper.hpp)
  ADD_DEFINITIONS(-DSDPA_BINARY="${SDPA_BINARY}")
  ADD_DEFINITIONS(-DHAVE_SDPA)
  ADD_DEFINITIONS(-DSDPA_PARAM_FILE="param.sdpa" )
  CONFIGURE_FILE(${SDPA_PATH}/param.sdpa ${CMAKE_CURRENT_BINARY_DIR}/param.sdpa COPYONLY)
  CONFIGURE_FILE(${SDPA_PATH}/param.sdpa ${CMAKE_CURRENT_BINARY_DIR}/tests/param.sdpa COPYONLY)
ENDIF(SDPA_FOUND)
IF(SDPLR_FOUND)
  MESSAGE(STATUS "Found SDPLR, wrapping it..." )
  SET(SRC ${SRC} sdplr_wrapper.cpp sdplr_wrapper.hpp)
  ADD_DEFINITIONS(-DSDPLR_BINARY="${SDPLR_BINARY}")
  ADD_DEFINITIONS(-DHAVE_SDPLR)
  ADD_DEFINITIONS(-DSDPLR_PARAM_FILE="sdplr.params" )
  CONFIGURE_FILE(${SDPLR_PATH}/sdplr.params ${CMAKE_CURRENT_BINARY_DIR}/sdplr.params COPYONLY)
  CONFIGURE_FILE(${SDPLR_PATH}/sdplr.params ${CMAKE_CURRENT_BINARY_DIR}/tests/sdplr.params COPYONLY)
ENDIF(SDPLR_FOUND)

# const-char* -> char* conversion, needed for Nana
IF( CMAKE_COMPILER_IS_GNUCXX )
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-write-strings -Wall" )
ENDIF()

include_directories( 
	. 
	${CMAKE_SOURCE_DIR}/configuration
	${CMAKE_SOURCE_DIR}/graph_gen
	${CMAKE_SOURCE_DIR}/prob_dt
	${CMAKE_SOURCE_DIR}/tools
	${CMAKE_SOURCE_DIR}/seriation_gen
	${CMAKE_SOURCE_DIR}/third_party/seqan
	)

#ADD_EXECUTABLE(          t_randadjmat ${SRC} )
#TARGET_LINK_LIBRARIES (  t_randadjmat ${Boost_LIBRARIES} ${NANA_LIBRARIES} hrconfig )
#ADD_TEST(                t_randadjmat ./t_randadjmat --catch_system_errors=no )

IF(SDPA_FOUND)                                                                                                                                   
	ADD_EXECUTABLE(t_sdpa_wrapper t_sdpa_wrapper.cpp)                                                                                  
	TARGET_LINK_LIBRARIES ( t_sdpa_wrapper ${Boost_LIBRARIES} ${NANA_LIBRARIES} hrconfig sdp_wrappers )                                   
	ADD_TEST(t_sdpa_wrapper ./t_sdpa_wrapper --catch_system_errors=no )                                                                
ENDIF(SDPA_FOUND) 
IF(SDPLR_FOUND)                                                                                                                                  
	ADD_EXECUTABLE(t_sdplr_wrapper t_sdplr_wrapper.cpp)                                                                                
	TARGET_LINK_LIBRARIES ( t_sdplr_wrapper ${Boost_LIBRARIES} ${NANA_LIBRARIES} hrconfig sdp_wrappers )                                  
	ADD_TEST(t_sdplr_wrapper ./t_sdplr_wrapper --catch_system_errors=no )                                                              
ENDIF(SDPLR_FOUND)
IF(DSDP_FOUND)                                                                                                                                   
	ADD_EXECUTABLE(t_dsdp_wrapper t_dsdp_wrapper.cpp)                                                                                  
	TARGET_LINK_LIBRARIES ( t_dsdp_wrapper ${Boost_LIBRARIES} ${NANA_LIBRARIES} hrconfig sdp_wrappers )                                   
	ADD_TEST(t_dsdp_wrapper ./t_dsdp_wrapper --catch_system_errors=no )                                                                
ENDIF(DSDP_FOUND)
IF(CSDP_FOUND)                                                                                                                                   
	ADD_EXECUTABLE(t_csdp_wrapper t_csdp_wrapper.cpp)                                                                                  
	TARGET_LINK_LIBRARIES ( t_csdp_wrapper ${Boost_LIBRARIES} ${NANA_LIBRARIES} hrconfig sdp_wrappers )                                   
	ADD_TEST(t_csdp_wrapper ./t_csdp_wrapper --catch_system_errors=no )                                                                
ENDIF(CSDP_FOUND)

ADD_EXECUTABLE(t_sdp_ser_gen t_sdp_seriation_gen.cpp)
TARGET_LINK_LIBRARIES ( t_sdp_ser_gen ${Boost_LIBRARIES} ${NANA_LIBRARIES} ${LAPACK_LIBRARIES} hrconfig hrtools prob_dt seriation_gen sdp_wrappers graph_gen )
ADD_TEST(t_sdp_ser_gen ./t_sdp_ser_gen --catch_system_errors=no )

ADD_EXECUTABLE(t_sdp_ser_prob_gen t_sdp_ser_prob_gen.cpp)                                                                                       
TARGET_LINK_LIBRARIES ( t_sdp_ser_prob_gen ${Boost_LIBRARIES} ${NANA_LIBRARIES} ${LAPACK_LIBRARIES} hrconfig hrtools prob_dt seriation_gen graph_gen )                  
ADD_TEST(t_sdp_ser_prob_gen ./t_sdp_ser_prob_gen) 

ADD_EXECUTABLE(t_config t_config.cpp)                                                                                       
TARGET_LINK_LIBRARIES ( t_config ${Boost_LIBRARIES} ${NANA_LIBRARIES} ${LAPACK_LIBRARIES} hrconfig )                  
ADD_TEST(t_config ./t_config) 

ADD_EXECUTABLE(t_jumbled_adjmat_gen t_jumbled_adjmat_gen.cpp)                                                                                       
TARGET_LINK_LIBRARIES ( t_jumbled_adjmat_gen ${Boost_LIBRARIES} ${NANA_LIBRARIES} ${LAPACK_LIBRARIES} hrconfig graph_gen prob_dt )                  
ADD_TEST(t_jumbled_adjmat_gen ./t_jumbled_adjmat_gen) 

ADD_EXECUTABLE(t_degreesort t_degreesort.cpp)                                                                                       
TARGET_LINK_LIBRARIES ( t_degreesort ${Boost_LIBRARIES} ${NANA_LIBRARIES} ${LAPACK_LIBRARIES} hrconfig graph_gen prob_dt )                  
ADD_TEST(t_degreesort ./t_degreesort) 

ADD_EXECUTABLE(t_origidx t_origidx.cpp)                                                                                       
TARGET_LINK_LIBRARIES ( t_origidx ${Boost_LIBRARIES} ${NANA_LIBRARIES} ${LAPACK_LIBRARIES} hrconfig graph_gen prob_dt )                  
ADD_TEST(t_origidx ./t_origidx) 

ADD_EXECUTABLE(t_seqan t_seqan.cpp)                                                                                       
TARGET_LINK_LIBRARIES ( t_seqan ${Boost_LIBRARIES} ${NANA_LIBRARIES} ${LAPACK_LIBRARIES} hrconfig graph_gen prob_dt )                  
ADD_TEST(t_seqan ./t_seqan) 
