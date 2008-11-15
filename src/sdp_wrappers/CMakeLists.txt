SET(CMAKE_MODULE_PATH  ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/../CMakeModules )

# const-char* -> char* conversion, needed for Nana
IF( CMAKE_COMPILER_IS_GNUCXX )
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-write-strings -Wall" )
ENDIF()

include_directories( 
	. 
	${CMAKE_SOURCE_DIR}/configuration
	${CMAKE_SOURCE_DIR}/seriation_gen
	)

SET(SRC ${SRC}  )
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
SET(SRC ${SRC} sdp_wrapper_cfg.cpp sdp_wrapper.cpp sdp_wrapper.hpp )

add_library( sdp_wrappers SHARED ${SRC} )
include_directories( ${CMAKE_CURRENT_BINARY_DIR})
add_precompiled_header( sdp_wrappers ${CMAKE_CURRENT_SOURCE_DIR}/sdp_wrappers_prec.hpp )
