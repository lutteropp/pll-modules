
file(GLOB NETWORK_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.c ${CMAKE_CURRENT_SOURCE_DIR}/../pllmod_common.c)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${PLLMOD_CFLAGS}")

add_pllmodules_lib(pllmodnetwork "${NETWORK_SOURCES}")

target_include_directories(pllmodnetwork_obj PRIVATE ${PLLMOD_INCLUDE_PATH}
  ${CMAKE_CURRENT_SOURCE_DIR}/../tree
  )

