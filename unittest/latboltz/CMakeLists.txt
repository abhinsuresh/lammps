set(FORCE_STYLE_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../force-styles)

#test for lattice boltzman multicomponent
add_executable(test_fix_lb_multicomponent 
               test_fix_lb_multicomponent.cpp
               test_lb_main.cpp
               ${FORCE_STYLE_TEST_DIR}/yaml_writer.cpp
               ${FORCE_STYLE_TEST_DIR}/test_config_reader.cpp)

if(NOT BUILD_SHARED_LIBS)
    target_compile_definitions(test_fix_lb_multicomponent PRIVATE USING_STATIC_LIBS=1)
endif()

target_include_directories(test_fix_lb_multicomponent PRIVATE
                           ${CMAKE_CURRENT_SOURCE_DIR}
                           ${FORCE_STYLE_TEST_DIR}
                           ${LAMMPS_SOURCE_DIR})

target_link_libraries(test_fix_lb_multicomponent PRIVATE lammps 
                      GTest::GTest Yaml::Yaml)


add_test(NAME FixLBMulticomponent:plain COMMAND test_fix_lb_multicomponent
            ${CMAKE_CURRENT_SOURCE_DIR}/tests/fix_lb_multicomponent.yaml)

#set_tests_properties(FixLBMulticomponent:plain PROPERTIES
#                    PASS_REGULAR_EXPRESSION "Total wall time")
