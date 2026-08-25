#test for lattice boltzman multicomponent
add_executable(test_fix_lb_multicomponent test_fix_lb_multicomponent.cpp)

if(NOT BUILD_SHARED_LIBS)
    target_compile_definitions(test_fix_lb_multicomponent PRIVATE USING_STATIC_LIBS=1)
endif()

target_link_libraries(test_fix_lb_multicomponent PRIVATE lammps GTest::GTest)

add_test(NAME FixLBMulticomponent:plain COMMAND test_fix_lb_multicomponent)

set_tests_properties(FixLBMulticomponent:plain PROPERTIES
                    PASS_REGULAR_EXPRESSION "Total wall time")
