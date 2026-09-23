set(FORCE_STYLE_TEST_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../force-styles)
set(TEST_INPUT_FOLDER ${CMAKE_CURRENT_SOURCE_DIR}/tests)

if (PKG_LATBOLTZ AND BUILD_MPI)
    #test for lattice boltzman multicomponent
    add_executable(test_fix_lb_multicomponent test_fix_lb_multicomponent.cpp)

    if(NOT BUILD_SHARED_LIBS)
        target_compile_definitions(test_fix_lb_multicomponent PRIVATE USING_STATIC_LIBS=1)
    endif()

    target_include_directories(test_fix_lb_multicomponent PRIVATE ${FORCE_STYLE_TEST_DIR}
                                ${LAMMPS_SOURCE_DIR}/LATBOLTZ)

    target_link_libraries(test_fix_lb_multicomponent PRIVATE style_tests)

    file(GLOB LB_FLUID_TESTS LIST_DIRECTORIES false ${TEST_INPUT_FOLDER}/test_lb_multicomponent*.yaml)
    foreach(TEST ${LB_FLUID_TESTS})
    add_mpi_test(NAME FixLBMulticomponent:plain NUM_PROCS 2 COMMAND
                 $<TARGET_FILE:test_fix_lb_multicomponent> ${TEST}
                 --gtest_filter=FixLBMulticomponent.plain)
    add_mpi_test(NAME FixLBMulticomponent:moments NUM_PROCS 1 COMMAND
                 $<TARGET_FILE:test_fix_lb_multicomponent> ${TEST}
                 --gtest_filter=LBMomentReconstruction.*)
    set_tests_properties(FixLBMulticomponent:plain PROPERTIES PROCESSORS 2 TIMEOUT 120)
    set_tests_properties(FixLBMulticomponent:moments PROPERTIES PROCESSORS 1 TIMEOUT 30)
    endforeach()
endif()
#set_tests_properties(FixLBMulticomponent:plain PROPERTIES
#                    PASS_REGULAR_EXPRESSION "Total wall time")
