add_test(
  NAME lb_multicomponent_smoke
  COMMAND $<TARGET_FILE:lmp> -log none -echo none -in
          ${CMAKE_CURRENT_SOURCE_DIR}/in.lb_multicomponent_smoke
)

set_tests_properties(lb_multicomponent_smoke PROPERTIES
  PASS_REGULAR_EXPRESSION "Total wall time"
)
