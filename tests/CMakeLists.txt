add_executable(unit_tests binomial_test.cpp)
target_include_directories(unit_tests PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(unit_tests
    gtest_main
    diceforge  
)
include(GoogleTest)
gtest_discover_tests(unit_tests)
