enable_testing()

# Binomial Test
add_executable(BinomialTest discrete/binomial_test.cpp)
target_include_directories(BinomialTest PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(BinomialTest gtest_main diceforge)
add_test(NAME BinomialTest COMMAND BinomialTest)

# Bernoulli Test
add_executable(BernoulliTest discrete/bernoulli_test.cpp)
target_include_directories(BernoulliTest PRIVATE ${PROJECT_SOURCE_DIR}/include)
target_link_libraries(BernoulliTest gtest_main diceforge)
add_test(NAME BernoulliTest COMMAND BernoulliTest)

