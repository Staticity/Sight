#include <catch2/catch.hpp>

#include <linear/vec.hpp>

TEST_CASE("Test vector creation")
{
    const sight::Vec<float, 5> x = {1.2f, 3.4f, 5.6f, 7.8f, 9};

    REQUIRE(x(0) == 1.2f);
    REQUIRE(x(1) == 3.4f);
    REQUIRE(x(2) == 5.6f);
    REQUIRE(x(3) == 7.8f);
    REQUIRE(x(4) == 9.0f);
}

TEST_CASE("Test vector addition")
{
    const sight::Vec2<int> a = {1, 2};
    const sight::Vec2<int> b = {5, 3};
    const auto c = a + b;

    REQUIRE(c(0) == 6);
    REQUIRE(c(1) == 5);
}

TEST_CASE("Test vector subtraction")
{
    const sight::Vec2<int> a = {17, 0};
    const sight::Vec2<int> b = {4, 2};
    const auto c = a - b;

    REQUIRE(c(0) == 13);
    REQUIRE(c(1) == -2);
}
