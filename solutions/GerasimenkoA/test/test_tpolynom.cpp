#include "tpolynom.h"

#include <gtest.h>

using namespace std;

TEST(TPolynom, can_set_polynom) {
    TPolynom p;
    ASSERT_NO_THROW(p.SetPolynom("3x^2y^4 + 4.5yz^2"));
}

TEST(TPolynom, can_set_empty_polynom)
{
    TPolynom p;
    ASSERT_NO_THROW(p.SetPolynom(""));
}

TEST(TPolynom, cant_set_polynom_err) {
    TPolynom p;
    ASSERT_ANY_THROW(p.SetPolynom("x^"));
}

TEST(TPolynom, cant_set_polynom_just_symbol) {
    TPolynom p;
    ASSERT_ANY_THROW(p.SetPolynom("+"));
}

TEST(TPolynom, is_correct_returns_true) {
    TPolynom p;
    p.SetPolynom("2x^2 + 3x");
    EXPECT_TRUE(p.IsCorrect());
}

TEST(TPolynom, is_correct_returns_true_for_zero) {
    TPolynom p;
    EXPECT_TRUE(p.IsCorrect());
}

TEST(TPolynom, plus_operator) {
    TPolynom p1, p2;
    p1.SetPolynom("3x^2y^4 + 4.5yz^2");
    p2.SetPolynom("3x^2y^4 + 4.5yz^2");
    TPolynom p3 = p1 + p2;
    double val = p3.Calculate(1.0, 1.0, 1.0);
    double expected = 6 * pow(1.0, 2) * pow(1.0, 4) + 9 * 1.0 * pow(1.0, 2);
    EXPECT_EQ(val, expected);
}

TEST(TPolynom, minus_operator) {
    TPolynom p1, p2;
    p1.SetPolynom("2x^2 + 2y");
    p2.SetPolynom("3x^2 + 3y");
    TPolynom p3 = p1 - p2;
    double val = p3.Calculate(1.0, 2.0, 0.0);
    EXPECT_EQ(val, -3);
}

TEST(TPolynom, multiply_polynom_by_polynom) {
    TPolynom p1, p2;
    p1.SetPolynom("x + y");
    p2.SetPolynom("x + y");
    TPolynom p3 = p1 * p2;
    double val = p3.Calculate(1.0, 1.0, 0.0);
    EXPECT_EQ(val, 4);
}

TEST(TPolynom, multiply_polynom_by_double) {
    TPolynom p1;
    p1.SetPolynom("2x + 3");
    TPolynom p2 = p1 * 2.0;
    double val = p2.Calculate(1.0, 0.0, 0.0);
    EXPECT_EQ(val, 10);
}

TEST(TPolynom, multiply_double_by_polynom) {
    TPolynom p1;
    p1.SetPolynom("2x + 3");
    TPolynom p2 = 2.0 * p1;
    double val = p2.Calculate(1.0, 0.0, 0.0);
    EXPECT_EQ(val, 10);
}

TEST(TPolynom, equal_polynoms_are_equal) {
    TPolynom p1, p2;
    p1.SetPolynom("3x^2y^4 + 4.5yz^2");
    p2.SetPolynom("3x^2y^4 + 4.5yz^2");
    EXPECT_TRUE(p1 == p2);
    EXPECT_FALSE(p1 != p2);
}

TEST(TPolynom, different_polynoms_are_not_equal) {
    TPolynom p1, p2;
    p1.SetPolynom("3x^2y^4 + 4.5yz^2");
    p2.SetPolynom("2x^2y^4 + 4.5yz^2");
    EXPECT_FALSE(p1 == p2);
    EXPECT_TRUE(p1 != p2);
}

TEST(TPolynom, empty_polynoms_are_equal) {
    TPolynom p1, p2;
    EXPECT_TRUE(p1 == p2);
    EXPECT_FALSE(p1 != p2);
}

TEST(TPolynom, zero_polynoms_are_equal) {
    TPolynom p1, p2;
    p1.SetPolynom("");
    p2.SetPolynom("0");
    EXPECT_TRUE(p1 == p2);
}

TEST(TPolynom, can_add_empty_monom) {
    TPolynom p;
    ASSERT_NO_THROW(p.Add(""));
}

TEST(TPolynom, can_delete_monom_by_position) {
    TPolynom p;
    p.SetPolynom("x^2 + y + z");
    p.Delete(1); 
    double val = p.Calculate(1.0, 1.0, 1.0);
    EXPECT_EQ(val, 2);
}

TEST(TPolynom, cant_delete_out_of_range_position) 
{
    TPolynom p;
    p.SetPolynom("x^2");
    ASSERT_ANY_THROW(p.Delete(5));
}