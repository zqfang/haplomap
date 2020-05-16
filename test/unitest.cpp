#include "gtest/gtest.h"

int add(int a, int b) {
    return a + b;
}
   
TEST(add_function, Two_Plus_Two_Equals_Four)
{
    EXPECT_EQ(add(2, 2), 4);
}