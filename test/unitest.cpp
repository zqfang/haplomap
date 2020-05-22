#include "gtest/gtest.h"

int add(int a, int b) {
    return a + b;
}
// use "DISABLED_Two_Plus_Two_Equals_Four" to skip test   
TEST(add_function, Two_Plus_Two_Equals_Four) 
{
    EXPECT_EQ(add(2, 2), 4);
}

// or on cmd line, add --gtest_filter