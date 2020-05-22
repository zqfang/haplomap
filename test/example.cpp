#include "gtest/gtest.h"
#include "gmock/gmock.h"

// Real class we want to mock
class MathMock {
public:
    virtual ~MathMock() {}

    virtual int getAnswer() {
        return multiply(3, 3);
    }

private:
    virtual int multiply(int a, int b) {
        return a * b;
    }

};

// Mock class
class MockMathMock : public MathMock {
public:
    MOCK_METHOD2(multiply, int(int a, int b));
};


using ::testing::Return;
using ::testing::_;

TEST(MathMockTest, Multiply) {
    MockMathMock MathMock;
    EXPECT_CALL(MathMock, multiply(_, _)).WillOnce(Return(6));
    EXPECT_EQ(MathMock.getAnswer(), 6);
}
