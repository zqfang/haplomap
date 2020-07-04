//
// Created by Zhuoqing Fang on 6/15/20.
//
#include <string>
//#include <filesystem>
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "utils.h"
#include "glob.h"



 //use "DISABLED_Two_Plus_Two_Equals_Four" to skip test
 //or on cmd line, add --gtest_filter
TEST(SPLIT_TRIM, split_and_trim) {
    std::string alleles("   \rsplit\ta\tb\tc\n\r");
    std::vector<std::string> results = split(alleles, '\t');
    for (auto &res: results)
        std::cout << res << std::endl;
    EXPECT_EQ(results[1], "a");

    //alleles.erase();
    EXPECT_EQ(trim(alleles, " \n\r\t"), "split\ta\tb\tc");
}

TEST(GET_CWD_GLOB, getcwd_glob_function) {
    std::string dir = GetCurrentWorkingDir();
    std::cout<<"MyFunction: __file__: "<< dir <<std::endl;
    //std::cout << "CPP17: Current path is " << std::__fs::filesystem::current_path() << '\n';
    //EXPECT_EQ(dir, std::__fs::filesystem::current_path());

    glob::glob globs(dir);
    while (globs) {
        std::cout << globs.current_match() << std::endl;
        globs.next();
    }

}




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

TEST(MathMockTest, DISABLED_Multiply) {
MockMathMock MathMock;
EXPECT_CALL(MathMock, multiply(_, _)).WillOnce(Return(6));
EXPECT_EQ(MathMock.getAnswer(), 6);
}

