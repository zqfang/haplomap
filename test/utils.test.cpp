//
// Created by Zhuoqing Fang on 6/15/20.
//
#include <string>
#include <string.h>
//#include <filesystem>
#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "utils.h"
#include "vep.h"



 //use "DISABLED_Two_Plus_Two_Equals_Four" to skip test
 //or on cmd line, add --gtest_filter
TEST(SPLIT_TRIM, DISABLED_split_and_trim) {
    std::string alleles("   \rsplit\ta\tb\tc\n\r");
    std::vector<std::string> results = split(alleles, '\t');
    for (auto &res: results)
        std::cout << res << std::endl;
    EXPECT_EQ(results[1], "a");

    //alleles.erase();
    EXPECT_EQ(trim(alleles, " \n\r\t"), "split\ta\tb\tc");
}

TEST(VEP_CWD, getcwd_glob_function) {
    std::string dir = GetCurrentWorkingDir();
    std::cout<<"MyFunction: __file__: "<< dir <<std::endl;
    //std::cout << "CPP17: Current path is " << std::__fs::filesystem::current_path() << '\n';
    //EXPECT_EQ(dir, std::__fs::filesystem::current_path());  

}

TEST(VEP_SUM, DISABLED_vep_fuction) {
    std::string dir = GetCurrentWorkingDir();
    std::cout<<"MyFunction: __file__: "<< dir <<std::endl;
    VarirantEeffectPredictor vep((char*)"../../data/chrX.pass.vep.txt", 
                                  (char*)"../../data/trait.000.txt");
    vep.writeVEPCsq((char*)"../../data/test.csq.txt");
    //vep.writeVEPSummay((char*)"../../data/test.impact.txt");
    // VarirantEeffectPredictor vep((char*)"../../data/chrX.indel.pass.vep.txt", 
    //                             NULL,
    //                             (char*)"../../data/test.vep.txt");


}

TEST(TOLOWER, tolower_function) {
    char *p1 = (char*)"All";
    char *p2 = strdup(p1);
    char*p = p2;
    for ( ; *p; ++p) *p = std::tolower(*p);
    std::cout<<(const char*)p1<<std::endl;
    std::cout<<(const char*)p2<<std::endl;
    // while (name) 
    // {
    //     *name = std::tolower(*name);
    //     name++;
    // }
    // std::cout<<*p<<std::endl;
    free(p2);
}




//// Real class we want to mock
//class MathMock {
//public:
//    virtual ~MathMock() {}
//
//    virtual int getAnswer() {
//        return multiply(3, 3);
//    }
//
//private:
//    virtual int multiply(int a, int b) {
//        return a * b;
//    }
//
//};
//
//// Mock class
//class MockMathMock : public MathMock {
//public:
//    MOCK_METHOD2(multiply, int(int a, int b));
//};
//
//
//using ::testing::Return;
//using ::testing::_;
//
//TEST(MathMockTest, DISABLED_Multiply) {
//MockMathMock MathMock;
//EXPECT_CALL(MathMock, multiply(_, _)).WillOnce(Return(6));
//EXPECT_EQ(MathMock.getAnswer(), 6);
//}

