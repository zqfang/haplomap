//
// Created by Zhuoqing Fang on 6/15/20.
//

#ifndef UTILS_H
#define UTILS_H


#include <stdio.h>  /* defines FILENAME_MAX */
#ifdef WIN32
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

#include <memory>
#include <vector>
#include <string>
#include <iostream>
#include <sstream>


std::string GetCurrentWorkingDir( void );
std::vector<std::string> split(const std::string& s, char delimiter);

// modern version of string parsing
inline std::string& ltrim(std::string& str, const std::string& delimiter = "\t\n\v\f\r ")
{
    str.erase(0, str.find_first_not_of(delimiter));
    return str;
}

inline std::string& rtrim(std::string& str, const std::string& delimiter = "\t\n\v\f\r ")
{
    str.erase(str.find_last_not_of(delimiter) + 1);
    return str;
}
inline std::string& trim(std::string& str, const std::string& delimiter = "\t\n\v\f\r ")
{
    return ltrim(rtrim(str, delimiter), delimiter);
}


class Singleton
{
private:
    /* Here will be the instance stored. */
    static Singleton* instance;
    /* Private constructor to prevent instancing. */
    Singleton() {};
public:
    /* Static access method. */
    static Singleton* getInstance();
};



//#include <filesystem>
//namespace fs = std::filesystem;
//int main()
//{
//    std::cout << "Current path is " << fs::current_path() << '\n';
//}


// Get Logger
class MyLogger
{
public:
    static std::shared_ptr<MyLogger> getInstance();
    //log4cplus::Logger m_rootLog;
    ~MyLogger();

private:
    MyLogger();
    static std::shared_ptr<MyLogger> m_logger;
};


#endif // UTILS_H
