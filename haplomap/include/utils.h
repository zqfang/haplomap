//
// Created by Zhuoqing Fang on 6/15/20.
//

#ifndef PLEIADES_UTILS_H
#define PLEIADES_UTILS_H


#include <stdio.h>  /* defines FILENAME_MAX */
#ifdef WIN
#include <direct.h>
#define GetCurrentDir _getcwd
#else
#include <unistd.h>
#define GetCurrentDir getcwd
#endif

#include <vector>
#include <string>
#include <iostream>
#include <sstream>

std::string GetCurrentWorkingDir( void );
std::vector<std::string> split(const std::string& s, char delimiter);
std::string trim(const std::string& str, const std::string delimiter);


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


inline int toInt(const std::string &s) {
    return std::stoi(s);
}

inline int toFloat(const std::string &s) {
    return std::stof(s);
}

#endif //PLEIADES_UTILS_H
