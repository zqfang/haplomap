//
// Created by Zhuoqing Fang on 6/15/20.
//

#include "utils.h"

std::string GetCurrentWorkingDir( void )
{
    char buff[FILENAME_MAX];
    GetCurrentDir(buff, FILENAME_MAX);
    std::string current_working_dir(buff);
    return current_working_dir;
}

std::vector<std::string> split(const std::string& s, char delimiter)
{
    std::vector<std::string> tokens;
    std::string token;
    std::istringstream tokenStream(s);
    while (std::getline(tokenStream, token, delimiter))
    {
        tokens.push_back(token);
    }
    return tokens;
}

std::string trim(const std::string& str, const std::string delimiter = " \n\r\t")
{
//    std::string s;
//    s.erase(s.find_last_not_of(" \n\r\t")+1);
    size_t first = str.find_first_not_of(delimiter);
    if (std::string::npos == first)
    {
        return str;
    }
    size_t last = str.find_last_not_of(delimiter);
    return str.substr(first, (last - first + 1));
}

/* NULL, because instance will be initialized on demand. */
Singleton* Singleton::instance = 0;
/* Static access method. */
Singleton* Singleton::getInstance()
{
    if (instance == 0)
        instance = new Singleton();
    return instance;
}


// define MyLogger
std::shared_ptr<MyLogger> MyLogger::m_logger = nullptr;

MyLogger::MyLogger() {}

MyLogger::~MyLogger()
{
    std::cout << "delete MyLogger" << std::endl;
}

std::shared_ptr<MyLogger> MyLogger::getInstance()
{
//    if (m_logger == nullptr)
//    {
//        std::cout << "construct MyLogger" << std::endl;
//        log4cplus::initialize();
//
//        log4cplus::PropertyConfigurator::doConfigure(LOG4CPLUS_TEXT(MY_LOG_FILE_PATH));
//        m_logger.reset(new MyLogger);
//
//        m_logger->m_rootLog = log4cplus::Logger::getRoot();
//    }

    return m_logger;
}



