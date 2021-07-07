#ifndef LOGGER_H_INLUDED
#define LOGGER_H_INLUDED

#include <iostream>

#include <string>

typedef enum {
    LOG_DEBUG_VERBOSE = 0,
    LOG_DEBUG_POINT,
    LOG_DEBUG,
    LOG_INFO,
    LOG_WARN,
    LOG_ERROR,
    LOG_FATAL
} LogLevel;

class Logger;

class Log {
    public:
        Log();
        ~Log();

        static void setLogThreshold(LogLevel level);
        static LogLevel getLogThreshold();

    private:
        static LogLevel threshold;
};

class Logger {
    public:
//	Logger()  {};
        Logger(const std::string& _name) : name(_name) {};
        std::ostream& operator<<(LogLevel level);
        void newLine(LogLevel level);

    private:
        LogLevel level;
        std::string name;
};

class NoopStream : public std::ostream {
    public:
        template<class T>
        NoopStream& operator<<(const T &arg) { return *this; };

        virtual ~NoopStream() {};
};

#endif
