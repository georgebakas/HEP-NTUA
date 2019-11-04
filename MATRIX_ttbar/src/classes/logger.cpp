#include "header/logger.h"

LogLevel Log::threshold = LOG_FATAL;

Log::Log() {}
Log::~Log() {}

void Log::setLogThreshold(LogLevel level) {
    threshold = level;
}

LogLevel Log::getLogThreshold() {
    return threshold;
}

std::ostream& Logger::operator<<(LogLevel level) {
    if(level >= Log::getLogThreshold()) {
//        std::clog << name << ": ";
//        return std::clog;
      std::cout << name << ": ";
      return std::cout;
    } else {
        static NoopStream s;
        return s;
    }
}

void Logger::newLine(LogLevel level) {
  if(level >= Log::getLogThreshold()) {
    std::cout << std::endl;
  }
}
