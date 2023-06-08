#ifndef TIME_H
#define TIME_H

#include <string>
#include <chrono>

using std::string;
extern std::string original_file_name;

extern std::chrono::high_resolution_clock::time_point t_start;

double subtractTimes(std::chrono::high_resolution_clock::time_point& t_end, std::chrono::high_resolution_clock::time_point& t_start);
double subtractTimes(std::chrono::high_resolution_clock::time_point& t_start);


#endif // TIME_H
