#include "algorithm/timer.hpp"


std::chrono::high_resolution_clock::time_point t_start, t_end;
std::string original_file_name;

double subtractTimes(std::chrono::high_resolution_clock::time_point& t_end, std::chrono::high_resolution_clock::time_point& t_start)
{
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(t_end - t_start);
    return elapsed.count() / 1000000.0;
}

double subtractTimes(std::chrono::high_resolution_clock::time_point& t_start)
{
    auto t_end = std::chrono::high_resolution_clock::now();
    return subtractTimes(t_end, t_start);
}
