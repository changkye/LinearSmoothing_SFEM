#pragma once

#include <chrono>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

#if defined(__unix__) || defined(__APPLE__)
#include <sys/resource.h>
#endif

namespace benchmark
{

struct Record
{
    std::string timestamp_utc;
    std::string problem_type;
    std::string method;
    int num_els_x = 0;
    int num_els_y = 0;
    int nstep = 0;
    int maxiter = 0;
    std::string solver;
    bool compute_exact_solution = true;
    bool write_vtu = true;
    bool write_postprocess = true;
    double relative_error = 0.0;
    double strain_energy = 0.0;
    double solve_wall_time_sec = 0.0;
    double total_wall_time_sec = 0.0;
    double peak_rss_mb = 0.0;
};

inline std::string utc_timestamp_now()
{
    const std::time_t now = std::time(nullptr);
    std::tm tm{};
#if defined(_WIN32)
    gmtime_s(&tm, &now);
#else
    gmtime_r(&now, &tm);
#endif
    std::ostringstream out;
    out << std::put_time(&tm, "%Y-%m-%dT%H:%M:%SZ");
    return out.str();
}

inline double peak_rss_mb()
{
#if defined(__unix__) || defined(__APPLE__)
    rusage usage{};
    if (getrusage(RUSAGE_SELF, &usage) != 0)
    {
        return 0.0;
    }
#if defined(__APPLE__)
    return static_cast<double>(usage.ru_maxrss) / (1024.0 * 1024.0);
#else
    return static_cast<double>(usage.ru_maxrss) / 1024.0;
#endif
#else
    return 0.0;
#endif
}

inline void append_csv(const std::filesystem::path &file_path, const Record &record)
{
    std::filesystem::create_directories(file_path.parent_path());
    const bool write_header = !std::filesystem::exists(file_path) || std::filesystem::file_size(file_path) == 0;

    std::ofstream out(file_path, std::ios::app);
    if (!out)
    {
        throw std::runtime_error("Failed to open benchmark log: " + file_path.string());
    }

    out << std::setprecision(16);
    if (write_header)
    {
        out << "timestamp_utc,problem_type,method,num_els_x,num_els_y,nstep,maxiter,solver,compute_exact_solution,"
               "write_vtu,write_postprocess,relative_error,strain_energy,solve_wall_time_sec,total_wall_time_sec,peak_rss_mb\n";
    }

    out << record.timestamp_utc << ','
        << record.problem_type << ','
        << record.method << ','
        << record.num_els_x << ','
        << record.num_els_y << ','
        << record.nstep << ','
        << record.maxiter << ','
        << record.solver << ','
        << (record.compute_exact_solution ? 1 : 0) << ','
        << (record.write_vtu ? 1 : 0) << ','
        << (record.write_postprocess ? 1 : 0) << ','
        << record.relative_error << ','
        << record.strain_energy << ','
        << record.solve_wall_time_sec << ','
        << record.total_wall_time_sec << ','
        << record.peak_rss_mb << '\n';
}

template <typename Clock = std::chrono::steady_clock>
inline double seconds_between(const typename Clock::time_point &start,
                              const typename Clock::time_point &end)
{
    return std::chrono::duration<double>(end - start).count();
}

} // namespace benchmark
