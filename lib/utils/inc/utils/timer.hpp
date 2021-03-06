#pragma once

#include <chrono>

namespace sight
{
    class Timer
    {
    public:
        Timer()
        {
            Reset();
        }

        void Reset()
        {
            start = std::chrono::high_resolution_clock::now();
        }

        double DurationSecs() const
        {
            const auto now = std::chrono::high_resolution_clock::now();
            const auto dur = std::chrono::duration_cast<std::chrono::duration<double>>(now - start);
            return dur.count();
        }

        double DurationAndReset()
        {
            const auto duration = DurationSecs();
            Reset();
            return duration;
        }
    
    private:
        std::chrono::high_resolution_clock::time_point start;
    };
}