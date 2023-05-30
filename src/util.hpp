#ifndef UTIL_H
#define UTIL_H

#include <sys/time.h>
#include <sys/resource.h>
#include <assert.h>

#include <random>
#include <chrono>

typedef long long vertex_id_t;
class RandNumGenerator
{
public:
    virtual vertex_id_t gen(vertex_id_t upper_bound) = 0;
    virtual float gen_float(float upper_bound) = 0;
    virtual ~RandNumGenerator() {}
};

class StdRandNumGenerator : public RandNumGenerator
{
    std::random_device *rd;
    std::mt19937 *mt;
public:
    StdRandNumGenerator()
    {
        rd = new std::random_device();
        mt = new std::mt19937((*rd)());
    }
    ~StdRandNumGenerator()
    {
        delete mt;
        delete rd;
    }
    vertex_id_t gen(vertex_id_t upper_bound)
    {
        std::uniform_int_distribution<vertex_id_t> dis(0, upper_bound - 1);
        return dis(*mt);
    }
    float gen_float(float upper_bound)
    {
        std::uniform_real_distribution<float> dis(0.0, upper_bound);
        return dis(*mt);
    }
};

//Timer is used for performance profiling
class Timer
{
    std::chrono::time_point<std::chrono::system_clock> _start = std::chrono::system_clock::now();
public:
    void restart()
    {
        _start = std::chrono::system_clock::now();
    }
    double duration()
    {
        std::chrono::duration<double> diff = std::chrono::system_clock::now() - _start;
        return diff.count();
    }
    static double current_time()
    {
        std::chrono::duration<double> val = std::chrono::system_clock::now().time_since_epoch();
        return val.count();
    }
};
#endif
