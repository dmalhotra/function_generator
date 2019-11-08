#include <chrono>
#include <cmath>
#include <functional>
#include <iostream>
#include <random>
#include <vector>

#include "function_generator.hpp"
#include <omp.h>

template class FunctionGenerator<8, 4096>;

int main(int argc, char *argv[]) {
    using namespace std::chrono;
    high_resolution_clock::time_point start = high_resolution_clock::now();

    double range[2] = {1e-10, 1000};
    std::function<double(double)> inFunc = log;
    FunctionGenerator<7, 4096> myLog(inFunc, range[0], range[1], 1E-12);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(range[0], range[1]);
    size_t n_el = 1000000;
    int n_loops = 100;
    std::vector<double> x(n_el);
    for (size_t i = 0; i < n_el; ++i)
        x[i] = dis(gen);
    high_resolution_clock::time_point finish = high_resolution_clock::now();
    duration<double> time_span =
        duration_cast<duration<double>>(finish - start);

    std::cout << "RNG generation took " << time_span.count() << " seconds.\n";

    std::cout << "First 10 test value deltas\n";
    for (auto i = 0; i < 10; ++i) {
        std::cout << "\t" << myLog(x[i]) - log(x[i]) << std::endl;
    }

    {
        std::vector<double> res(x.size());
        high_resolution_clock::time_point start = high_resolution_clock::now();

        double sum = 0.0;
        for (int irun = 0; irun < n_loops; ++irun) {
            for (size_t i = 0; i < n_el; ++i) {
                res[i] = log(x[i]);
                sum += res[i];
            }
        }
        std::cout << sum << std::endl;
        high_resolution_clock::time_point finish = high_resolution_clock::now();
        duration<double> time_span =
            duration_cast<duration<double>>(finish - start);

        std::cout << "Efficiency: " << n_el * n_loops / time_span.count() / 1e6
                  << " million evaluations / sec\n";
        std::cout << "System 'log' took " << time_span.count() << " seconds.\n";
    }

    {
        std::vector<double> res(x.size());
        high_resolution_clock::time_point start = high_resolution_clock::now();

        double sum = 0.0;
        for (int irun = 0; irun < n_loops; ++irun) {
            for (size_t i = 0; i < n_el; ++i) {
                res[i] = myLog(x[i]);
                sum += res[i];
            }
        }
        std::cout << sum << std::endl;
        high_resolution_clock::time_point finish = high_resolution_clock::now();
        duration<double> time_span =
            duration_cast<duration<double>>(finish - start);

        std::cout << "Chebyshev 'log' took " << time_span.count()
                  << " seconds.\n";
        std::cout << "Efficiency: " << n_el * n_loops / time_span.count() / 1e6
                  << " million evaluations / sec\n";
    }

    return 0;
}
