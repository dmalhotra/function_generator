#include <chrono>
#include <cmath>
#include <complex>
#include <fstream>
#include <functional>
#include <iostream>
#include <random>
#include <vector>

#include "function_generator.hpp"
#include <omp.h>

#include <gsl/gsl_sf.h>

typedef std::chrono::high_resolution_clock clk;
typedef std::chrono::duration<double> duration;

template <typename T> struct func_t {
    std::string name;
    double low;
    double high;
    std::function<T(double)> f;
};

template <int n_, int table_size_, typename T>
void timeit(func_t<T> func, std::vector<double> &xraw) {
    FunctionGenerator<n_, table_size_, T> f(func.f, func.low, func.high, 1E-12);
    const double range = func.high - func.low;

    std::vector<double> x = xraw;
    for (auto &xel : x)
        xel = xel * range + func.low;

    std::cout << "First 5 test value deltas\n";
    for (auto i = 0; i < 5; ++i) {
        T yExp = f(x[i]);
        T yLib = func.f(x[i]);
        std::cout << "\t" << std::abs((yExp - yLib) / std::abs(yLib))
                  << std::endl;
    }

    int n_loops = 10;
    std::vector<T> res(x.size());

    {
        clk::time_point start = clk::now();

        for (int irun = 0; irun < n_loops; ++irun) {
            for (size_t i = 0; i < x.size(); ++i) {
                res[i] = func.f(x[i]);
            }
        }

        clk::time_point finish = clk::now();
        duration time_span =
            std::chrono::duration_cast<duration>(finish - start);

        std::cout << "System Efficiency: "
                  << x.size() * n_loops / time_span.count() / 1e6
                  << " million evaluations / sec\n";

        std::ofstream f_out(func.name + "_system.dat");
        for (size_t i = 0; i < 1000; ++i) {
            double a = func.low + i * range / 1000;
            f_out << a << "\t" << func.f(a) << std::endl;
        }
    }

    {
        clk::time_point start = clk::now();

        for (int irun = 0; irun < n_loops; ++irun) {
            for (size_t i = 0; i < x.size(); ++i) {
                res[i] = f(x[i]);
            }
        }

        clk::time_point finish = clk::now();
        duration time_span =
            std::chrono::duration_cast<duration>(finish - start);

        std::cout << "FunctionGenerator Efficiency: "
                  << x.size() * n_loops / time_span.count() / 1e6
                  << " million evaluations / sec\n";
        std::ofstream f_out(func.name + "_expansion.dat");
        for (size_t i = 0; i < 1000; ++i) {
            double a = func.low + i * range / 1000;
            f_out << a << "\t" << f(a) << std::endl;
        }
    }
}

int main(int argc, char *argv[]) {
    size_t n_el = 1000000;

    clk::time_point start = clk::now();

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(0, 1);
    std::vector<double> x(n_el);
    for (size_t i = 0; i < n_el; ++i)
        x[i] = dis(gen);
    clk::time_point finish = clk::now();
    duration time_span = std::chrono::duration_cast<duration>(finish - start);

    std::cout << "RNG generation took " << time_span.count() << " seconds.\n";

    std::vector<func_t<double>> double_funcs{
        func_t<double>{"log", 1e-15, 1000, [](double x) { return log(x); }},
        func_t<double>{"gsl_sf_bessel_J0", 0, 100, gsl_sf_bessel_J0},
        func_t<double>{"gsl_sf_bessel_I0", 0, 100, gsl_sf_bessel_I0},
        func_t<double>{"std::erfc", -2, 2,
                       [](double x) { return std::erfc(x); }},
        func_t<double>{"std::erf", -2, 2, [](double x) { return std::erf(x); }},
        func_t<double>{
            "gsl_sf_airy_Ai", -20, 5,
            [](double x) { return gsl_sf_airy_Ai(x, GSL_PREC_DOUBLE); }},
        func_t<double>{"sin", 0, 2 * M_PI, [](double x) { return sin(x); }}};

    for (auto func : double_funcs) {
        std::cout << std::endl << func.name << std::endl;
        timeit<8, 4096, double>(func, x);
    }

    using std::complex;
    std::vector<func_t<complex<double>>> complex_funcs{
        func_t<complex<double>>{
            "log", 1e-15, 1000,
            [](complex<double> x) { return log(x * complex<double>(0, 1)); }},
        func_t<complex<double>>{"sin", 0, 2 * M_PI, [](complex<double> x) {
                                    return sin(x * complex<double>(0, 1));
                                }}};

    for (auto func : complex_funcs) {
        std::cout << std::endl << func.name << std::endl;
        timeit<8, 4096, complex<double>>(func, x);
    }
    return 0;
}
