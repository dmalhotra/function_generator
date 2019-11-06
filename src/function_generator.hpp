#ifndef FUNCTIONGENERATOR_H
#define FUNCTIONGENERATOR_H
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/LU>
#include <functional>
#include <iostream>
#ifdef PYTHON_MODULE
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#endif
#include <vector>
typedef Eigen::VectorXd vecxd;

#ifdef PYTHON_MODULE
namespace py = pybind11;
#endif

class FunctionGenerator {

  public:
    FunctionGenerator(std::function<double(double)> &f, double a, double b,
                      double tol = 1E-10, int n = 12, double mw = 1E-15);

#ifdef PYTHON_MODULE
    FunctionGenerator(py::function fpy, double a, double b, double tol = 1e-10,
                      int n = 12, double mw = 1e-15);

    // py::array_t<double> arr_call_old(py::array_t<double> x) {
    //     py::buffer_info xbuf = x.request();
    //     auto res = py::array_t<double>(xbuf.size);
    //     py::buffer_info rbuf = res.request();
    //     double *xptr = (double *)xbuf.ptr;
    //     double *rptr = (double *)rbuf.ptr;

    //     for (auto i = 0; i < xbuf.size; ++i) {
    //         rptr[i] = (*this)(xptr[i]);
    //     }

    //     return res;
    // };

    py::array_t<double> arr_call(py::array_t<double> x);
#endif

    double operator()(double x);

  private:
    std::function<double(double)> f_;
    const double a_;
    const double b_;
    const double tol_;
    const int n_;
    const double mw_;

    Eigen::MatrixXd V_;
    Eigen::PartialPivLU<Eigen::MatrixXd> VLU_;

    std::vector<double> lbs_;
    std::vector<double> ubs_;
    std::vector<double> coeffs_;
    std::vector<std::pair<uint16_t, uint16_t>> bounds_table_;

    void init();

    double chbevl(double x, double *c);
    
    void init_vandermonde();

    vecxd get_chebyshev_nodes(double lb, double ub, int order);

    void fit(double a, double b);
    void init_lookup();

    int bisect_bracketed(double x, int n1, int n2) {
        while (n2 - n1 > 1) {
            const int m = n1 + (n2 - n1) / 2;
            if (x < lbs_[m])
                n2 = m;
            else
                n1 = m;
        }

        return n1;
    };

    int bisect(double x) { return bisect_bracketed(x, 0, lbs_.size()); }

    int bisect_cache(double x) {
        constexpr int half_width_cache = 4;
        static int n1 = half_width_cache;
        n1 -= half_width_cache;
        n1 = std::max(n1, 0);
        int n2 = n1 + 2 * half_width_cache;

        while (lbs_[n1] > x)
            n1 /= 2;
        while (n2 < (int)lbs_.size() && lbs_[n2] < x)
            n2 *= 2;
        n2 = std::min(n2, (int)lbs_.size());

        n1 = bisect_bracketed(x, n1, n2);
        return n1;
    };

    int bisect_lookup(double x) {
        int table_index = x / (b_ - a_) * bounds_table_.size();
        auto bisect_bounds = bounds_table_[table_index];
        return bisect_bracketed(x, bisect_bounds.first, bisect_bounds.second);
    };

    double standard_error(vecxd &coeffs) {
        double maxcoeff = 0.0;
        for (auto i = n_ - 2; i < n_; ++i)
            maxcoeff = std::max(fabs(coeffs[i]), maxcoeff);

        return maxcoeff / std::max(1.0, fabs(coeffs[0]));
    };
};

#endif
