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

//! Namespace for error calculation helper functions.
namespace FGError {
//! Model used to calculate error in approximation the input function.
enum ErrorModel : uint16_t { standard = 0, relative = 1 };

inline double standard_error(vecxd &coeffs) {
    double maxcoeff = 0.0;
    int n = coeffs.size();
    for (auto i = n - 2; i < n; ++i)
        maxcoeff = std::max(fabs(coeffs[i]), maxcoeff);

    return maxcoeff / std::max(1.0, fabs(coeffs[0]));
};

inline double relative_error(vecxd &coeffs) {
    double maxcoeff = 0.0;
    int n = coeffs.size();
    for (auto i = n - 2; i < n; ++i)
        maxcoeff = std::max(fabs(coeffs[i]), maxcoeff);

    return maxcoeff / fabs(coeffs[0]);
}

inline double calc_error(vecxd &coeffs, ErrorModel error_model) {
    switch (error_model) {
    case relative:
        return FGError::relative_error(coeffs);
    case standard:
        return FGError::standard_error(coeffs);
    default:
        return standard_error(coeffs);
    }
}

} // namespace FGError

//! Class to approximate functions on a given range by sub-divisions of
//! Chebyshev polynomial expansions.
/*!

This is a template "header-only" class which aims to be virtually
indistinguishable, but often orders of magnitude faster, than special functions
on a given range. This goal is accomplished by recursively dividing the input
function into sub-regions which are fit to Chebyshev polynomials of a given
order until the error is within some given tolerance.

To aid in optimization for C++, this class is provided as a template library,
which can result in a 50% performance increase vs. a shared library
implementation. See README for more information.

@tparam n_ Order of expansion to use for Chebyshev subunits
@tparam table_size_ Number of elements in lookup table which is used to assist
in finding the appropriate Chebyshev subunit.
*/
template <int n_, int table_size_> class FunctionGenerator {
    // TODO: Add complex/vector function support
  public:
    //!  Default constructor for FunctionGenerator.
    //!
    /*!

    @param f One dimensional function to approximate.
    @param a Lower bound of function domain to approximate.
    @param b Upper bound of function domain to approximate.
    @param tol Desired accuracy of approximation according to error_model.
    @param mw Minimum width of a sub-division.
    @param error_model Which model to use to calculate error for 'tol'
    parameter.
    */
    FunctionGenerator(
        std::function<double(double)> &f, double a, double b,
        double tol = 1E-12, double mw = 1E-15,
        FGError::ErrorModel error_model = FGError::ErrorModel::standard)
        : f_(f), a_(a), b_(b), tol_(tol), mw_(mw),
          scale_factor_(table_size_ / (b_ - a_)), bounds_table_(table_size_),
          error_model_(error_model) {

        init();
    };

#ifdef PYTHON_MODULE
    FunctionGenerator(py::function fpy, double a, double b, double tol = 1e-12,
                      double mw = 1e-15,
                      uint16_t error_model = FGError::ErrorModel::standard)
        : a_(a), b_(b), tol_(tol), mw_(mw),
          scale_factor_(table_size_ / (b_ - a_)), bounds_table_(table_size_),
          error_model_((FGError::ErrorModel)error_model) {
        f_ = [fpy](double x) { return fpy(x).cast<double>(); };
        init();
    };

    py::array_t<double> arr_call(py::array_t<double> x) {
        auto xin = x.unchecked<1>();
        auto res = py::array_t<double>(xin.shape(0));
        auto out = res.mutable_unchecked<1>();

        for (int i = 0; i < xin.shape(0); ++i) {
            out(i) = (*this)(xin(i));
        }

        return res;
    };
#endif

    double operator()(double x) {
        const int index = bisect_lookup(x);
        const double a = lbs_[index];
        const double b = lbs_[index + 1];
        const double xinterp = 2 * (x - a) / (b - a) - 1.0;

        return chbevl(xinterp, &coeffs_[index * n_]);
    };

  private:
    std::function<double(double)> f_;
    const double a_;
    const double b_;
    const double tol_;
    const double mw_;
    const double scale_factor_;

    std::vector<double> lbs_;
    std::vector<double> coeffs_;
    std::vector<std::pair<uint16_t, uint16_t>> bounds_table_;

    const FGError::ErrorModel error_model_;

    double chbevl(const double x, const double *c) {
        const double x2 = 2 * x;

        double c0 = c[0];
        double c1 = c[1];
        for (int i = 2; i < n_; ++i) {
            double tmp = c1;
            c1 = c[i] - c0;
            c0 = tmp + c0 * x2;
        }

        return c1 + c0 * x;
    };

    void init() {
        Eigen::MatrixXd V = calc_vandermonde();
        Eigen::PartialPivLU<Eigen::MatrixXd> VLU =
            Eigen::PartialPivLU<Eigen::MatrixXd>(V);

        fit(a_, b_, VLU);
        init_lookup();

        // Add one element to bounds_table to handle call on largest upper bound
        bounds_table_.push_back(std::make_pair(lbs_.size() - 1, lbs_.size()));

        // Fixup lbs_ to include last upper bound in case we are in last
        // subdivision for bisect_lookup
        lbs_.push_back(b_);
    }

    Eigen::MatrixXd calc_vandermonde() {
        Eigen::MatrixXd V(n_, n_);

        vecxd x = get_chebyshev_nodes(-1, 1, n_);
        for (int j = 0; j < n_; ++j) {
            V(0, j) = 1;
            V(1, j) = x(j);
        }

        for (int i = 2; i < n_; ++i) {
            for (int j = 0; j < n_; ++j) {
                V(i, j) = V(i - 1, j) * 2 * x(j) - V(i - 2, j);
            }
        }
        V = V.transpose().eval();
        return V;
    }

    vecxd get_chebyshev_nodes(double lb, double ub, int order) {
        vecxd res(order);
        for (int k = 0; k < order; ++k)
            res[order - k - 1] =
                0.5 * (lb + ub) +
                0.5 * (ub - lb) * cos(M_PI * (k + 0.5) / order);

        return res;
    }

    void fit(double a, double b, Eigen::PartialPivLU<Eigen::MatrixXd> &VLU) {
        double m = 0.5 * (a + b);
        vecxd fx = get_chebyshev_nodes(a, b, n_);

        for (int i = 0; i < n_; ++i)
            fx[i] = f_(fx[i]);

        vecxd coeffs = VLU.solve(fx);
        double tail_energy = FGError::calc_error(coeffs, error_model_);

        if (tail_energy < tol_ || (b - a) < mw_) {
            lbs_.push_back(a);
            std::vector<double> coeffstmp(n_);
            // Reverse list for cache performance reasons
            for (int i = 0; i < n_; ++i)
                coeffstmp[i] = coeffs[n_ - i - 1];

            // append to single coeffs_ vector, for cache reasons.
            for (int i = 0; i < n_; ++i)
                coeffs_.push_back(coeffstmp[i]);
        } else {
            fit(a, m, VLU);
            fit(m, b, VLU);
        }
    }

    void init_lookup() {
        for (int i = 0; i < table_size_; ++i) {
            double x0 = a_ + i * (b_ - a_) / table_size_;
            double x1 = a_ + (i + 1) * (b_ - a_) / table_size_;
            bounds_table_[i].first = bisect(x0);
            bounds_table_[i].second = bisect(x1);
        }
    }

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

    int bisect_lookup(double x) {
        int table_index = (x - a_) * scale_factor_;
        auto bisect_bounds = bounds_table_[table_index];
        return bisect_bracketed(x, bisect_bounds.first, bisect_bounds.second);
    };
};

#endif
