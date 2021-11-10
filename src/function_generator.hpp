#ifndef FUNCTIONGENERATOR_H
#define FUNCTIONGENERATOR_H
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/LU>
#include <functional>
#include <vector>

#ifdef HAVE_SCTL
#include <sctl.hpp>
#endif

#ifdef __AVX__
#include <immintrin.h>
#endif

extern "C" {
typedef double (*cfunc_double)(double);
typedef std::complex<double> (*cfunc_complex)(double);
};

//! Namespace for error calculation helper functions.
namespace FGError {
//! Model used to calculate error in approximation the input function.
enum ErrorModel : int8_t { standard = 0, relative = 1 };
using Eigen::Dynamic;
using Eigen::Matrix;

template <typename T>
inline double standard_error(Matrix<T, Dynamic, 1> &coeffs) {
    double maxcoeff = 0.0;
    int n = coeffs.size();
    for (auto i = n - 2; i < n; ++i)
        maxcoeff = std::max(std::abs(coeffs[i]), maxcoeff);

    return maxcoeff / std::max(1.0, std::abs(coeffs[0]));
}

template <typename T>
inline double relative_error(Matrix<T, Dynamic, 1> &coeffs) {
    double maxcoeff = 0.0;
    int n = coeffs.size();
    for (auto i = n - 2; i < n; ++i)
        maxcoeff = std::max(std::abs(coeffs[i]), maxcoeff);

    return maxcoeff / std::abs(coeffs[0]);
}

template <typename T>
inline double calc_error(Matrix<T, Dynamic, 1> &coeffs,
                         ErrorModel error_model) {
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
template <uint16_t n_, uint16_t table_size_, typename T = double> class FunctionGenerator {
    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> FGmatrix;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1> FGvec;

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
        cfunc_double fin, double a, double b, double tol = 1E-12, double mw = 1E-15,
        FGError::ErrorModel error_model = FGError::ErrorModel::standard)
        : a_(a), b_(b), tol_(tol), mw_(mw),
          scale_factor_(table_size_ / (b_ - a_)), bounds_table_(table_size_),
          error_model_(error_model) {
        init(fin);
    }

    T eval_xk(double x) {
        const int index = bisect_lookup(x);
        const double x0 = exp_centers_[index];
        const double xinterp = x - x0;
        return polyeval(xinterp, &coeffs_new_[index * n_]);
    }

    T eval_cheb(double x) {
        const int index = bisect_lookup(x);
        const double a = lbs_[index];
        const double b = lbs_[index + 1];
        const double xinterp = 2 * (x - a) / (b - a) - 1.0;
        return chbevl(xinterp, &coeffs_[index * n_]);
    }

    T operator()(double x) {
        return eval_xk(x);
    }

    void eval_batched_sorted(std::vector<T>& fx, const std::vector<double>& x) {
        const size_t N = x.size();
        if (fx.size() != N) fx.resize(N);

        size_t i = 0;
        if (n_ == 8) {
          #ifdef HAVE_SCTL
          constexpr int VLen = sctl::Vec<T>::Size();
          const size_t N_ = N-(VLen-1);

          auto forward_binary_search = [this,&x](const int index0, const size_t i_start, const size_t i_end, int step) {
            while (i_start+step<i_end && bisect_lookup(x[i_start+step]) == index0) step *= 2;
            step/=2;

            size_t i = i_start;
            for (; step>0; step/=2) {
              if (i+step<i_end && bisect_lookup(x[i+step])==index0) i += step;
            }
            return i+1;
          };

          while (i < N_) {
              const int index = bisect_lookup(x[i]);
              size_t b = forward_binary_search(index, i, N_, VLen);

              const auto c0 = sctl::Vec<T>(coeffs_new_[index * n_ + 0]);
              const auto c1 = sctl::Vec<T>(coeffs_new_[index * n_ + 1]);
              const auto c2 = sctl::Vec<T>(coeffs_new_[index * n_ + 2]);
              const auto c3 = sctl::Vec<T>(coeffs_new_[index * n_ + 3]);
              const auto c4 = sctl::Vec<T>(coeffs_new_[index * n_ + 4]);
              const auto c5 = sctl::Vec<T>(coeffs_new_[index * n_ + 5]);
              const auto c6 = sctl::Vec<T>(coeffs_new_[index * n_ + 6]);
              const auto c7 = sctl::Vec<T>(coeffs_new_[index * n_ + 7]);
              const auto x0 = sctl::Vec<T>(exp_centers_[index]);

              for (size_t j = i; j < b; j+= VLen) {
                const auto x1 = sctl::Vec<T>::Load(&x[j]) - x0;
                const auto x2 = x1*x1;
                const auto x4 = x2*x2;
                const auto fx_ = FMA(x4, FMA(x2, FMA(x1, c7, c3), FMA(x1, c6, c2)), FMA(x2, FMA(x1, c5, c1), FMA(x1, c4, c0)));
                fx_.Store(&fx[j]);
              }
              i = b;
          }
          #endif
        }
        for (; i < N; i++) {
            fx[i] = (*this)(x[i]);
        }
    }

  private:
    const double a_;
    const double b_;
    const double tol_;
    const double mw_;
    const double scale_factor_;

    std::vector<double> lbs_;
    std::vector<double> exp_centers_;
    std::vector<T> coeffs_;
    std::vector<T> coeffs_new_;
    std::vector<std::pair<uint16_t, uint16_t>> bounds_table_;

    const FGError::ErrorModel error_model_;

    T chbevl(const double x, const T *c) {
        const double x2 = 2 * x;

        T c0 = c[0];
        T c1 = c[1];
        for (int i = 2; i < n_; ++i) {
            T tmp = c1;
            c1 = c[i] - c0;
            c0 = tmp + c0 * x2;
        }

        return c1 + c0 * x;
    }
    T polyeval(const double x, const T *c) {
        if (n_ == 8) {
            #ifdef __AVX__
            alignas(sizeof(double)*4) double y[4];
            __m256d c0_ = _mm256_loadu_pd(c+0);
            __m256d c1_ = _mm256_loadu_pd(c+4);
            __m256d x_ = _mm256_set1_pd(x);
            __m256d c0_c1x = _mm256_add_pd(c0_, _mm256_mul_pd(c1_, x_));
            _mm256_store_pd(y, c0_c1x);

            const double x2 = x*x;
            const double x4 = x2*x2;
            return (y[0] + y[1]*x2) + (y[2] + y[3]*x2)*x4;
            #else
            const double x2 = x*x;
            const double x4 = x2*x2;
            return (c[0]+c[4]*x + (c[1]+c[5]*x)*x2) + (c[2]+c[6]*x + (c[3]+c[7]*x)*x2) * x4;
            #endif
        } else { // generic n_
            const double* c0_ = c;
            const double* c1_ = c + n_/2;
            double x2 = x*x, x_ = 1, sum = 0;
            for (int i = 0; i < n_/2; i++) {
                sum += (c0_[i] + c1_[i] * x) * x_;
                x_ *= x2;
            }
            return sum;
        }
    }

    void init(cfunc_double f) {
        FGmatrix V = calc_vandermonde();
        Eigen::PartialPivLU<FGmatrix> VLU = Eigen::PartialPivLU<FGmatrix>(V);

        fit(f, a_, b_, VLU);
        init_lookup();

        // Add one element to bounds_table to handle call on largest upper bound
        bounds_table_.push_back(std::make_pair(lbs_.size() - 1, lbs_.size()));

        // Fixup lbs_ to include last upper bound in case we are in last
        // subdivision for bisect_lookup
        lbs_.push_back(b_);
    }

    FGmatrix calc_vandermonde() {
        FGmatrix V(n_, n_);

        auto x = get_chebyshev_nodes(-1, 1, n_);
        for (int j = 0; j < n_; ++j) {
            V(0, j) = 1;
            V(1, j) = x(j);
        }

        for (int i = 2; i < n_; ++i) {
            for (int j = 0; j < n_; ++j) {
                V(i, j) = T(2) * V(i - 1, j) * x(j) - V(i - 2, j);
            }
        }
        V = V.transpose().eval();
        return V;
    }

    Eigen::VectorXd get_chebyshev_nodes(double lb, double ub, int order) {
        Eigen::VectorXd res(order);
        for (int k = 0; k < order; ++k)
            res[order - k - 1] =
                0.5 * (lb + ub) +
                0.5 * (ub - lb) * cos(M_PI * (k + 0.5) / order);

        return res;
    }

    void fit(cfunc_double f, double a, double b, Eigen::PartialPivLU<FGmatrix> &VLU) {
        double m = 0.5 * (a + b);
        auto xvec = get_chebyshev_nodes(a, b, n_);
        FGvec yvec(n_);

        for (int i = 0; i < n_; ++i)
            yvec[i] = f(xvec[i]);

        FGvec coeffs = VLU.solve(yvec);
        double tail_energy = FGError::calc_error<T>(coeffs, error_model_);

        if (tail_energy < tol_ || (b - a) < mw_) {
            lbs_.push_back(a);
            std::vector<T> coeffstmp(n_);
            // Reverse list for cache performance reasons
            for (int i = 0; i < n_; ++i)
                coeffstmp[i] = coeffs[n_ - i - 1];

            // append to single coeffs_ vector, for cache reasons.
            for (int i = 0; i < n_; ++i)
                coeffs_.push_back(coeffstmp[i]);

            { // set coeff_new_
              exp_centers_.push_back((a+b)*0.5);

              auto getVandermonde = [this]() {
                FGmatrix V(n_, n_);
                auto x = get_chebyshev_nodes(-1, 1, n_);
                for (int i = 0; i < n_; ++i) {
                    for (int j = 0; j < n_; ++j) {
                        V(i, j) = pow(x(i),j);
                    }
                }
                return V;
              };
              static const auto V = getVandermonde();
              static const auto VLU = Eigen::PartialPivLU<FGmatrix>(V);

              FGvec coeffs_new = VLU.solve(yvec);
              for (int i = 0; i < n_; ++i)
                  coeffs_new[i] *= pow((b-a)*0.5,-i); // include scaling here to avoid remapping xinterp to (-1,1)

              // store even and odd coefficients separately for vectorization
              static_assert(n_%2 == 0, "n_ must be even.");
              for (int i = 0; i < n_/2; ++i) coeffs_new_.push_back(coeffs_new[i*2+0]);
              for (int i = 0; i < n_/2; ++i) coeffs_new_.push_back(coeffs_new[i*2+1]);
            }

            // Require that number of subdivisions is less than maximum possible value in
            // lookup table. If you are failing here, then use a domain of your function that
            // is more well behaved.
            assert((coeffs_.size() / n_ < UINT16_MAX) &&
                   "Too many subdivisions. See comment here for details.");
        } else {
            fit(f, a, m, VLU);
            fit(f, m, b, VLU);
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

    uint16_t bisect_bracketed(double x, int n1, int n2) {
        while (n2 - n1 > 1) {
            const uint16_t m = n1 + (n2 - n1) / 2;
            if (x < lbs_[m])
                n2 = m;
            else
                n1 = m;
        }

        return n1;
    }

    uint16_t bisect(double x) { return bisect_bracketed(x, 0, lbs_.size()); }

    uint16_t bisect_lookup(double x) {
        int table_index = (x - a_) * scale_factor_;
        auto bisect_bounds = bounds_table_[table_index];
        return bisect_bracketed(x, bisect_bounds.first, bisect_bounds.second);
    }
};


#endif
