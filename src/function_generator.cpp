#include <algorithm>
#include <chrono>
#include <cmath>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/LU>
#include <iostream>
#include <random>
#include <string>
#include <vector>

typedef Eigen::VectorXd vecxd;

class FunctionGenerator {
  public:
    FunctionGenerator(double (*f)(double), double a, double b,
                      double tol = 1e-10, int n = 12, double mw = 1e-15)
        : f_(f), a_(a), b_(b), tol_(tol), n_(n), mw_(mw) {

        init_vandermonde();
        VLU_ = Eigen::PartialPivLU<Eigen::MatrixXd>(V_);

        fit(a_, b_);
        init_lookup();
    }

    double operator()(double x) {
        int index = bisect_lookup(x);

        double a = lbs_[index];
        double b = ubs_[index];
        double xinterp = 2 * (x - a) / (b - a) - 1.0;

        return chbevl(xinterp, coeffs_[index]);
    }

  private:
    double (*f_)(double);
    const double a_;
    const double b_;
    const double tol_;
    const int n_;
    const double mw_;

    Eigen::MatrixXd V_;
    Eigen::PartialPivLU<Eigen::MatrixXd> VLU_;

    std::vector<double> lbs_;
    std::vector<double> ubs_;
    std::vector<vecxd> coeffs_;
    std::vector<std::pair<uint16_t, uint16_t>> bounds_table_;

    inline double chbevl(double x, vecxd &c) {
        const double x2 = 2 * x;
        double c0 = c[c.size() - 2];
        double c1 = c[c.size() - 1];

        for (int i = 3; i < n_ + 1; ++i) {
            double tmp = c0;
            c0 = c[c.size() - i] - c1;
            c1 = tmp + c1 * x2;
        }

        return c0 + c1 * x;
    }

    void init_vandermonde() {
        V_.resize(n_, n_);
        vecxd x = get_chebyshev_nodes(-1, 1, n_);
        for (int j = 0; j < n_; ++j) {
            V_(0, j) = 1;
            V_(1, j) = x(j);
        }

        for (int i = 2; i < n_; ++i) {
            for (int j = 0; j < n_; ++j) {
                V_(i, j) = V_(i - 1, j) * 2 * x(j) - V_(i - 2, j);
            }
        }
        V_ = V_.transpose().eval();
    }

    vecxd get_chebyshev_nodes(double lb, double ub, int order) {
        vecxd res(order);
        for (int k = 0; k < order; ++k)
            res[order - k - 1] =
                0.5 * (lb + ub) +
                0.5 * (ub - lb) * cos(M_PI * (k + 0.5) / order);

        return res;
    }

    void fit(double a, double b) {
        double m = 0.5 * (a + b);
        vecxd fx = get_chebyshev_nodes(a, b, n_);

        for (int i = 0; i < n_; ++i)
            fx[i] = f_(fx[i]);

        vecxd coeffs = VLU_.solve(fx);
        double tail_energy = standard_error(coeffs);

        if (tail_energy < tol_ || (b - a) < mw_) {
            lbs_.push_back(a);
            ubs_.push_back(b);
            coeffs_.push_back(coeffs);
        } else {
            fit(a, m);
            fit(m, b);
        }
    }

    void init_lookup() {
        constexpr int table_size = 512;
        bounds_table_.resize(table_size);

        // FIXME: This screws up sometimes. Probably an off by one error or
        // rounding issue
        for (int i = 0; i < table_size; ++i) {
            double x0 = a_ + i * (b_ - a_) / table_size;
            double x1 = a_ + (i + 1) * (b_ - a_) / table_size;
            bounds_table_[i].first = bisect(x0);
            bounds_table_[i].second = bisect(x1);
        }
    }

    inline int bisect_bracketed(double x, int n1, int n2) {
        while (n2 - n1 > 1) {
            const int m = n1 + (n2 - n1) / 2;
            if (x < lbs_[m])
                n2 = m;
            else
                n1 = m;
        }

        return n1;
    }

    inline int bisect(double x) { return bisect_bracketed(x, 0, lbs_.size()); }

    inline int bisect_cache(double x) {
        constexpr int half_width_cache = 4;
        static int n1 = half_width_cache;
        n1 -= half_width_cache;
        n1 = std::max(n1, 0);
        int n2 = n1 + 2 * half_width_cache;

        while (lbs_[n1] > x)
            n1 /= 2;
        while (n2 < (int) lbs_.size() && lbs_[n2] < x)
            n2 *= 2;
        n2 = std::min(n2, (int)lbs_.size());

        n1 = bisect_bracketed(x, n1, n2);
        return n1;
    }

    int bisect_lookup(double x) {
        int table_index = x / (b_ - a_) * bounds_table_.size();
        auto bisect_bounds = bounds_table_[table_index];
        return bisect_bracketed(x, bisect_bounds.first, bisect_bounds.second);
    }

    double standard_error(vecxd &coeffs) {
        double maxcoeff = 0.0;
        for (auto i = n_ - 2; i < n_; ++i)
            maxcoeff = std::max(fabs(coeffs[i]), maxcoeff);

        return maxcoeff / std::max(1.0, fabs(coeffs[0]));
    }
};

int main(int argc, char *argv[]) {
    using namespace std::chrono;
    high_resolution_clock::time_point start = high_resolution_clock::now();

    double range[2] = {1e-10, 1000};
    FunctionGenerator myLog(log, range[0], range[1]);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<> dis(range[0], range[1]);
    size_t n_el = 100000000;
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
        high_resolution_clock::time_point start = high_resolution_clock::now();

        for (size_t i = 0; i < n_el; ++i) {
            log(x[i]);
        }
        high_resolution_clock::time_point finish = high_resolution_clock::now();
        duration<double> time_span =
            duration_cast<duration<double>>(finish - start);

        std::cout << "System 'log' took " << time_span.count() << " seconds.\n";
    }

    {
        high_resolution_clock::time_point start = high_resolution_clock::now();

        for (size_t i = 0; i < n_el; ++i) {
            myLog(x[i]);
        }
        high_resolution_clock::time_point finish = high_resolution_clock::now();
        duration<double> time_span =
            duration_cast<duration<double>>(finish - start);

        std::cout << "Chebyshev 'log' took " << time_span.count()
                  << " seconds.\n";
        std::cout << "Efficiency: " << n_el / time_span.count() / 1e6
                  << " million evaluations / sec\n";
    }

    return 0;
}
