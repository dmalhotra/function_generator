#include <algorithm>
#include <cmath>
#include <eigen3/Eigen/Eigen>
#include <eigen3/Eigen/LU>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include <boost/math/special_functions/chebyshev.hpp>

typedef Eigen::VectorXd vecxd;

double Vandermonde[8 * 8] = {
    1.,          1.,          1.,          1.,          1.,
    1.,          1.,          1.,          -0.98078528, -0.83146961,
    -0.55557023, -0.19509032, 0.19509032,  0.55557023,  0.83146961,
    0.98078528,  0.92387953,  0.38268343,  -0.38268343, -0.92387953,
    -0.92387953, -0.38268343, 0.38268343,  0.92387953,  -0.83146961,
    0.19509032,  0.98078528,  0.55557023,  -0.55557023, -0.98078528,
    -0.19509032, 0.83146961,  0.70710678,  -0.70710678, -0.70710678,
    0.70710678,  0.70710678,  -0.70710678, -0.70710678, 0.70710678,
    -0.55557023, 0.98078528,  -0.19509032, -0.83146961, 0.83146961,
    0.19509032,  -0.98078528, 0.55557023,  0.38268343,  -0.92387953,
    0.92387953,  -0.38268343, -0.38268343, 0.92387953,  -0.92387953,
    0.38268343,  -0.19509032, 0.55557023,  -0.83146961, 0.98078528,
    -0.98078528, 0.83146961,  -0.55557023, 0.19509032};

class FunctionGenerator {
  public:
    FunctionGenerator(double (*f)(double), double a, double b,
                      double tol = 1e-10, int n = 8, double mw = 1e-15)
        : f_(f), a_(a), b_(b), tol_(tol), n_(n), mw_(mw) {

        double r = b - a;
        double a1 = a + r / 3;
        double a2 = a + 2 * r / 3;

        x_ = get_chebyshev_nodes(a, b, n);
        V_ = Eigen::Map<Eigen::MatrixXd>(Vandermonde, 8, 8);
        VLU_ = Eigen::PartialPivLU<Eigen::MatrixXd>(V_);

        fit(a_, b_);
    };

    double operator()(double x) {
        auto el = std::lower_bound(lbs_.begin(), lbs_.end(), x);
        size_t index = std::distance(lbs_.begin(), el);
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

    Eigen::VectorXd x_;
    Eigen::MatrixXd V_;
    Eigen::PartialPivLU<Eigen::MatrixXd> VLU_;

    std::vector<double> lbs_;
    std::vector<double> ubs_;
    std::vector<vecxd> coeffs_;

    inline double chbevl(double x, vecxd &c) {
        const double x2 = 2 * x;
        double c0 = c[c.size() - 2];
        double c1 = c[c.size() - 1];

        for (int i = 3; i < c.size() + 1; ++i) {
            double tmp = c0;
            c0 = c[c.size() - i] - c1;
            c1 = tmp + c1 * x2;
        }

        return c0 + c1 * x;
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
        std::cout << "[" << a << ", " << b << "]\n";

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

    double standard_error(vecxd &coeffs) {
        double maxcoeff = 0.0;
        for (auto i = n_ - 2; i < n_; ++i)
            maxcoeff = std::max(fabs(coeffs[i]), maxcoeff);

        return maxcoeff / std::max(1.0, fabs(coeffs[0]));
    }
};

int main(int argc, char *argv[]) {
    FunctionGenerator myFunc(exp, 0, 1);

    std::cout << myFunc(0.5) << std::endl;

    return 0;
}
