#include "function_generator.hpp"

#ifdef PYTHON_MODULE
namespace py = pybind11;
#endif

#ifdef PYTHON_MODULE
py::array_t<double> FunctionGenerator::arr_call(py::array_t<double> x) {
    auto xin = x.unchecked<1>();
    auto res = py::array_t<double>(xin.shape(0));
    auto out = res.mutable_unchecked<1>();

    for (ssize_t i = 0; i < xin.shape(0); ++i) {
        out(i) = (*this)(xin(i));
    }

    return res;
};
#endif

double FunctionGenerator::operator()(double x) {
    int index = bisect_lookup(x);

    double a = lbs_[index];
    double b = ubs_[index];
    double xinterp = 2 * (x - a) / (b - a) - 1.0;

    return chbevl(xinterp, coeffs_[index]);
};

inline double FunctionGenerator::chbevl(double x, std::vector<double> &c) {
    const double x2 = 2 * x;
    double c0 = c[c.size() - 2];
    double c1 = c[c.size() - 1];

    for (size_t i = 3; i < c.size() + 1; ++i) {
        double tmp = c0;
        c0 = c[c.size() - i] - c1;
        c1 = tmp + c1 * x2;
    }

    return c0 + c1 * x;
};

FunctionGenerator::FunctionGenerator(std::function<double(double)> &f, double a,
                                     double b, double tol, int n, double mw)
    : f_(f), a_(a), b_(b), tol_(tol), n_(n), mw_(mw) {
    init();
}

void FunctionGenerator::init() {
    init_vandermonde();
    VLU_ = Eigen::PartialPivLU<Eigen::MatrixXd>(V_);

    fit(a_, b_);
    init_lookup();
}

void FunctionGenerator::init_vandermonde() {
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

vecxd FunctionGenerator::get_chebyshev_nodes(double lb, double ub, int order) {
    vecxd res(order);
    for (int k = 0; k < order; ++k)
        res[order - k - 1] =
            0.5 * (lb + ub) + 0.5 * (ub - lb) * cos(M_PI * (k + 0.5) / order);

    return res;
}

void FunctionGenerator::fit(double a, double b) {
    double m = 0.5 * (a + b);
    vecxd fx = get_chebyshev_nodes(a, b, n_);

    for (int i = 0; i < n_; ++i)
        fx[i] = f_(fx[i]);

    vecxd coeffs = VLU_.solve(fx);
    double tail_energy = standard_error(coeffs);

    if (tail_energy < tol_ || (b - a) < mw_) {
        lbs_.push_back(a);
        ubs_.push_back(b);
        std::vector<double> coefftmp(n_);
        for (int i = 0; i < n_; ++i)
            coefftmp[i] = coeffs[i];
        coeffs_.push_back(coefftmp);
    } else {
        fit(a, m);
        fit(m, b);
    }
}

void FunctionGenerator::init_lookup() {
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

#ifdef PYTHON_MODULE
FunctionGenerator::FunctionGenerator(py::function fpy, double a, double b,
                                     double tol, int n, double mw)
    : a_(a), b_(b), tol_(tol), n_(n), mw_(mw) {
    f_ = [fpy](double x) { return fpy(x).cast<double>(); };
    init();
}

PYBIND11_MODULE(FunctionGenerator, m) {
    py::class_<FunctionGenerator>(m, "FunctionGenerator")
        .def(py::init<py::function, double, double, double, int, double>(),
             py::arg("fpy") = "", py::arg("a") = 0.0, py::arg("b") = 1.0,
             py::arg("tol") = 1e-10, py::arg("n") = 12, py::arg("mw") = 1e-15)
        // .def("__call__", &FunctionGenerator::operator())
        // .def("__call__", py::vectorize(&FunctionGenerator::operator()))
        .def("__call__", &FunctionGenerator::arr_call);
    ;
}
#endif
