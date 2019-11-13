#include "function_generator.hpp"

#ifdef PYTHON_MODULE
namespace py = pybind11;

#define ORDER 8
#define TABLE_SIZE 4096

PYBIND11_MODULE(FunctionGenerator, m) {
    py::class_<FunctionGenerator<ORDER, TABLE_SIZE, double>>(m, "FunctionGenerator")
        .def(py::init<py::function, double, double, double, double, uint16_t>(),
             py::arg("fpy") = "", py::arg("a") = 0.0, py::arg("b") = 1.0,
             py::arg("tol") = 1e-10, py::arg("mw") = 1e-15,
             py::arg("error_model") = (uint16_t) FGError::ErrorModel::standard)
        .def("__call__", &FunctionGenerator<ORDER, TABLE_SIZE, double>::operator())
        .def("__call__", &FunctionGenerator<ORDER, TABLE_SIZE, double>::arr_call);
    ;
}
#undef ORDER
#undef TABLE_SIZE
#endif
