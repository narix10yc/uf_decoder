#include "pybind11/pybind11.h"
#include "pybind11/stl.h"

#include "UF/UFDecoder.h"

#include <sstream>

namespace py = pybind11;

namespace {

template <typename T> void bind_BitArray(py::module_& m, const char* name) {
  py::class_<T>(m, name)
      .def(py::init<>())
      .def("set", &T::set)
      .def("unset", &T::unset)
      .def("flip", &T::flip)
      .def("test", &T::test)
      .def("clear", &T::clear)
      .def("set_all", &T::set_all)
      .def("set_word", &T::set_word)
      .def("num_bits",
           [](const T&) {
             return T::num_bits();
           }) // num_bits is consteval
      .def("num_words",
           [](const T&) {
             return T::num_words();
           }) // num_words is consteval
      .def("__str__", [](const T& self) {
        std::ostringstream oss;
        self.display(oss);
        return oss.str();
      });
}

template <int L> void bind_decoder(py::module_& m) {
  using SyndromeType = typename uf::ToricCode<L>::Syndrome;
  using ErrorType = typename uf::ToricCode<L>::Error;
  std::string syndromeTypeName = "Syndrome" + std::to_string(L);
  std::string errorTypeName = "Error" + std::to_string(L);
  bind_BitArray<SyndromeType>(m, syndromeTypeName.c_str());
  bind_BitArray<ErrorType>(m, errorTypeName.c_str());

  std::string methodName = "decode_" + std::to_string(L);
  m.def(methodName.c_str(), &uf::decode<L>, py::arg("syndromes"));
  std::string computeLogicalOpName = "compute_logical_op_" + std::to_string(L);
  m.def(computeLogicalOpName.c_str(),
        &uf::compute_logical_op<L>,
        py::arg("error"),
        py::arg("r0") = 0, // default reference row for X_h
        py::arg("c0") = 0  // default reference column for X_v
  );
}
} // anonymous namespace

PYBIND11_MODULE(uf_python, m) {
  m.doc() = "UF Python bindings";
  // Because we made grid size a compile-time constant (for performance),
  // we need to bind each grid size needed.
  bind_decoder<3>(m);
  bind_decoder<5>(m);
  bind_decoder<7>(m);
}