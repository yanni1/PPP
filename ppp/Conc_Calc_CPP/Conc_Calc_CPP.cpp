/*
 *  C++ source file for module ppp.Conc_Calc_CPP
 */


// See http://people.duke.edu/~ccc14/cspy/18G_C++_Python_pybind11.html for examples on how to use pybind11.
// The example below is modified after http://people.duke.edu/~ccc14/cspy/18G_C++_Python_pybind11.html#More-on-working-with-numpy-arrays
#include <nanobind/nanobind.h>
#include <nanobind/ndarray.h> // add support for multi-dimensional arrays
#include <vector>
#include <omp.h>
namespace nb = nanobind;




void
add
  ( nb::ndarray<double> a // in
  , nb::ndarray<double> b // in
  , nb::ndarray<double> c // inout
  )
{
    size_t n = 1;
    for (size_t idim = 0; idim < a.ndim(); ++idim) {
        n *= a.shape(idim);
    }
    double const * a_ = a.data();
    double const * b_ = b.data();
    double       * c_ = c.data();
    for(size_t i = 0; i < n; ++i) {
        c_[i] = a_[i] + b_[i];
    }
}


NB_MODULE(Conc_Calc_CPP, m) {
    m.doc() = "A simple example python extension";

    m.def("add", &add, "Add two np.float64 arrays.");

    m.def("inspect"
         , [](nb::ndarray<> a)
           {
                printf("Array data pointer : %p\n", a.data());
                printf("Array dimension : %zu\n", a.ndim());
                for (size_t i = 0; i < a.ndim(); ++i) {
                    printf("Array dimension [%zu] : %zu\n", i, a.shape(i));
                    printf("Array stride    [%zu] : %zd\n", i, a.stride(i));
                }
                printf("Device ID = %u (cpu=%i, cuda=%i)\n"
                , a.device_id()
                , int(a.device_type() == nb::device::cpu::value)
                , int(a.device_type() == nb::device::cuda::value)
                );
                printf("Array dtype: int16=%i, uint32=%i, float32=%i, float64=%i\n"
                , a.dtype() == nb::dtype<int16_t>()
                , a.dtype() == nb::dtype<uint32_t>()
                , a.dtype() == nb::dtype<float>()
                , a.dtype() == nb::dtype<double>()
                );
            }
          , "inspect an array"
    );
}