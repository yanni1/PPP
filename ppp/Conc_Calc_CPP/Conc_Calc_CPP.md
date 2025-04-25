This file documents a python module built from C++ code with pybind11.
You should document the Python interfaces, *NOT* the C++ interfaces.

Module ppp.Conc_Calc_CPP
**************************

Module :py:mod:`Conc_Calc_CPP` built from C++ code in ppp.Conc_Calc_CPP
cpp`.

.. function:: add(x,y,z)
   :module: ppp.Conc_Calc_CPP
   
   Compute the sum of *x* and *y* and store the result in *z* (overwrite).

   :param x: 1D Numpy array with ``dtype=numpy.float64`` (input)
   :param y: 1D Numpy array with ``dtype=numpy.float64`` (input)
   :param z: 1D Numpy array with ``dtype=numpy.float64`` (output)
   