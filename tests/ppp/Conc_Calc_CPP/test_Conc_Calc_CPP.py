#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Tests for C++ module ppp.Conc_Calc_CPP.
"""

import sys
sys.path.insert(0,'.')

import numpy as np

import ppp.Conc_Calc_CPP as Conc_Calc_CPP


def test_cpp_add():
    x = np.array([0,1,2,3,4],dtype=float)
    shape = x.shape
    y = np.ones (shape,dtype=float)
    z = np.zeros(shape,dtype=float)
    expected_z = x + y
    result = Conc_Calc_CPP.add(x,y,z)
    assert (z == expected_z).all()


#===============================================================================
# The code below is for debugging a particular test in eclipse/pydev.
# (normally all tests are run with pytest)
#===============================================================================
if __name__ == "__main__":
    the_test_you_want_to_debug = test_cpp_add

    print(f"__main__ running {the_test_you_want_to_debug} ...")
    the_test_you_want_to_debug()
    print('-*# finished #*-')
#===============================================================================
