# -*- coding: utf-8 -*-

"""Tests for ppp package."""

import sys
sys.path.insert(0,'.')

import ppp


def test_hello_default_arg():
    result = ppp.hello()
    assert result == "Hello world!"


def test_hello_me():
    result = ppp.hello('me')
    assert result == "Hello me!"


# ==============================================================================
# The code below is for debugging a particular test in eclipse/pydev.
# (otherwise all tests are normally run with pytest)
# Make sure that you run this code with the project directory as CWD, and
# that the source directory is on the path
# ==============================================================================
if __name__ == "__main__":
    the_test_you_want_to_debug = test_hello_default_arg

    print("__main__ running", the_test_you_want_to_debug)
    the_test_you_want_to_debug()
    print('-*# finished #*-')

# eof