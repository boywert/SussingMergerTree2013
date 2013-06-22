#!/usr/bin/env python

from numpy.distutils.core import setup, Extension
import os
setup(name = "gadgetPyIO",version = "1.0",ext_modules = [Extension("gadgetPyIO", ["gadgetmodule.c"])])


