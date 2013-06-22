#!/usr/bin/env python

from numpy.distutils.core import setup, Extension
import os
setup(name = "gadgetPydensity",version = "1.0",ext_modules = [Extension("gadgetPydensity", ["densitymake.c"])])


