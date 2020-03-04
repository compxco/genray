#!/usr/bin/env python

import os, sys, platform
import multiprocessing as mp
import pyctest.pyctest as pyctest
import pyctest.helpers as helpers

directory = "./"

# these are required

pyctest.PROJECT_NAME = "GENRAY"
pyctest.SOURCE_DIRECTORY = directory
pyctest.BINARY_DIRECTORY = directory

args = helpers.ArgumentParser(pyctest.PROJECT_NAME,
                              pyctest.SOURCE_DIRECTORY,
                              pyctest.BINARY_DIRECTORY).parse_args()

pyctest.BUILD_COMMAND = "make -f makefile_gfortran64"

test = pyctest.test()
test.SetName("Cmod_LH_edge")
test.SetProperty("WORKING_DIRECTORY","00_Genray_Regression_Tests/test-CMod-LH-edge")
test.SetCommand(["./test.sh"])

pyctest.run()
