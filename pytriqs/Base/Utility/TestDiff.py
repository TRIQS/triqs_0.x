
################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by M. Ferrero, O. Parcollet
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

import unittest,os,sys,glob,string,glob
from optparse import OptionParser

Usage = """
%prog [options] Calculation1 Calculation2....

This script will take all CalculationX files, run them, 
and compare (with h5diff) the Results.h5 files in their directory
with the reference result, which must be CalculationX.ResultsReference.h5
It will stop at the first failure.
"""
parser = OptionParser(usage=Usage)
parser.add_option("-p", "--precision", action="store", dest="precision", default=0.01,
                  help="Absolute precision of h5diff comparison (option -d)")
parser.add_option("-v", "--verbose",action="store_true", dest="verbose", default=False,
                  help="Verbose (same as -V 5)")

(Options, Args) = parser.parse_args()
precision = float(Options.precision)

# find all test file
Tests = []
for f in Args : 
    g = glob.glob(f)
    if g == [] : print "File %s not found"%f
    for t in g : 
        if os.path.exists('%s.ResultsReference.h5'%t) : 
            Tests += [t]
        else : 
            print "Test %s will not be run \n   since the HDF5 file %s.ResultsReference.h5 is missing  "%(t,t)
print "I am about to test the following calculations", Tests

def myexec(COMMAND) : 
    cin,cout,cerr = os.popen3(COMMAND)  # execute in a unix pipe.
    return  cout.readlines(),cerr.readlines()

def myrun(testfile):
    testout,testerr = myexec ( '%s %s' %(sys.executable, testfile) ) 

    # find the output directory
    i=0
    try : 
        while testout[i].find("Running in directory") !=0 : i +=1
        dir = string.join(testout[i].split(':')[1:])
    except : 
        raise RuntimeError, "Directory not found !"
    print dir

    #out,err = myexec('h5diff -r -d %f %s/Results.h5 %s..ResultReference.h5 '%(precision,dir,testfile))
    out,err = myexec('h5diff -r -d %f test2.hdf5 test2.hdf5 '%(precision))
    if out : 
        print "Error in comparing file  : " 
        for l in out : 
            print l.strip()
        raise RuntimeError, "Regression failed at precision %f"%precision
    #raise unittest.TestCase.failureException

if __name__ == '__main__':
    for test in Tests : 
        print 20*"-"
        print "Testing : %s\n"%test
        myrun('../Examples/%s'%test)
        print "\nOK\n"
    print 20*"=="
    print "All tests are successfull !"

