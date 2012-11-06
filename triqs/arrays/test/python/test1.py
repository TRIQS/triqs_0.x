
################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2011 by O. Parcollet
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

import numpy
import sys

a=numpy.array([[1,2],[3,4]]) #,[5,6]])

aa = numpy.array( [[0,0,0,0], [0,1,2,0],[0,3,4,0],[0,0,0,0]])
a = aa[1:3,1:3]

print a

from _array_tests import *
print test1(a)

print print_array_i(a)


al=[[1,2],[3,4]] 
a=numpy.array(al)
#print_array_f(al)
#print_array_f(a)

try: 
    print_array_view_f(al)
except:
    pass
try :
    print_array_view_f(a)
except : pass

print make_array()

#print " ============="


#sys.exit(0)

print "testing range conversion"
print test3( slice (0,4,2))

def run_test(a,s1,s2,test_f) : 
      #print a[s1,s2] 
      #print test_f(a,s1,s2)
      if not ((a[s1,s2] == test_f(a,s1,s2)).all()) :
        print s1, s2
        print "orig : ", a[s1,s2]
        print "C ++ : ",test_f(a,s1,s2)
        raise RunTimeError, "test failed"

b = numpy.array([[1,2,3,4],[10,20,30,40],[15,25,35,45],[17,27,37,47]])
sli_list = [ slice(i,j,s) for i in range (0,4) for j in range (i,4) for s in [1,2]]
print "testing slicing return"
#print b

print "Running 2d slice tests "
for s1 in sli_list : 
  for s2 in sli_list : 
      #test1(b)
      #print "b,ref in ", b, sys.getrefcount(b)
      run_test(b,s1,s2,test2)
      #print "b,ref out ", b, sys.getrefcount(b)
print " ... ok"

print "Running 1d slice tests "
for s1 in sli_list : 
  for j in range(4) :
    run_test(b,s1,j,test4)
    run_test(b,j,s1,test5)
print " ... ok"

a,v = test6()
print a
print v
v[0,0] = 19
print a
print v

