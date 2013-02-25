
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


from pytriqs.tools.test.pytriqs_boost_converter import *

f1( (1,2,3,4) )
f11([ [1,2,3,4], [10,20] ])

print f2()
print f3()

g( { "a" : 12, "b" : 102} )


print "g2", g2()

print "g3", g3()


p2 = p2()
print p2

p1(p2)

print "p3 = ",p3()

print "a1 = ",a1()

