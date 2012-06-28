from pytriqs.demo.mymodule import *

class Solver:
  def Solve(self):
    myClass(self).solvecpp()



S = Solver()
S.U = 11
S.Solve()

inc (10)

d = dict ( i = 2, s = "a nice string ...", l = [1,2,3])
d['C'] =  myClass (S)


read_dict(d)

print modify_dict (d)
print d

C= make_myclass(19)
print C, C.U
