from pytriqs.demo.mymodule import myClass, inc

class Solver:

  def Solve(self):

    myClass(self).solvecpp()


S = Solver()
S.U = 11
S.Solve()

inc (10)

