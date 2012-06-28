from pytriqs.demo.mymodule import myClass

class Solver:

  def Solve(self):

    myClass(self).solvecpp()


S = Solver()
S.U = 11
S.Solve()
