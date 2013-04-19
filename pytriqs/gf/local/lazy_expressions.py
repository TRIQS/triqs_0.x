
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

__all__ = ['LazyExpr', 'LazyExprTerminal', 'eval_expr_with_context', 'lazy', 'lazy_function', 'transform', 'eval_expr']

class __aux(object): 

    def __add__(self, y): return LazyExpr("+", LazyExpr(self), LazyExpr(y))
    def __sub__(self, y): return LazyExpr("-", LazyExpr(self), LazyExpr(y))
    def __mul__(self, y): return LazyExpr("*", LazyExpr(self), LazyExpr(y))
    def __div__(self, y): return LazyExpr("/", LazyExpr(self), LazyExpr(y))

    def __radd__(self, y): return LazyExpr("+", LazyExpr(y), LazyExpr(self))
    def __rsub__(self, y): return LazyExpr("-", LazyExpr(y), LazyExpr(self))
    def __rmul__(self, y): return LazyExpr("*", LazyExpr(y), LazyExpr(self))
    def __rdiv__(self, y): return LazyExpr("/", LazyExpr(y), LazyExpr(self))

    def __iadd__(self, y): return self.set_from(self+y)
    def __isub__(self, y): return self.set_from(self-y)
    def __imul__(self, y): return self.set_from(self*y)
    def __idiv__(self, y): return self.set_from(self/y)

    def __call__(self, *args): return LazyExpr("F", make_lazy(self), *map(make_lazy, args))

class LazyExprTerminal (__aux):
    pass 

class LazyExpr (__aux): 
    """
    """

    def __init__ (self, *args): 
        if len(args) == 1: 
            a0 = args[0]
            self.tag, self.childs = (a0.tag, a0.childs) if isinstance(a0, self.__class__) else ("T", [a0])
        elif len(args) >1: 
            self.tag, self.childs = args[0], args[1:]
        else: raise ValueError, "too few arguments"
    
    def copy(self): 
        """ Deep copy"""
        return LazyExpr(self.tag, self.childs)

    def set_from(self, y):
        """ self:= y """
        self.tag, self.childs = tmp.tag, tmp.childs
        return self

    def is_terminal(self): 
        """Returns true iif the expression is a terminal  """
        return self.tag=="T"
 
    def get_terminal(self): 
        """Returns the terminal if the expression is a terminal else None """
        return self.childs[0] if self.tag=="T" else None
 
    def __aux_print(self, F): 
        op_priority = {'T': 100, "+": 1 , '-': 1, '*': 2, '/': 2}
        if self.tag == "T":  return F(self.childs[0])
        if self.tag == "F":  
            return reduce (lambda s, e: s+ F(e), self.childs[1:], self.childs[0].get_terminal()[0] + "(" ) + ')'
        par = lambda op, e: "%s"%e if op_priority[e.tag] >= op_priority[op] else "(%s)"%e
        return "%s %s %s "%(par(self.tag , self.childs[0]), self.tag , par(self.tag , self.childs[1]))

    def __str__(self): return self.__aux_print(str)
    def __repr__(self): return self.__aux_print(repr)
   
    #def __call__ (self, *args, **kwargs): 


#-----------------------------------------------------
 
def eval_expr_with_context(eval_term, expr ): 

    if expr.tag == "T": return eval_term(expr.childs[0]) #eval the terminals

    if expr.tag == "F":
        f= expr.childs[0].get_terminal()[1]
        return f (*map(lambda e:eval_expr_with_context(eval_term, e) , expr.childs[1:]) ) 
      
    # Binary operations: 
    ops = { "+": lambda x, y: x + y, "-": lambda x, y: x - y, "*": lambda x, y: x * y, "/": lambda x, y: x / y }
    return ops[expr.tag] (*map(lambda e:eval_expr_with_context(eval_term, e) , expr.childs) ) 

#-----------------------------------------------------

def make_lazy(x): return LazyExpr(x)

#-----------------------------------------------------

def lazy_function(name, F): 
    return LazyExpr("T", (name, F))  
 
#-----------------------------------------------------

def transform (expr, Fnode, Fterm = lambda x: x ): 
    """Given two functions   
           Fnode(tag, childs) -> (tag, childs)
           Fterm(x) -> x'
           it transforms the expression recursively
    """
    if expr.tag == "T": return LazyExpr("T", Fterm(expr.childs[0]))
    tag, ch = Fnode (expr.tag, map(lambda e: transform (e, Fnode), expr.childs))
    ch = map(lambda x: LazyExpr(x), ch)
    return LazyExpr (tag, *ch)

#-----------------------------------------------------

def all_terminals (expr): 
    """Generate all terminals of an expression"""
    if expr.tag == "T": 
        yield expr.childs[0]
    else: 
        for ch in expr.childs: 
            for t in all_terminals(ch): 
                yield t

def eval_expr (expr):
#def eval_expr_or_pass (expr):
    """ 
    If expr is not a LazyExpr: returns expr unchanged.
    Otherwise, tries to eval it by looking for some element in the tree that can create the evaluation context and is not purely abstract
    """
    if not isinstance (expr, LazyExpr): return expr # do nothing 
    # first take all terminals
    C = [ t.__lazy_expr_eval_context__() for t in all_terminals(expr) if hasattr(t, "__lazy_expr_eval_context__") ]
    if C == []: raise ValueError, "Evaluation impossible: expression is purely abstract"
    all_equal = reduce (lambda x, y: x and y , [ C[0] == x for x in C ])
    if not all_equal: raise ValueError, "Evaluation impossible: various terminals lead to incompatible evaluation contexts: their type are not compatible for binary ops"
    C = C[0]
    return eval_expr_with_context(C, expr) 

#--------------- TEST --------------------------------------

if __name__ == "__main__":

    class T (LazyExprTerminal): 
        def __init__(self, n): self.name = n
        def __repr__(self): return self.name

    g1, g2 = T("g1"), T("g2")
    a = 2 * g1 + g2/3 + 5
    print a

    def e_t(x): 
        d =  { "g1": 10, "g2": 100}
        return d[x.name] if isinstance(x, T) else x

    assert eval_expr_with_context(e_t, a) ==58

    def find_sca(tag, childs): 
        if tag =="+":
            t = childs[1].get_terminal()
            if t: childs[1] =  T("$$" + str(childs[1]))
        return (tag, childs)

    b= transform(a, find_sca) 
    print b

    # a lazy function: 
    def f (x): return -x

    fa = lazy_function("f", f) (a)
    print fa, eval_expr_with_context(e_t, fa) 
   
    assert eval_expr_with_context(e_t, fa) == f(eval_expr_with_context(e_t, a) )
    

