from pytriqs.solvers.ctqmc_hyb import Solver
import string
import textwrap

kl = max ( len(k) for k in Solver.Optional) + 4
dl = max ( len(str(defaut)) for k, (doc,defaut,typ) in Solver.Optional.items() )  + 4
docl = 60

l = (kl-2)*"=" + "  " + (dl-2)*"=" + "  " + 10 *"=" + "  " + docl*"=" +  '\n'
s = l
def add (k, defaut, typ, doc) : 
    L = textwrap.wrap( doc.strip(), docl)
    r =  k + (kl-len(k))*" " +  defaut + (dl - len(defaut) ) *" " +  typ  + (12 - len(typ))*" " + L[0] +  '\n' 
    for l in L[1:] : 
        r += (kl + dl + 12)*" " + l + '\n'
    return r

s += add("Key", "Default", "Type", "Documentation")
s += l
for k, (doc,defaut,typ) in Solver.Optional.items(): 
    s += add( k, str(defaut) , typ.__name__, doc)
s += '\n' + l

print s


kl = max ( len(k) for k in Solver.Required ) + 4
docl = 60

l = (kl-2)*"=" + "  " +  10 *"=" + "  " + docl*"=" +  '\n'
s = l
def add2 (k, typ, doc) : 
    L = textwrap.wrap( doc.strip(), docl)
    r =  k + (kl-len(k))*" " +   typ  + (12 - len(typ))*" " + L[0] +  '\n' 
    for l in L[1:] : 
        r += (kl + dl + 12)*" " + l + '\n'
    return r

s += add2("Key", "Type", "Documentation")
s += l
for k, (doc,typ) in Solver.Required.items(): 
    s += add2( k,  typ.__name__, doc)
s += '\n' + l

print s






