
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

__all__ = ['plt','oplot','subplots','figsize_default','use_amsmath']

import numpy, matplotlib as mpl, matplotlib.pylab as plt
from protocol import plot_protocol_apply
from pytriqs.Base.GF_Local.lazy_expressions import eval_expr_or_pass
from matplotlib import rc

try:
  subplots = mpl.pyplot.subplots
except:
  def subplots(nrows=1, ncols=1, sharex=False, sharey=False, squeeze=True, subplot_kw=None, **fig_kw):
    print "subplots not supported"
    return plt.figure(1), [plt.subplot(nrows,ncols,x+1) for x in range(nrows*ncols)]

figsize_default = (12,8)


def oplot (*ob_list, **OptionsDict) : 
    """
    A thin layer above pyplot.plot function that allows plotting objects with
    plot protocol as well as arrays. 
    Options are the same as for the pyplot.plot function.
    """
    plt.figure(1, figsize = OptionsDict.pop('figsize', figsize_default ) )
    __oplot_impl(plt.plot, plt.xlabel,plt.ylabel,plt.legend, *ob_list,**OptionsDict)
    if hasattr(plt.figure(1), "show") : plt.figure(1).show()

mpl.axes.Axes.oplot = lambda self, *ob_list, **OptionsDict : __oplot_impl(self.plot,self.set_xlabel, self.set_ylabel, self.legend, *ob_list,**OptionsDict)

def __oplot_impl (plotFnt,xlabelFnt, ylabelFnt, legendFnt, *ob_list, **OptionsDict) : 
    """
    A thin layer above pyplot.plot function that allows plotting objects with
    plot protocol as well as arrays. 
    Options are the same as for the pyplot.plot function.
    """

    def objs() : # filter the arguments for the format strings ... 
        i, l = 0, []
        while i< len (ob_list) : 
            if i < len(ob_list) - 1 and type(ob_list[i+1]) == type("") : 
                res = ob_list[i], [ ob_list[i+1] ]
                i+=2
            else :  
                res =  ob_list[i], [ ]
                i+=1
            yield res

    for ob, OptionsList in objs() :
        opt = OptionsDict.copy() # the plot protocol will consume the dict....
        #ob2 = eval_expr_or_pass (ob) # if it is a lazy_expr, it is time to evaluate it !
        for curvedata in plot_protocol_apply(ob,opt, plt.xlim ) : 
            X,Y = curvedata['xdata'],curvedata['ydata']
            d = { 'label' :  curvedata['label'] } 
            d.update(opt)
            try : 
                plotFnt(X,Y,*OptionsList,**d)
            except TypeError, e:
                import re
                m = re.search('(?<=There is no line property )"(.*)"', str(e) )
                if m : 
                   raise RuntimeError, "Option %s is not understood in plot function : it is not an option of the object to be plotted, nor a matplotlib option"%m.group(0)
                else : 
                   raise 
            if 'xlabel' in curvedata : xlabelFnt(curvedata['xlabel'], fontsize=20) 
            if 'ylabel' in curvedata : ylabelFnt(curvedata['ylabel'], fontsize=20) 

    legendFnt(loc = 1) #legend is built from the label

def use_amsmath():
  rc('text', usetex=True)
  rc('text.latex', preamble="\usepackage{amsmath}")


#------------------- OBSOLETE ?? ------------------------

class Plotter_OneGraph (object): 
    """ This is one of the plots appearing in the Plotter figure. """
    def __init__(self, *coord) :
        """ Constructs the graph at coordinate (i,j) in the panel """
        self.coord = coord
        for method in ['axis','xlim','ylim','grid','text','draw','setp', 'xlabel', 'ylabel'] : 
            setattr(self,method, getattr(plt,method))
        plt.subplot(*coord)

    def _focuss(self) : 
        plt.subplot(*self.coord)
        return self

    def apply(self, F) : 
        """ Apply x,y -> x, F(x,y) to all datas"""
        plt.subplot(*self.coord)
        for line in self.gca().get_lines():
            X,Y = line.get_data()
            line.set_data(*F(X,Y))

    def clear(self):
        """ Clear the graph """
        plt.subplot(*self.coord).clear()

    def plot(self, *args, **kwargs) :
        """ 
        Focuss on this graph and call the plot function
        """
        plt.subplot(*self.coord)
        plot(*args, **kwargs)

#############################

class Plotter(object) :
    """ 
       A quick interface to a multiple matplotlib figure.

       Each graph is an instance of a *Plotter_OneGraph* and can be accessed with the brackets: 
       
       * self [i] is the i-th graph. 
       * self [i,j] is the graph on row i, col j.
     """

    def __init__(self, MultiGraph=(1,1)) :
        """
        Starts the plotter and open a window with the figure.
        :param MultiGraph: a tuple (n,m) defining how many graphs are in the plot (n rows, m columns).
       """
        for method in ['savefig','xlim','ylim','grid','text','draw','setp']: setattr(self,method, getattr(plt,method))
        self.__Commands = []
        plt.figure(1, figsize = (12,8))
        self.plt = plt
        self._grid = MultiGraph
        self.__graphlist = []
        n=1
        for i in range (self._grid[0]):
            for j in range (self._grid[1]):
                self.__graphlist.append(Plotter_OneGraph ( self._grid[0],self._grid[1], n))
                n += 1
        plt.figure(1).show()

    def __getitem__( self, item ):
        """Access a specific graph.  Can use either p[num] or p[row, col]."""
        def convert(i,j):
            if i >= self._grid[0] or j >= self._grid[1]:
                raise IndexError, 'graph index out of range'
            return i*self._grid[1] + j

        if type(item) == type(1):
            return self.__graphlist[item]._focuss()
        elif type(item) == type( () ) and len(item) <= 2:
            return self.__graphlist[convert(item[0],item[1])]._focuss()
        else:
            raise TypeError, 'graph index must be integer or two integers'


