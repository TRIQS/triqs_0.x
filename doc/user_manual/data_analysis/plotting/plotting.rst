
.. index:: matplotlib plotter

.. module:: pytriqs.plot

.. _plotting:

Plotting TRIQS objects
################################

TRIQS objects can be easily plotted, for example with the standard python plot toolkit `matplotlib
<http://matplotlib.sourceforge.net/>`_.
In fact, TRIQS introduces a simple :ref:`plot protocol <plot_protocol>`, which allows to plot 
objects which have a graphical representation.

A thin layer above matplotlib
=================================

TRIQS defines a function *oplot*, similar to the standard matplotlib pyplot.plot function,
but that can plot TRIQS objects (in fact *any* object, see below).

We can reproduce the first example of the Green function tutorial :

.. plot:: green/example.py
   :include-source:
   :scale: 70

The *oplot* function takes :

* as arguments any object that implements the :ref:`plot protocol <plot_protocol>`, 
  for example Green function, Density of state : in fact any object where plotting is reasonable and has been defined ...

* string formats following objects, as in regular matplotlib, like in the example above.

* regular options of the matplotlib *pyplot.plot* function 

* options specific to the object to be plotted : here the `x_window` tells the Green function to plot itself in a reduced window of :math:`\omega_n`.
  
Multiple panels figures
=================================

`Only valid for matplotlib v>=1.0`.

While one can use the regular matplotlib subfigure to make multi-panel figures, 
subplots makes it a bit more pythonic :

.. plot:: data_analysis/plotting/example.py
   :include-source:
   :scale: 70

.. index:: plotting protocol

.. _plot_protocol:

Plot protocol [Advanced]
===========================

What do we need to implement to plot an object ? 
Simply a little `_plot_` method that reduces the object to a set of curves.
This section describes the conventions on this function.

As soon as an object defines this method, it can be plotted by the `oplot` function of `pytriqs.plot.mpl_interface`.
See example below.

.. function:: _plot_( OptionsDict )

  * OptionDict is a dictionnary of options.

  .. warning:: 
     * The method _plot_ must consume the options it uses (using  e.g. the pop method of dict).
     * Other options will be passed to matplotlib, so leaving spurious options here will lead to errors.

  :rtype: a list of dict representing one curve each. These dict must have the following fields:

    * *xdata* : A 1-dimensional numpy array describing the x-axis points
    * *ydata* : A 1-dimensional numpy array describing the y-axis points
    * *label* : Label of the curve for the legend of the graph
    * *type* : a string : currently "XY" [ optional] 

   and optionally : 
    
    * *xlabel* : a label for the x axis. The last object plotted will overrule the previous ones.
    * *ylabel* : a label for the y axis. The last object plotted will overrule the previous ones.

Example
-------

Here's a simple example to illustrate the protocol:

.. plot:: data_analysis/plotting/myobject.py
   :include-source:
   :scale: 70

Example with options
---------------------------

A little bit more complex, with options. 
Note the use of the `pop method of dict <http://docs.python.org/library/stdtypes.html#dict>`_, 
which returns and removes the entry from the dict (with a default value).

.. plot:: data_analysis/plotting/myobject2.py
   :include-source:
   :scale: 70


