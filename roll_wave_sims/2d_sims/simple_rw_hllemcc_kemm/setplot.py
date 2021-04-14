
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
"""

from __future__ import absolute_import

import numpy as np
from matplotlib import pylab
import pandas

hn = 0.00798
domain_x = 40.0

# --------------------------


def setplot(plotdata=None):
    # --------------------------
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of clawpack.visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.

    """

    if plotdata is None:
        from clawpack.visclaw.data import ClawPlotData
        plotdata = ClawPlotData()

    from clawpack.visclaw import colormaps

    plotdata.clearfigures()  # clear any old figures,axes,items data

    # Figure for pcolor plot
    plotfigure = plotdata.new_plotfigure(name='q[0]', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'Depth Contour'
    plotaxes.scaled = False

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 0
    plotitem.pcolor_cmap = colormaps.red_yellow_blue
    plotitem.pcolor_cmin = 0.00
    plotitem.pcolor_cmax = 2.80*hn
    plotitem.add_colorbar = True
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 1
    plotitem.show = True       # show on plot?

    # -----------------------------------------
    # Figure for cross section at y=0
    # -----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='cross-section', figno=3)
    plotfigure.kwargs = {'figsize': [10, 10], 'facecolor': 'white'}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0.0, domain_x]
    plotaxes.ylimits = [0.0, hn*14.0]
    plotaxes.title = 'Cross section at y=0'

    def plot_topo_xsec(current_data):
        from pylab import plot, cos, sin, where, legend, nan
        t = current_data.t

        x = np.linspace(0.0, domain_x, 201)
        #y = 0.
        B = where(x > 40.0, where(x < 40.40, 7.0, 0.0), 0.0)
        plot(x, B, 'g', label="internal walls")
        legend()
        pylab.legend(fontsize=18)

        pylab.xticks(fontsize=18, fontname="Tex Gyre Pagella")
        pylab.yticks(fontsize=18, fontname="Tex Gyre Pagella")
        t = current_data.t
        pylab.title("Run-up at time t = %10.4e" % t, fontsize=16)

    plotaxes.afteraxes = plot_topo_xsec

    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')

    def xsec(current_data):
        # Return x value and surface depth at this point, along y=0
        from pylab import where, ravel
        x = current_data.x
        y = ravel(current_data.y)
        dy = current_data.dy
        q = current_data.q
        frameno= str(current_data.frameno)

        ij = where((y <= dy/1.) & (y > -dy/1.))
        x_slice = ravel(x)[ij]
        ij1 = where((x_slice > 40.40) | (x_slice < 40.0))
        x_slice = x_slice[ij1]
        depth_slice = ravel(q[0, :, :])[ij]
        depth_slice = depth_slice[ij1]
        profile_arr = np.array([np.transpose(x_slice), np.transpose(depth_slice)])
        data_frame = pandas.DataFrame(profile_arr.T)
        data_frame.to_csv(frameno+'_profile.csv') 
        return x_slice, depth_slice

    plotitem.map_2d_to_1d = xsec
    plotitem.plotstyle = 'k-o'  # need to be able to set amr_plotstyle
    plotitem.kwargs = {'markersize': 4}
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via clawpack.visclaw.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    # plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata
