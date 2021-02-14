
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

from __future__ import absolute_import

import numpy as np

hn = 0.00798
domain_x = 10.00
domain_y = 4.0

#--------------------------
def setplot(plotdata=None):
#--------------------------
    
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
    

    #-----------------------------------------
    # Figure for pcolor plot
    #-----------------------------------------
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
    plotitem.pcolor_cmap = colormaps.yellow_red_blue # not the default colormap
    plotitem.pcolor_cmin = 0.00
    plotitem.pcolor_cmax = 2.80*hn
    plotitem.add_colorbar = True
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 1
    plotitem.show = True       # show on plot?
    
    #-----------------------------------------
    # Figure for zoomed-in pcolor plot
    #-----------------------------------------    
    plotfigure = plotdata.new_plotfigure(name='q[0]_zoomed', figno=1)
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [8.60,domain_x]
    plotaxes.ylimits = [-domain_y/2.0,domain_y/2.0]
    plotaxes.title = 'Zoomed-in Depth Contour'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 0
    plotitem.pcolor_cmap = colormaps.yellow_red_blue # not the default colormap
    plotitem.pcolor_cmin = 0.00
    plotitem.pcolor_cmax = 5.00*hn
    plotitem.add_colorbar = True
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 1
    plotitem.show = True       # show on plot?    
    
    #-----------------------------------------
    # Figure for momentum pcolor plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='q[1]', figno=2)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'
    plotaxes.title = 'qx Contour'
    plotaxes.scaled = False

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = 1
    plotitem.pcolor_cmap = colormaps.yellow_red_blue # not the default colormap
    plotitem.pcolor_cmin = 0.00
    plotitem.pcolor_cmax = 'auto'
    plotitem.add_colorbar = True
    plotitem.celledges_show = 0
    plotitem.patchedges_show = 1
    plotitem.show = True       # show on plot?

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

    
