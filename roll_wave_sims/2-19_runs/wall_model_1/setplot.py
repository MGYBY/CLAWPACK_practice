
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

from __future__ import absolute_import

import numpy as np
from matplotlib import pylab

hn = 0.00798
domain_x = 42.00
domain_y = 2.0
x1 = 40.0
y1 = -0.20

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
    def change_fonts(current_data):
        pylab.xticks(fontsize=21, fontname="Tex Gyre Pagella")
        pylab.yticks(fontsize=21, fontname="Tex Gyre Pagella")
        "Fill the step area with a black rectangle."
        import matplotlib.pyplot as plt
        rectangle = plt.Rectangle((x1,y1),0.4,0.4,color="k",fill=True)
        plt.gca().add_patch(rectangle)        
        
    def change_fonts2(current_data):
        pylab.xticks(fontsize=17, fontname="Tex Gyre Pagella")
        pylab.yticks(fontsize=17, fontname="Tex Gyre Pagella")   
        "Fill the step area with a black rectangle."
        import matplotlib.pyplot as plt
        rectangle = plt.Rectangle((x1,y1),0.4,0.4,color="k",fill=True)
        plt.gca().add_patch(rectangle)  
        
    plotfigure = plotdata.new_plotfigure(name='q[0]', figno=0)
    plotfigure.kwargs = {'figsize':[10,10],'facecolor':'white'}

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
    plotaxes.afteraxes = change_fonts
    
    #-----------------------------------------
    # Figure for zoomed-in pcolor plot
    #-----------------------------------------    
    plotfigure = plotdata.new_plotfigure(name='q[0]_zoomed', figno=1)
    plotfigure.kwargs = {'figsize':[10,10],'facecolor':'white'}
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [39.4,domain_x]
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
    plotaxes.afteraxes = change_fonts2
    
    #-----------------------------------------
    # Figure for momentum pcolor plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='q[1]', figno=2)
    plotfigure.kwargs = {'figsize':[10,10],'facecolor':'white'}

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
    plotaxes.afteraxes = change_fonts2
    
    #-----------------------------------------
    # Figure for cross section at y=0
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='cross-section', figno=3)
    plotfigure.kwargs = {'figsize':[10,10],'facecolor':'white'}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0.0,domain_x]
    plotaxes.ylimits = [0.0,hn*9.0]
    plotaxes.title = 'Cross section at y=0'
    def plot_topo_xsec(current_data):
        from pylab import plot, cos,sin,where,legend,nan
        t = current_data.t

        x = np.linspace(0.0,domain_x,201)
        #y = 0.
        B = where(x>40.0, where(x<40.40,7.0,0.0), 0.0)
        plot(x, B, 'g', label="bathymetry")
        legend()
        
        pylab.xticks(fontsize=17, fontname="Tex Gyre Pagella")
        pylab.yticks(fontsize=17, fontname="Tex Gyre Pagella")
        
    plotaxes.afteraxes = plot_topo_xsec

    plotitem = plotaxes.new_plotitem(plot_type='1d_from_2d_data')

    def xsec(current_data):
        # Return x value and surface depth at this point, along y=0
        from pylab import where,ravel
        x = current_data.x
        y = ravel(current_data.y)
        dy = current_data.dy
        q = current_data.q

        ij = where((y <= dy/1.) & (y > -dy/1.))
        x_slice = ravel(x)[ij]
        ij1 = where((x_slice > 40.40) | (x_slice < 40.0))
        x_slice = x_slice[ij1]
        depth_slice = ravel(q[0,:,:])[ij]
        depth_slice = depth_slice[ij1]
        return x_slice, depth_slice

    plotitem.map_2d_to_1d = xsec
    plotitem.plotstyle = 'k-o'     ## need to be able to set amr_plotstyle
    plotitem.kwargs = {'markersize':4}

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

    
