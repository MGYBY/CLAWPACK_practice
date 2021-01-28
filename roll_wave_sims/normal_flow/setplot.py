
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

from __future__ import absolute_import
import numpy

So = 0.05011
hn = 0.00798
un = 1.03774
grav = 9.81
cf = grav * So * 2 * hn / (un**2)

domain_x = 40.0
domain_y = 1.0
centerline_index = 8
cmax1 = hn*3.0

#--------------------------
def setplot(plotdata=None):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    from clawpack.visclaw import colormaps, geoplot

    if plotdata is None:
        from clawpack.visclaw.data import ClawPlotData
        plotdata = ClawPlotData()


    plotdata.clearfigures()  # clear any old figures,axes,items data

    def set_drytol(current_data):
        # The drytol parameter is used in masking land and water and
        # affects what color map is used for cells with small water depth h.
        # The cell will be plotted as dry if h < drytol.
        # The best value to use often depends on the application and can
        # be set here (measured in meters):
        current_data.user["drytol"] = 1.e-5

    plotdata.beforeframe = set_drytol

    #-----------------------------------------
    # Figure for pcolor plot
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='pcolor', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes('pcolor')
    plotaxes.title = 'Depth'
    plotaxes.scaled = True
    
    def max_cmap(current_data):
        q = current_data.q
        return ((q[1,:,:].max())*1.1)

    # Water
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.depth
    plotitem.pcolor_cmap = geoplot.tsunami_colormap
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = cmax1
    plotitem.add_colorbar = True
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 1

    # Land
    plotitem = plotaxes.new_plotitem(plot_type='2d_pcolor')
    plotitem.plot_var = geoplot.land
    plotitem.pcolor_cmap = geoplot.land_colors
    plotitem.pcolor_cmin = 0.0
    plotitem.pcolor_cmax = 100.0
    plotitem.add_colorbar = False
    plotitem.amr_celledges_show = [0,0,0]
    plotitem.patchedges_show = 1
    plotaxes.xlimits = [0.0,domain_x]
    plotaxes.ylimits = [-domain_y/2.0,domain_y/2.0]

    # # Add contour lines of bathymetry:
    # plotitem = plotaxes.new_plotitem(plot_type='2d_contour')
    # plotitem.plot_var = geoplot.topo
    # from numpy import arange, linspace
    # plotitem.contour_levels = linspace(-.1, 0.5, 20)
    # plotitem.amr_contour_colors = ['k']  # color on each level
    # plotitem.kwargs = {'linestyles':'solid'}
    # plotitem.amr_contour_show = [1]  
    # plotitem.celledges_show = 0
    # plotitem.patchedges_show = 0
    # plotitem.show = True

    #-----------------------------------------
    # Figure for centerline slice
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='Centerline_slice', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0.0,domain_x]
    plotaxes.ylimits = [0.0,hn*1.5]
    plotaxes.title = 'Centerline slice'
    #def depth_plot(current_data):
        #from pylab import plot, cos,sin,where,legend,nan
        #t = current_data.t
        #q = current_data.q
        #q1 = q[1,:,centerline_index]
        #x = numpy.linspace(0,domain_x,len(q1))

        #plot(x, q1, 'k.', label="true solution", linewidth=2)
        
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 0
    plotitem.plotstyle = 'k--'     ## need to be able to set amr_plotstyle
    #plotitem.kwargs = {'markersize':3}
    #plotitem.amr_show = [1]  # plot on all levels
    #plotaxes.afteraxes = depth_plot


    #-----------------------------------------
    # Figure for grids alone
    #-----------------------------------------
    plotfigure = plotdata.new_plotfigure(name='grids', figno=2)
    plotfigure.show = True

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [0.0,domain_x]
    plotaxes.ylimits = [-domain_y/2.0,domain_y/2.0]
    plotaxes.title = 'grids'
    plotaxes.scaled = True

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='2d_patch')
    plotitem.amr_patch_bgcolor = ['#ffeeee', '#eeeeff', '#eeffee']
    plotitem.amr_celledges_show = [1,1,0]   
    plotitem.amr_patchedges_show = [1]


    #-----------------------------------------
    
    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_gaugenos = []             # list of gauges to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'   # pointer for top of index
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?
    plotdata.parallel = True                 # make multiple frame png's at once

    return plotdata

    
