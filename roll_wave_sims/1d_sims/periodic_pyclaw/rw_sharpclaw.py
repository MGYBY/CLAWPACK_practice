#!/usr/bin/env python
# encoding: utf-8

r"""
Shallow water flow
==================

Solve the one-dimensional shallow water equations:

.. math::
    h_t + (hu)_x & = 0 \\
    (hu)_t + (hu^2 + \frac{1}{2}gh^2)_x & = -gb_xh-\frac{1}{2}c_fu^2

Roll wave simulation
Periodic BC for now
"""

from __future__ import absolute_import
import numpy as np
from clawpack import riemann
from clawpack.riemann.shallow_roe_with_efix_1D_constants import depth, momentum, num_eqn
import os

# to mimic Brock's setup
# parameter list
xlower = 0.0
xupper = 1.285
mx = 300

normal_depth = 0.00798
normal_velocity = 1.038
dist_amp = 0.01 # for two kinds of stability problems
wave_length = xupper
sim_time = 24.0
output_interval = 1.0 # in seconds
output_flag = 1.0 #avoid ouput overwrite
channel_slope = 0.05011
gravity_val = 9.81
cf = gravity_val * channel_slope * 2 * normal_depth / (normal_velocity**2)



def step_swe(solver,state,dt):
    """
    source terms for channel slope and bed-friction for classic solvers
    RK2 TVD time integration
    """

    q = state.q

    qstar = np.empty(q.shape) # note that only hu component is used

    # qstar[0,:,:] = q[0,:,:] - dt2/rad * q[2,:,:]
    qstar[1,:] = q[1,:] + dt*(channel_slope*gravity_val*q[0,:]-cf/2.0*(q[1,:])**2/((q[0,:])**2))

    # q[0,:,:] = q[0,:,:] - dt/rad * qstar[2,:,:]
    q[1,:] = 0.50*q[1,:]+0.50*qstar[1,:]+0.50*dt*(channel_slope*gravity_val*q[0,:]-cf/2.0*(qstar[1,:])**2/((q[0,:])**2))


def dq_swe(solver,state,dt):
    """
    source terms for channel slope and bed-friction for WENO solvers
    forward Euler
    RK2 TVD time integration
    """
    q   = state.q

    dq = np.empty(q.shape)

    dq[0,:] = 0.0
    dq[1,:] = dt*(channel_slope*gravity_val*q[0,:]-cf/2.0*(q[1,:])**2/((q[0,:])**2))

    return dq

def b4step(solver,state):
    """
    1. output XYZ files
    2. obtain maximum values
    """
    q  = state.q
    height = q[0,:]
    mom = q[1,:]
    t = state.t
    dt = solver.dt
    # X = state.grid.p_nodes
    dx = state.grid.delta[0]
    # output maximum values every time step
    max_h = np.max(height)
    max_h_ind = np.argmax(height)
    f1 = open("maxHU.txt", "a+")
    f1.write("%18.8e %18.8e %18.8e %18.8e \n" % (t, (dx*(max_h_ind-0.50)), max_h, mom[max_h_ind]))

    # output XYZ files every output_interval
    if t%output_interval<dt:
        format_string_time = f"{t:.1f}"
        file_name = 'outXYZ_%s' % format_string_time
        # check the existence of the file
        if (not os.path.exists('./%s' % file_name)):
            f = open(file_name, "w+")
            for k in range(len(height)):
                # f.write("%18.8e %18.8e %18.8e \n" % (X[k], height[k], mom[k]))
                f.write("%18.8e %18.8e %18.8e \n" % ((dx*(k+0.50)), height[k], mom[k]))
            print("Writing XYZ output for %f" % t)




def setup(use_petsc=False,kernel_language='Fortran',outdir='./_output',solver_type='sharpclaw',
          riemann_solver='roe', disable_output=False):
    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if kernel_language == 'Python':
        if riemann_solver.lower() == 'roe':
            raise Exception('Python Roe solver not implemented.')
        elif riemann_solver.lower() == 'hlle':
            rs = riemann.shallow_1D_py.shallow_hll_1D
    elif kernel_language == 'Fortran':
        if riemann_solver.lower() == 'roe':
            rs = riemann.shallow_roe_with_efix_1D
        elif riemann_solver.lower() == 'hlle':
            rs = riemann.shallow_hlle_1D
 
    if solver_type == 'classic':
        solver = pyclaw.ClawSolver1D(rs)
        # solver.limiters = pyclaw.limiters.tvd.vanleer
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver1D(rs)    
    if use_petsc:
        import clawpack.petclaw as pyclaw
    else:
        from clawpack import pyclaw

    if solver_type=='sharpclaw':
        solver.dq_src = dq_swe
        # solver.call_before_step_each_stage = True
        solver.weno_order = 5
        solver.lim_type   = 2 # weno resonstruction
        solver.cfl_max = 0.36
        solver.cfl_desired = 0.35
    else:
        # solver = pyclaw.ClawSolver2D(riemann.euler_5wave_2D)
        solver.step_source = step_swe
        solver.source_split = 1 # Godunov splitting
        # solver.limiters = [11, 11] # 11 for A-R limiter
        solver.limiters = [4, 4] # 4 for MC limiter
        solver.cfl_max = 0.36
        solver.cfl_desired = 0.35



    solver.kernel_language = kernel_language

    solver.bc_lower[0] = pyclaw.BC.periodic
    solver.bc_upper[0] = pyclaw.BC.periodic

    solver.before_step = b4step


    x = pyclaw.Dimension(xlower,xupper,mx,name='x')
    domain = pyclaw.Domain(x)
    state = pyclaw.State(domain,num_eqn)

    # Gravitational constant
    state.problem_data['grav'] = 9.81
    state.problem_data['dry_tolerance'] = 1e-5
    state.problem_data['sea_level'] = 0.0
    
    xc = state.grid.x.centers


    # I.C.: spatially varying disturbance
    state.q[depth,:] = normal_depth*(1.0+dist_amp*np.sin(2.0*np.pi*xc/wave_length))
    state.q[momentum,:] = normal_velocity*normal_depth*(1.0+dist_amp*np.sin(2.0*np.pi*xc/wave_length))

    claw = pyclaw.Controller()
    claw.keep_copy = True
    if disable_output:
        claw.output_format = None
    claw.output_style = 1
    claw.tfinal = sim_time
    claw.num_output_times = int(sim_time/output_interval) # conversion between two output styles
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir
    claw.setplot = setplot

    return claw


# #--------------------------
# def setplot(plotdata):
# #--------------------------
#     """ 
#     Specify what is to be plotted at each frame.
#     Input:  plotdata, an instance of visclaw.data.ClawPlotData.
#     Output: a modified version of plotdata.
#     """ 
#     plotdata.clearfigures()  # clear any old figures,axes,items data

#     # Figure for depth
#     plotfigure = plotdata.new_plotfigure(name='Water height', figno=0)

#     # Set up for axes in this figure:
#     plotaxes = plotfigure.new_plotaxes()
#     # plotaxes.xlimits = [-5.0,5.0]
#     plotaxes.title = 'Water height'
#     plotaxes.axescmd = 'subplot(211)'

#     # Set up for item on these axes:
#     plotitem = plotaxes.new_plotitem(plot_type='1d')
#     plotitem.plot_var = depth
#     plotitem.plotstyle = '-'
#     plotitem.color = 'b'
#     plotitem.kwargs = {'linewidth':3}

#     # Figure for momentum[1]
#     #plotfigure = plotdata.new_plotfigure(name='Momentum', figno=1)

#     # Set up for axes in this figure:
#     plotaxes = plotfigure.new_plotaxes()
#     plotaxes.axescmd = 'subplot(212)'
#     # plotaxes.xlimits = [-5.0,5.0]
#     plotaxes.title = 'Momentum'

#     # Set up for item on these axes:
#     plotitem = plotaxes.new_plotitem(plot_type='1d')
#     plotitem.plot_var = momentum
#     plotitem.plotstyle = '-'
#     plotitem.color = 'b'
#     plotitem.kwargs = {'linewidth':3}
    
#     return plotdata


#--------------------------
def setplot(plotdata):
#--------------------------
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    """ 
    from matplotlib import animation, pylab
    from IPython.display import HTML
    from clawpack.visclaw.data import ClawPlotData
    from clawpack.visclaw import plotpages
    def change_fonts(current_data):
        pylab.xticks(fontsize=15)
        pylab.yticks(fontsize=15)     

    # Figure for depth and momentum:

    plotdata.clearfigures()  # clear any old figures,axes,items data

    # Figure for depth
    plotfigure = plotdata.new_plotfigure(name='Water height', figno=0)
    plotfigure.kwargs = {'figsize':[9,12],'facecolor':'white'}

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    # plotaxes.xlimits = [-5.0,5.0]
    plotaxes.title = 'Water height'
    plotaxes.axescmd = 'subplot(211)'
    plotaxes.afteraxes = change_fonts

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = depth
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.kwargs = {'linewidth':3}

    # Figure for momentum[1]
    #plotfigure = plotdata.new_plotfigure(name='Momentum', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(212)'
    # plotaxes.xlimits = [-5.0,5.0]
    plotaxes.title = 'Momentum'
    plotaxes.afteraxes = change_fonts

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = momentum
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.kwargs = {'linewidth':3}
    
    plotdata.printfigs = True                # print figures False to supress png output
    plotdata.plotdir = './_plots'
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = False                     # create html files of plots?
    plotdata.html_homelink = '../README.html'
    plotdata.latex = False                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?
    plotpages.plotclaw_driver(plotdata)

    
    return plotdata


if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)
