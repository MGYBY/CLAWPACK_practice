# %matplotlib inline
import numpy as np
import matplotlib.pyplot as plt
def wave_maker_bc(state,dim,t,qbc,auxbc,num_ghost):
    "Generate waves at left boundary as if there were a moving wall there."
    if dim.on_lower_boundary:
        qbc[0,:num_ghost,:]=qbc[0,num_ghost,:] 
        t=state.t;
        amp = state.problem_data['amp'];
        if t<=state.problem_data['t1']: 
            vwall = amp*(np.sin(t*np.pi/1.5))
        else: 
            vwall=0.
        for ibc in range(num_ghost-1):
            qbc[1,num_ghost-ibc-1,:] = 2*vwall - qbc[1,num_ghost+ibc,:]

def qinit(state):
    "Gaussian surface perturbation"
    x0=0.
    y0=0.

    b = state.aux[0,:,:] # Bathymetry

    X,Y = state.grid.p_centers
    xleft = X.min()
    surface = ambient_surface_height+pulse_amplitude*np.exp(-(X-(xleft+2.))**2/pulse_width)
    state.q[0,:,:] = surface - b
    state.q[1,:,:] = 0.
    state.q[2,:,:] = 0.
    
def bathymetry(x):
    "Flat bottom for x<3; then a steep slope to x=5, followed by a gentle slope."
    return (-x*0.30)

def bed_friction(solver,state,dt):
    dt2 = dt/2.
    q = state.q
    qstar = np.empty(q.shape)
    
    qstar[1,:,:] = q[1,:,:] - dt2*cf/2 * q[1,:,:]*np.sqrt(q[1,:,:]**2+q[2,:,:]**2)/q[0,:,:]**2
    qstar[2,:,:] = q[2,:,:] - dt2*cf/2 * q[2,:,:]*np.sqrt(q[1,:,:]**2+q[2,:,:]**2)/q[0,:,:]**2
    
    q[1,:,:] = q[1,:,:] - dt*cf/2 * qstar[1,:,:]*np.sqrt(qstar[1,:,:]**2+qstar[2,:,:]**2)/q[0,:,:]**2
    q[2,:,:] = q[2,:,:] - dt*cf/2 * qstar[2,:,:]*np.sqrt(qstar[1,:,:]**2+qstar[2,:,:]**2)/q[0,:,:]**2

        
def setup(num_cells=500,tfinal=30,solver_type='classic',num_output_times=150):

    from clawpack import riemann
    from clawpack import pyclaw

    if solver_type == 'classic':
        solver = pyclaw.ClawSolver2D(riemann.sw_aug_2D)
        solver.step_source = bed_friction
        solver.dimensional_split=True
        solver.limiters = pyclaw.limiters.tvd.minmod
        solver.cfl_max     = 0.45
        solver.cfl_desired = 0.3
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver2D(riemann.sw_aug_2D)

    solver.bc_lower[0] = pyclaw.BC.custom 
    solver.user_bc_lower = wave_maker_bc
    solver.bc_upper[0] = pyclaw.BC.extrap
    solver.bc_lower[1] = pyclaw.BC.periodic
    solver.bc_upper[1] = pyclaw.BC.periodic

    solver.aux_bc_lower[0] = pyclaw.BC.extrap
    solver.aux_bc_upper[0] = pyclaw.BC.extrap
    solver.aux_bc_lower[1] = pyclaw.BC.periodic
    solver.aux_bc_upper[1] = pyclaw.BC.periodic

    solver.fwave = True

    # Domain:
    xlower = -0.0;  xupper =  40.
    ylower = -0.5;  yupper =  0.5

    mx = num_cells
    my = 2

    x = pyclaw.Dimension(xlower,xupper,mx,name='x')
    y = pyclaw.Dimension(ylower,yupper,my,name='y')
    domain = pyclaw.Domain([x,y])

    num_aux = 1
    state = pyclaw.State(domain,solver.num_eqn,num_aux)
    state.aux[:,:,:] = bathymetry(state.p_centers[0])

    state.problem_data['grav'] = 9.810   # Gravitational force
    state.problem_data['t1']   = 50.0  # Stop generating waves after this time
    state.problem_data['amp']  = 0.3   # Amplitude of incoming waves
    qinit(state)

    #===========================================================================
    # Set up controller and controller parameters
    #===========================================================================
    claw = pyclaw.Controller()
    claw.tfinal = tfinal
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.num_output_times = num_output_times
    claw.keep_copy = True
    claw.output_format = None

    return claw

ambient_surface_height  = 0.02
pulse_amplitude         = 0. # Use this to add an initial Gaussian wave
pulse_width             = 1.
cf=0.003

claw = setup(num_cells=1000,tfinal=30.)

#claw.verbosity=0 # Use this to suppress output during the run
claw.run()

#===========================================================================
# Plotting Results
#===========================================================================
from matplotlib import animation
from IPython.display import HTML

def plot_waves(claw,ylim=(0,1.2),save_plots=True):
    fig = plt.figure(figsize=[12,4])
    ax1 = fig.add_subplot(111)
    fills = []
    frame = claw.frames[0]
    b = frame.aux[0,:,:]
    h = frame.q[0,:,:]
    surface = np.maximum(b,h+b)

    x, y = frame.state.grid.p_centers    
    # save_plots = True
    slice = 1
    #line, = ax1.plot(x[:,0],surface[:,slice],'-k',linewidth=3)
    fill = ax1.fill_between(x[:,0],b[:,slice],surface[:,slice],facecolor='blue')
    fill2 = ax1.fill_between(x[:,0],0*b[:,slice],b[:,slice],facecolor='brown')
    fills = [fill,fill2]
    ax1.set_xlim(0.0,40.0)
    if ylim: ax1.set_ylim(ylim)

    def fplot(frame_number):
        fills[-2].remove()
        fills[-1].remove()
        frame = claw.frames[frame_number]
        b = frame.aux[0,:,:]
        h = frame.q[0,:,:]
        surface = np.maximum(b,h+b)
        #line.set_data(x[:,0],surface[:,slice])
        fill = ax1.fill_between(x[:,0],b[:,slice],surface[:,slice],facecolor='blue',where=b[:,slice]<surface[:,slice])
        fill2 = ax1.fill_between(x[:,0],0*b[:,slice],b[:,slice],facecolor='brown')
        fills.append(fill)
        fills.append(fill2)
        if save_plots:
            fname = 'frame'+str(frame_number).zfill(4)+'.eps'
            fig.savefig(fname)   
        return fill,

    anim = animation.FuncAnimation(fig, fplot, frames=len(claw.frames), interval=100, repeat=False)
    plt.close()
    return HTML(anim.to_jshtml())

plot_waves(claw, ylim=(0,0.05))
