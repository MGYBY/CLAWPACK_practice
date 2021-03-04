from matplotlib import animation
from IPython.display import HTML
import numpy as np
import matplotlib.pyplot as plt

def plot_waves(claw,ylim=(0,1.2),save_plots=False):
    fig = plt.figure(figsize=[12,4])
    ax1 = fig.add_subplot(111)
    fills = []
    frame = claw.frames[0]
    b = frame.aux[0,:,:]
    h = frame.q[0,:,:]
    surface = np.maximum(b,h+b)

    x, y = frame.state.grid.p_centers    
    save_plots = True
    slice = 1
    #line, = ax1.plot(x[:,0],surface[:,slice],'-k',linewidth=3)
    fill = ax1.fill_between(x[:,0],b[:,slice],surface[:,slice],facecolor='blue')
    fill2 = ax1.fill_between(x[:,0],0*b[:,slice],b[:,slice],facecolor='brown')
    fills = [fill,fill2]
    ax1.set_xlim(-15,15)
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