
"""
Module to create topo and qinit data files for this example.
"""

from __future__ import absolute_import
from clawpack.geoclaw.topotools import Topography
from numpy import *

#from pyclaw.data import Data
#probdata = Data('setprob.data')

So = 0.05011
hn = 0.00798
un = 1.03774
grav = 9.81
cf = grav * So * 2 * hn / (un**2)

domain_x = 42.0
domain_y = 2.0

def maketopo():
    """
    Output topography file for the entire domain
    """
    nxpoints=2520
    nypoints=150
    xupper=domain_x
    yupper=domain_y/2.0
    xlower = 0.e0
    ylower = -domain_y/2.0
    outfile= "channel.topotype2"

    topography = Topography(topo_func=topo)
    topography.x = linspace(xlower,xupper,nxpoints)
    topography.y = linspace(ylower,yupper,nypoints)
    topography.write(outfile, topo_type=2, Z_format="%22.15e")


def topo(x,y):
    """
    simple channel
    """
    z = -So*x+where((y>-0.20), where(y<0.20, where(x<40.40, where(x>40.0, 4.0, 0.0), 0.0), 0.0), 0.0)
    return z


if __name__=='__main__':
    maketopo()
