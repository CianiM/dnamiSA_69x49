import numpy as np
import matplotlib.pyplot as plt
from pltdnami import loadrstax
from plot_result import build_q


def main(fbeg,fend,nstep,bc,scale=0,scale2=0):
    field = { 'rho':1,
              'u'  :2,
              'v'  :3,
              'et' :4,
              'nut':5,
              'tau':6}
    var=field.tau
    with open('x_coord.dat','r') as f:
        x_axis = []
        coord=f.readlines()
        x_axis = [ float(s) for s in coord]
        x_axis = np.array(x_axis)

    q = build_q('./output/','outputoutput',scale,bc)
    tau_wall=q[:,:,0,var]
    fig = plt.figure(figsize=(10,5))
    plt.plot(tau_wall,x_axis)
    plt.show()