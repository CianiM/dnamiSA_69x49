import numpy as np
import matplotlib.pyplot as plt
from pltdnami import loadrstax
from plot_result import build_q
import sys
import scipy.integrate as integrate

def main(fbeg,fend,nstep,bc,scale=0,scale2=0):
    axis  = build_q('./grid/','output',0,bc)
    q = build_q('./output/','outputoutput',scale,bc)
    ksi = axis[:,:,0,0]
    x0=12
    field = { 'rho':1,
              'u'  :2,
              'v'  :3,
              'et' :4,
              'nut':5,
              'tau':6}
    var=field['tau']-1
    x=ksi[x0:,0]
    dx=np.zeros(len(x)-1)
    integ=dx
    tau_wall=q[x0:,0,0,var]
    for i in np.arange(0,len(x)-1):
        dx[i]=x[i+1]-x[i]
        integ[i]=dx[i]*(tau_wall[i+1]+tau_wall[i])/2

    cf=2*2*np.sum(integ)
#    with open('x_coord.dat','r') as f:
#        x_axis = []
#        coord=f.readlines()
#        x_axis = [ float(s) for s in coord]
#        x_axis = np.array(x_axis)
    print(cf)
    fig = plt.figure(figsize=(10,5))
    ax = fig.add_subplot(111)
    im = ax.plot(x,tau_wall)
    plt.show();sys.exit()


if __name__ == '__main__':
    name,fbeg,fend,nstep,bc,scale,scale2 = sys.argv
    if bc=="True":
        bc = True
    else:
        bc = False    
    main(np.int(fbeg),np.int(fend),np.int(nstep),bc,np.int(scale),np.int(scale2))