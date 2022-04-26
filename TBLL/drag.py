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
    eta = axis[:,:,0,1]
    y=eta[0,:]
    x0=12
    field = { 'rho':1,
              'u'  :2,
              'v'  :3,
              'et' :4,
              'nut':5,
              'tau':6}
    var_rho=field['rho']-1
    var_u=field['u']-1
    var_nut=field['nut']-1
    var_tau=field['tau']-1
    x=ksi[x0:,0]
    dx=np.zeros(len(x)-1)
    integ=dx
    rho = q[56,:,0,var_rho]
    tau_wall=q[x0:,0,0,var_tau]
    u=q[56,:,0,var_u]
    tw = q[56,:,0,var_tau]
    uw=np.sqrt(tw/rho)
    up = u/uw

    for i in np.arange(0,len(x)-1):
        dx[i]=x[i+1]-x[i]
        integ[i]=dx[i]*(tau_wall[i+1]+tau_wall[i])/2
    
    cf = 2*q[56,0,0,var_u]
    cf_tot=2*2*np.sum(integ)
#    with open('x_coord.dat','r') as f:
#        x_axis = []
#        coord=f.readlines()
#        x_axis = [ float(s) for s in coord]
#        x_axis = np.array(x_axis)
    with open('flatplate_u.dat','r') as f:
        dat=f.readlines()[2:387]
        data = [word for line in dat for word in line.split()] 
        cc = [ float(s) for s in data]
        cc = np.array(cc)
        cc = cc.reshape((385,2))
    u_nasa = cc[:,0]
    y_nasa = cc[:,1]
    print( cf )
    print(cf_tot)
    fig1 = plt.figure(figsize=(5,5))
    ax = fig1.add_subplot(111)
    im = ax.plot(u_nasa,y_nasa,'r')
    ax.plot(u,y,'k--')
    #fig2 = plt.figure(figsize=(5,5))
    plt.show();sys.exit()


if __name__ == '__main__':
    name,fbeg,fend,nstep,bc,scale,scale2 = sys.argv
    if bc=="True":
        bc = True
    else:
        bc = False    
    main(np.int(fbeg),np.int(fend),np.int(nstep),bc,np.int(scale),np.int(scale2))
