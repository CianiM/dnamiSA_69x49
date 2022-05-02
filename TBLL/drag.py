import numpy as np
import matplotlib.pyplot as plt
from pltdnami import loadrstax
from plot_result import build_q
import sys
import scipy.integrate as integrate
import math

def main(fbeg,fend,nstep,bc,scale=0,scale2=0):
    axis  = build_q('./grid/','output',0,bc)
    q = build_q('./output/','outputoutput',scale,bc)
    Re = 5000000
    ksi = axis[:,:,0,0]
    eta = axis[:,:,0,1]
    y=eta[0,:]
    x0=13
    field = { 'rho':1,
              'u'  :2,
              'v'  :3,
              'et' :4,
              'nut':5,
              'tau':6,
              'visc_SA':7,
              'visc_turb':8,
              'Pressure':9,
              'chi_coeff':10,
              'Production':11,
              'deltaxI':12,
              'deltayI':13}
    var_rho=field['rho']-1
    var_u=field['u']-1
    var_nut=field['nut']-1
    var_tau=field['tau']-1
    var_eddy=field['visc_turb']-1
    var_SA = field['visc_SA']-1
    x=ksi[x0:,0]
    
    dx=np.zeros(len(x)-1)
    integ=dx
    rho = q[56,:,0,var_rho]
    tau_wall=q[x0:,0,0,var_tau]
    u=q[56,:,0,var_u]
    tw = 1/Re*q[56,0,0,var_tau]
    eddy_visc=q[56,0:50,0,var_eddy]
    SA_visc = q[56,0:50,0,var_SA]
    uw=np.sqrt(tw/rho[1:40])
    #print(tw)
    #print(uw)
    
    up = np.sqrt(1/Re)*u[1:40]/uw
    
    
    #dstar = 1.0/(rho[1:40]*uw)
    yp = np.sqrt(Re)*y[1:40]/(1.0/np.sqrt(tw/rho[1:40]))
    #print(yp)
    
    for i in np.arange(0,len(yp)):
        yp[i]=math.log10(yp[i])
    
    for i in np.arange(0,len(x)-1):
        dx[i]=x[i+1]-x[i]
        integ[i]=dx[i]*(tau_wall[i+1]+tau_wall[i])/2
    
    cf = 2*q[x0:,0,0,var_tau]
    cfx = 2*q[56,0,0,var_tau]
    cf_tot=2*2*np.sum(integ)
#    with open('x_coord.dat','r') as f:
#        x_axis = []
#        coord=f.readlines()
#        x_axis = [ float(s) for s in coord]
#        x_axis = np.array(x_axis)
    with open('flatplate_u.dat','r') as f:
        dat=f.readlines()[2:387]
        data = [word for line in dat for word in line.split()] 
        uc = [ float(s) for s in data]
        uc = np.array(uc)
        uc = uc.reshape((385,2))
 
    with open('flatplate_u+y+.dat','r') as f:
        dat=f.readlines()[2:385]
        data = [word for line in dat for word in line.split()] 
        ypup = [ float(s) for s in data]
        ypup = np.array(ypup)
        ypup = ypup.reshape((383,2))

    with open('cf_plate.dat','r') as f:
        dat=f.readlines()[2:449]
        data = [word for line in dat for word in line.split()] 
        xcf = [ float(s) for s in data]
        xcf = np.array(xcf)
        xcf = xcf.reshape((447,2))

    u_nasa = uc[:,0]
    y_nasa = uc[:,1]

    yp_nasa= ypup[:,0]
    up_nasa= ypup[:,1]
    #up_nasa = [math.log10(s) for s in up_nasa]
    #for i in np.arange(0,len(up_nasa)):
    #    up_nasa[i]=math.log10(up_nasa[i])
    x_nasa = xcf[:,0]
    cf_nasa = xcf[:,1]
    
    
    fig1 = plt.figure(figsize=(5,5))
    ax = fig1.add_subplot(111)
    velocity = False
    veloclog = False
    friction = True
    visc     = False
    SA       = False
    if velocity:
        im = ax.plot(u_nasa,y_nasa,'r')
        ax.plot(u,y,'k--')
    if veloclog:
        #print(up)
        #print(yp)
        im = ax.plot(yp_nasa,up_nasa,'r')
        ax.plot(yp,up,'k--')
        #ax.set_xscale('log')
    if friction:
        im = ax.plot(x,tau_wall,'k--')
        ax.plot(x_nasa,cf_nasa,'r')
        print(cf_tot)
    if visc :
        im = ax.plot(eddy_visc[0:35],y[0:35],'k--')
    if SA :
        im = ax.plot(SA_visc,y,'k--')
    #fig2 = plt.figure(figsize=(5,5))
    plt.show();sys.exit()


if __name__ == '__main__':
    name,fbeg,fend,nstep,bc,scale,scale2 = sys.argv
    if bc=="True":
        bc = True
    else:
        bc = False    
    main(np.int(fbeg),np.int(fend),np.int(nstep),bc,np.int(scale),np.int(scale2))
