# Plot the 2D results 

import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm

import pickle

from scipy.integrate import simps

from pltdnami import loadrstax

def main(fbeg,fend,nstep,bc,scale=0,scale2=0):
    
    from rhsinfo import hlo_rhs

    dt =  1.0-9
    Ma = 0.2
    gamma = 1.4
    P0 = 1.0/(Ma**2*gamma)

    wp = 'float64' # working precision
    nvar = 4
    iVDW = False
    ndim  = 2
    ivar = 1
    varcont = 4

    from pltdnami import loadrstax

    # [axis]  = loadrstax(0,rdir='./grid/',fname='output',noaxis=True)
    axis  = build_q('./grid/','output',0,bc)

    ksi = axis[:,:,0,0]
    eta = axis[:,:,0,1]
    
    nx = np.shape(ksi)[0]
    ny = np.shape(ksi)[1]

    # Plot Mesh:
    PlotMesh = False
    if PlotMesh:
        
        ibeg = hlo_rhs
        iend = nx-hlo_rhs

        jbeg = hlo_rhs
        jend = ny-hlo_rhs

        ibeg = 0
        iend = nx

        jbeg = 0
        jend = ny

        fig = plt.figure(figsize=(10,5))
        ax  = fig.add_subplot(111)       
        for i in range(iend):
           ax.plot( ksi[i,jbeg:jend], eta[i,jbeg:jend], 'k--', lw=0.1  )
        for j in range(jend):
           ax.plot( ksi[ibeg:iend,j], eta[ibeg:iend,j], 'k--', lw=0.1  )

        theta = np.linspace(0.,2.0*np.pi,100)
        ax.plot(np.cos(theta),np.sin(theta),'b-',lw=0.2)
           
        plt.axis('scaled')
        plt.show()



    vardic = ['rho','u','v','et','nut','tau_wall']


    # lvls = np.arange(-270.,270.,20)
    # lvls = np.arange(-50.,50.,2)    
    # lvls = np.arange(-0.04,0.08,0.005)*1e-5     
    # lvls = np.arange(-1.2,0.8,0.2)
    # lvls = np.arange(69,70,0.1)    
    
    lvls = np.arange(-4.,0.1,0.1) 

    norm  = 1.0#P0
    norm2 = P0
    Normalize = False

    cmap = plt.get_cmap('bwr')
    # cmap = plt.get_cmap('gray_r')
    # normColor = BoundaryNorm(np.arange(-0.1,0.002,0.001), ncolors=cmap.N, clip=True)
    # normColor = BoundaryNorm(np.arange(-1e-6,1e-6,1e-7), ncolors=cmap.N, clip=True)
    # normColor = BoundaryNorm(np.arange(-0.004,0.004,0.001), ncolors=cmap.N, clip=True)
    normColor = BoundaryNorm(np.arange(-0.9974297337335298,-0.9946505541624183,0.0001), ncolors=cmap.N, clip=True)

    # sys.exit()
    nn  = np.arange(fbeg, fend, nstep)
    
    nsnap = scale
    iAnim = True 
    iCut  = False
    compare = False


    time      = []
    enstrophy = []
    kinetic   = []

    ibeg = hlo_rhs
    iend = nx-hlo_rhs

    jbeg = hlo_rhs
    jend = ny-hlo_rhs

    ibeg = 0
    iend = nx

    jbeg = 0
    jend = 20

    if iAnim :

        fig = plt.figure(figsize=(10,5))
        ax  = fig.add_subplot(111)


        qdat = build_q('./restarts/'  ,'restart',scale,bc)
        #qdatO = build_q('./output/','output' ,scale,bc)
        qdatO = build_q('./output/','outputoutput',scale,bc)

        # qdatcmp = build_q('./output/omg/','output' ,scale2,bc)
        # qrst_mpi = build_q('./restarts_MPI/restarts/'  ,'restart',nsnap,bc)
        # qdat_mpi = build_q('./output/omg/','output' ,nsnap,bc)

        P00  = qdatO[int(nx/2),int(ny/2),0, ivar   ]
        Pinf = qdatO[0,0,0, ivar   ]      

        if Normalize:
            vcont  = qdatO[:,:,0,varcont]
            
            qp     = (qdatO[:,:,0,ivar   ] -Pinf)/abs(P00 - Pinf)

        else:
            vcont  = qdatO[:,:,0,0]
            # qp     = qdat[:,:,0,5] 
            # qp  = -  qdat[:,:,0,7]* qdat[:,:,0,0]  +  qdat[:,:,0,5] *  qdat[:,:,0,4]  
            qp    =     qdatO[:,:,0,ivar]   
            qp    =     qdatO[:,:,0,varcont]         
            # qpcmp =  qdatcmp[:,:,0,2]

            # qpmpi     = qrst_mpi[:,:,0,2]  

        print(np.min(qp),np.max(qp))

        iplt,jplt = 69-1,49-1
        TimeSerie = False
        if TimeSerie:

            datbase = []
            
            oneDLinex = []
            oneDLiney = []

            for n in range(scale,scale2,nstep):
                
                oneDLinex.append(n)

                datbase.append(build_q('./output/omg/','output' ,n,bc))
            
            for i,q in enumerate(datbase[1:]):
                print(np.amax(q[:,:,0,4]),np.unravel_index(np.argmax(q[:,:,0,4],axis=None),q[:,:,0,4].shape))
                # oneDLiney.append(np.amin(q[:,:,0,4]))
                oneDLiney.append(q[iplt,jplt,0,4]-datbase[i-1][iplt,jplt,0,4])


                
            ax.plot(oneDLinex[1:],oneDLiney,'ko-')    
            plt.show()
            sys.exit()    


        # qp = qp - qpcmp
        # qp[np.isnan(qp)] = 10000.
        
        # print(qp[0,:])
        # ax.plot(qp[40,:],eta[40,:],'ro-')
        # ax.plot(qp[40,:],eta[40,:],'ro-')
        # print(qp.shape)
        # ax.plot(qp[120,jbeg:jend],eta[100,jbeg:jend],'bo-')
        # ax.plot(qp[540,jbeg:jend],eta[540,jbeg:jend],'ro-')

        # ax.plot(qp[:,100],'ro-')

        # ax.plot(ksi[:,0],qp[:,0],'ro-')

        # ax.plot(qpcmp[40,:],'ko--')
        # ax.plot(qpcmp[80,:],'co--')

        # print(qpcmp[40,:])
        # print('\n\n')
        # print(qp[40,:])
        # print(np.amax(qp[:,:]),np.unravel_index(np.argmax(qp[:,:],axis=None),qp[:,:].shape))

        # plt.show();sys.exit()

        # print("Lines min max"  ,np.amin(vcont), np.amax(vcont))      
        # print("Coutour min max rho ",np.min(qdat[:,:,0,1]-qdatcmp[:,:,0,1])    , np.max(qdat[:,:,0,1]-qdatcmp[:,:,0,1])    )
        # print("Coutour min max u   ",np.min(qdat[:,:,0,0]-qdatcmp[:,:,0,0])    , np.max(qdat[:,:,0,0]-qdatcmp[:,:,0,0])    )
        # print("Coutour min max v   ",np.min(qdat[:,:,0,4]-qdatcmp[:,:,0,4])    , np.max(qdat[:,:,0,4]-qdatcmp[:,:,0,4])    )
        # print("Coutour min max P   ",np.min(qdat[:,:,0,2]-qdatcmp[:,:,0,2])    , np.max(qdat[:,:,0,2]-qdatcmp[:,:,0,2])    )

        # min = np.amin(qp)/10.

        # max = np.amax(qp)/1000.  
        # min = -max
        # max = -min

        # max = 1.e-3
        # min =-1.e-3

        max = np.amax(qp)
        min = np.amin(qp)
        print(min,max)
        # max = 1.00e-14
        # min = -1.0e-14

        normColor = BoundaryNorm(np.arange(min    , max,abs(max-min)/50.), ncolors=cmap.N, clip=True)
        
        # normColor = BoundaryNorm(np.arange(-0.5 , 0.5,1.0/10.), ncolors=cmap.N, clip=True)
        lvls = np.arange(0.0   , 0.1,abs(np.amax(vcont)-np.amin(vcont))/2.)

        im  = ax.pcolormesh(ksi[ibeg:iend,jbeg:jend],eta[ibeg:iend,jbeg:jend],   qp[ibeg:iend,jbeg:jend], cmap=cmap,norm=normColor);fig.colorbar(im,ax=ax)     
        # cnt = ax.contour(   ksi[ibeg:iend,jbeg:jend],eta[ibeg:iend,jbeg:jend],vcont[ibeg:iend,jbeg:jend],levels=lvls,colors='k')   
        # plt.axis('scaled')
        # im  = ax.pcolormesh(   qp[ibeg:iend,jbeg:jend], cmap=cmap,norm=normColor);fig.colorbar(im,ax=ax)     
        # cnt = ax.contour(   vcont[ibeg:iend,jbeg:jend],levels=lvls,colors='k')   


        # im  = ax.pcolormesh(ksi[ibeg:iend,jbeg:jend],eta[ibeg:iend,jbeg:jend],qp[ibeg:iend,jbeg:jend], cmap='Blues_r');fig.colorbar(im,ax=ax)

#           ax.scatter(ksi[iplt,jplt],eta[iplt,jplt]) 
        for i in range(iend):
           ax.plot( ksi[i,jbeg:jend], eta[i,jbeg:jend], 'k--', lw=0.1  )
        for j in range(jend):
           ax.plot( ksi[ibeg:iend,j], eta[ibeg:iend,j], 'k--', lw=0.1  )
        plt.show();sys.exit()

        

        if compare:    
            # im  = ax.pcolormesh(ksi[ibeg:iend,jbeg:jend],eta[hlo_rhs:ny-hlo_rhs,:],omg[ibeg:iend,jbeg:jend], cmap='Blues_r')
            axis2  = build_q('./grid480/grid/','output',0,bc)

            ksi2 = axis2[:,:,0,0]
            eta2 = axis2[:,:,0,1]

            nx2 = np.shape(ksi2)[0]
            ny2 = np.shape(ksi2)[1]
            qdat2 = build_q('./grid480/output/omg/','output',nsnap,bc)
            omg2  = qdat2[:,:,0,2]       
            cnt2 = ax.contour(ksi2[hlo_rhs:nx2-hlo_rhs,:],eta2[hlo_rhs:nx2-hlo_rhs,:],omg2[hlo_rhs:nx2-hlo_rhs,:],levels=lvls,colors='r') 
    
        # ax.axhline(y=0)
        # ax.axvline(x=0)
        
        # -- Add the grid
        # nx = np.shape(ksi)[0]
        # ny = np.shape(ksi)[1]
        # print(nx,ny)

        # for i in range(nx-hlo_rhs):
        #    ax.plot( ksi[i+hlo_rhs,hlo_rhs:ny-hlo_rhs], eta[i+hlo_rhs,hlo_rhs:ny-hlo_rhs], 'k--', lw=0.1  )
        # for j in range(ny-hlo_rhs):
        #    ax.plot( ksi[hlo_rhs:nx-hlo_rhs,j+hlo_rhs], eta[hlo_rhs:nx-hlo_rhs,j+hlo_rhs], 'k--', lw=0.1  )

        for i in range(iend):
           ax.plot( ksi[i,jbeg:jend], eta[i,jbeg:jend], 'k--', lw=0.1  )
        for j in range(jend):
           ax.plot( ksi[ibeg:iend,j], eta[ibeg:iend,j], 'k--', lw=0.1  )
        plt.show();sys.exit()

        #  ----- Plot colormesh
        for k,n in enumerate(nn):

            print('Making fig at n:', n)

            qrst  = build_q('./restarts/'  ,'restart',n,bc)
            qdat  = build_q('./output/omg/','output' ,n,bc)

            if Normalize:
                vcont  = qdat[:,:,0,varcont]
                qp     = (qdat[:,:,0,ivar   ] -Pinf)/abs(P00 - Pinf)
            else:
                vcont  = qdat[:,:,0,varcont]
                qp     = qdat[:,:,0,ivar   ]  

            #  Enstrophy
            datI = vcont*vcont
            
            # Perform integrals:
            Iy = []
            for j in np.arange(0,ny,1):
                Iy.append(simps(datI[hlo_rhs:nx-hlo_rhs,j],ksi[hlo_rhs:nx-hlo_rhs,j]))    
 
            Omega = simps(Iy[:],eta[hlo_rhs,:])

            datI = qrst[:,:,0,1]*qrst[:,:,0,1] + qrst[:,:,0,2]*qrst[:,:,0,2]

            Iy = []
            for j in np.arange(0,ny,1):
                Iy.append(simps(datI[hlo_rhs:nx-hlo_rhs,j],ksi[hlo_rhs:nx-hlo_rhs,j]))    
 
            E = simps(Iy[:],eta[hlo_rhs,:])
        
            print("Integral of enstrophy:",Omega/2.,"Integral of kinetic energy:",E/2.,"Time:",n*dt)

            time.append(     n*dt)
            enstrophy.append(Omega/2.)
            kinetic.append(  E/2.)

            # qp  = qrst[:,:,0,ivar]
            # print(np.amin(vcont), np.amax(vcont), np.shape(vcont))
            
            print("Lines min max"  ,np.amin(vcont), np.amax(vcont))      
            print("Coutour min max",np.min(qp)    , np.max(qp)    )            
                       
            
            im  = ax.pcolormesh(ksi[ibeg:iend,jbeg:jend],
                                eta[ibeg:iend,jbeg:jend],
                                 qp[ibeg:iend,jbeg:jend], cmap=cmap,norm=normColor)

            # im  = ax.pcolormesh(qp[ibeg:iend,jbeg:jend], cmap=cmap,norm=normColor)
             
            for coll in cnt.collections: 
                    plt.gca().collections.remove(coll)
            
            cnt = ax.contour(ksi[ibeg:iend,jbeg:jend],
                             eta[ibeg:iend,jbeg:jend],
                             vcont[ibeg:iend,jbeg:jend],levels=lvls,colors='k')

            # cnt = ax.contour(vcont[ibeg:iend,jbeg:jend],levels=lvls,colors='k')

            
            ax.set(title='Var: ' + vardic[ivar] + ' \n min/max: {:.4} {:.4}'.format(np.amin(qp), np.amax(qp)))
            plt.axis('scaled')
            figname = 'pics/omg_' + str(np.int(n/nstep)).zfill(4) + '.png'
            
            plt.savefig(figname,format='png',dpi=200)
    
    # -- Comparison
    if iCut:
        ivar = 2
        fig = plt.figure(figsize=(8,8))
        ax  = fig.add_subplot(111)
        datname = './restarts/restart_' + str(0).zfill(8) 
        nt,t,q = read_restart(datname)
        datname = './out/omg/output_' + str(0).zfill(8) 
        omg0 = read_coords(datname)
        qp0 = q[:,:,ivar]
    
        datname = './restarts/restart_' + str(595000).zfill(8) 
        nt,t,q = read_restart(datname)
        qp1 = q[:,:,ivar]
        datname = './out/omg/output_' + str(595000).zfill(8) 
        omg1 = read_coords(datname)
    
        nymid = int(np.size(qp0[0,:])/2.) 
        x = ksi[:,nymid] 
        uinf =  0.5916079783099616
        ax.plot(x, qp0[:,nymid]/uinf, ls='', marker='o', color='k', label=r'Initial condition')
        ax.plot(x, qp1[:,nymid]/uinf, ls='-', color='k', label=r'(9pt, 8th ord) - 11 travel times ')
        ax.set_xlabel(r'x/r_v')
        ax.set_ylabel(r'v/u_inf')
        ax.set_xlim([-5,5])
        ax.legend(loc='best')
    
        figname = 'pics/cut.png'
        plt.savefig(figname,format='png',dpi=200)

    pickle.dump([time,enstrophy,kinetic],open("average_ens_kin.bin","wb"))    
        


def build_qwbc(path,fname,nit):

    from rhsinfo import bc_info,wp
    from pltdnami import loadrstax

    bc = bc_info[1]

    sizebc_i = 0
    sizebc_j = 0

    if ('j1' in bc) or ('jmax' in bc):
        [qj1]    = loadrstax(nit,rdir=path,fname=fname+'shell_j1',noaxis=True)
        [qjmax]  = loadrstax(nit,rdir=path,fname=fname+'shell_jmax',noaxis=True)
        sizebc_j = 2*qj1.shape[1]
 
    if ('i1' in bc) or ('imax' in bc):   
        [qi1]    = loadrstax(nit,rdir=path,fname=fname+'shell_i1',noaxis=True)
        [qimax]  = loadrstax(nit,rdir=path,fname=fname+'shell_imax',noaxis=True)
        sizebc_i = 2*qi1.shape[0]
 
    [qcore]  = loadrstax(nit,rdir=path,fname=fname,noaxis=True)

    nxgb,nygb,nzgb,nvar = qcore.shape[0],qcore.shape[1],qcore.shape[2],qcore.shape[3]

    q = np.empty(shape=(nxgb+sizebc_i,nygb+sizebc_j,nzgb,nvar),dtype = wp)

    if ('i1' in bc) or ('imax' in bc):  
        hlo = int(sizebc_i/2)
        q[0:hlo              ,0:nygb+2*hlo,0:nvar] = qi1
        q[hlo+nxgb:2*hlo+nxgb,0:nygb+2*hlo,0:nvar] = qimax

    if ('j1' in bc) or ('jmax' in bc):   
        hlo = int(sizebc_j/2)
        # print(q[0:nxgb+2*hlo,0:hlo              ,0:nvar].shape,qj1.shape,nxgb)
        q[0:nxgb+2*hlo,0:hlo              ,0:nvar] = qj1
        q[0:nxgb+2*hlo,hlo+nygb:2*hlo+nygb,0:nvar] = qjmax
         
    q[int(sizebc_i/2):nxgb+int(sizebc_i/2),int(sizebc_j/2):int(sizebc_j/2)+nygb       ,0:nvar] = qcore   

    # q[hlo:nxgb+hlo,hlo:hlo+nygb       ,0:nvar] = qcore

    return q    

def build_q(path,fname,nit,bc=False):

    from rhsinfo import wp

    if bc:
        [qi1]    = loadrstax(nit,rdir=path,fname=fname+'shell_i1',noaxis=True)
        [qimax]  = loadrstax(nit,rdir=path,fname=fname+'shell_imax',noaxis=True)
    
        [qj1]    = loadrstax(nit,rdir=path,fname=fname+'shell_j1',noaxis=True)
        [qjmax]  = loadrstax(nit,rdir=path,fname=fname+'shell_jmax',noaxis=True)

    [qcore]  = loadrstax(nit,rdir=path,fname=fname,noaxis=True)

    nxgb,nygb,nzgb,nvar = qcore.shape[0],qcore.shape[1],qcore.shape[2],qcore.shape[3]
    
    if bc:

        hlo = qi1.shape[0]
    
        q = np.empty(shape=(nxgb+2*hlo,nygb+2*hlo,nzgb,nvar),dtype = wp)
    
        q[0:hlo              ,0:nygb+2*hlo,0:nvar] = qi1
        q[hlo+nxgb:2*hlo+nxgb,0:nygb+2*hlo,0:nvar] = qimax
    
        q[0:nxgb+2*hlo,0:hlo              ,0:nvar] = qj1
        q[0:nxgb+2*hlo,hlo+nygb:2*hlo+nygb,0:nvar] = qjmax
    
        q[hlo:nxgb+hlo,hlo:hlo+nygb       ,0:nvar] = qcore
    else:
        q = qcore    
    return q




if __name__ == '__main__':
    name,fbeg,fend,nstep,bc,scale,scale2 = sys.argv
    if bc=="True":
        bc = True
    else:
        bc = False    
    main(np.int(fbeg),np.int(fend),np.int(nstep),bc,np.int(scale),np.int(scale2))
