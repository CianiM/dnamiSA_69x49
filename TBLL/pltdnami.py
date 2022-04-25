import numpy as np
import matplotlib.pyplot as plt
import sys
# import mkl_fft as mklfft
# import tools_fft as tl
# from dnami import cst
from rhsinfo import wp
import math as m

# from owned_io import write_data

def main():
	
	nit = 14000

	# hlo = 2	

	# first data 
	# path = './restarts/'
	path = './output/'
	fname = 'output'
	# fname = 'restart'

	q  = build_q(path,fname,nit)

	# qc = build_q(path+'run1/',fname,nit)

    # second data	

	# [qi1_2]    = loadrstax(nit,rdir='./restarts/run1/',fname='restartshell_i1',noaxis=True)
	# [qimax_2]  = loadrstax(nit,rdir='./restarts/run1/',fname='restartshell_imax',noaxis=True)

	# [qj1_2]    = loadrstax(nit,rdir='./restarts/run1/',fname='restartshell_j1',noaxis=True)
	# [qjmax_2]  = loadrstax(nit,rdir='./restarts/run1/',fname='restartshell_jmax',noaxis=True)

	# [qcore_2]  = loadrstax(nit,rdir='./restarts/run1/',fname='restart',noaxis=True)

	# nxgb,nygb,nzgb,nvar = qcore_2.shape[0],qcore_2.shape[1],qcore_2.shape[2],qcore_2.shape[3]
	# hlo = qi1_2.shape[0]

	# print(nxgb,nygb,nzgb)

	# q_2 = np.empty(shape=(nxgb,nygb+2*hlo,nzgb,nvar),dtype = wp)

	# q_2[0:nxgb,0:hlo,0:nvar]               = qj1_2
	# q_2[0:nxgb,hlo:hlo+nygb,0:nvar]        = qcore_2
	# q_2[0:nxgb,hlo+nygb:2*hlo+nygb,0:nvar] = qjmax_2


	# [X,Y,Z,qtest] = loadrstax(100,rdir='/restarts/',noaxis=False)


	# k   = []
	# spc = []
	# label = ['Skew FD05 130**3','Skew FD11 130**3','Skew + DivDif FD11 130^3','Skew. + Flt (200 it) FD11' ,'Cons. + Flt FD11' ,'Skew FD11 260**3']
	# color = ['red','blue','gray','cyan','orange','k']
	# for q in [qfd11,qfd05,qfd11visc, qfd11flt, qfd11flt0, qref]:

	# 	u = q[:,:,:,1]
	# 	v = q[:,:,:,2]
	# 	w = q[:,:,:,3]
	
	
	# 	[out1,out2] = tl.power_spectra(u,v,w)	
	# 	k.append(out1[2:])
	# 	spc.append(out2[2:])

	# wref   = qref[:,:,:,3]
	# wstore = qstored[:,:,:,3]
	# print(wstore-wref)
	fig, ax = plt.subplots(figsize=(12,8))	


	print(np.transpose(q).shape)

	Lx = 1.0
	Ly = 1.0
	nxgb = q.shape[0]
	nygb = q.shape[0]

	ax.contour(np.transpose(q[:,:,0,2]),np.arange(-200.0,200.0,10.))
	# ax.contour(np.transpose(q[:,:,0,2]),np.arange(-2.0,2.0,0.1))

	# print(q[:,:,0,0])

	# ax.contour(np.transpose(q[:,:,0,0]),200)

	# print(q[:,-5,0,0])
	


	# ax.plot(np.transpose(qc[:,1,0,1]),'r+-')	

	# ax.contour(np.transpose(q[:,:,0,1]),40,colors='r',linestyles='--')	

	# ax.contour(np.transpose(q_2[:,:,0,1]-q[:,:,0,1]),np.arange(-5.0,5.0,0.4))

	# print(q[10,:,0,1])
	# for line in np.arange(0,5):
	# 	ax.plot(np.transpose(q[:,line,0,0]),'+-')
	# ax.plot(np.transpose(q_2[45,:,0,1]-q[45,:,0,1]),'r+-:')
	
	# ax.plot(np.transpose(q[:,10,0,0]),'+-')

	# ax.plot(xloc,-np.cos(xloc),'r+')
	# ax.plot(xloc,-np.sin(xloc),'r+')
	# print(qstored[:,52,89,0]-qref[:,52,89,0])
	# ax.plot(qstored[:,52,89,0]-qref[:,52,89,0])

	# T0 = q[:,:,0,0]
	# T1 = np.float64(1.)
	# T2 = np.float64(1.)
	
	# # - Analytical solution
	# x = np.linspace(cst(0.0),Lx,nxgb,dtype=np.float64)
	# y = np.linspace(cst(0.0),Ly,nygb,dtype=np.float64)

	# Tan = np.zeros((np.size(x), np.size(y)), dtype=np.float64)
	# Lx = np.amax(x)
	# Ly = np.amax(y)
	
	# for i,xv in enumerate(x):
	#     for j,yv in enumerate(y):
	#         Tan[i,j] = T1 + T2 * np.sinh(np.pi*yv/Lx)*np.sin(np.pi*xv/Lx)/np.sinh(np.pi*Ly/Lx) 
	
	# print('Max error:')
	# max_error = np.amax( np.abs(Tan[:, :] - T0 )   )

	# ax.plot(Tan[20, :],'r+-')
	# ax.plot(q[20,:,0,0],'bo')

	# print(max_error)
	

	plt.show()



def loadrstax(n,axis='./out/',rdir='./restarts/',fname=None,noaxis=False):
	if not noaxis:
		xyzsize    = np.fromfile(axis+'axes.bin', dtype=np.float64, count=3)
		nxgb       = np.int(xyzsize[0])
		nygb       = np.int(xyzsize[1])
		nzgb       = np.int(xyzsize[2])
		
		# print("[grid]    nxgb, nygb, nzgb, nvar: ",nxgb, nygb, nzgb)
	restartdir = rdir 

	if fname != None:
		fname = restartdir + fname + '_' + str(n).zfill(8)
	else:	
		fname = restartdir + 'restart_' + str(n).zfill(8)	
	headsize = 7
	with open(fname,"rb") as fh:
				head = np.fromfile(fh,dtype='float64',count=headsize)
				nxgb, nygb, nzgb, nvar = int(head[3]),int(head[4]),int(head[5]),int(head[6])
				# print("[restart] nxgb, nygb, nzgb, nvar: ",nxgb, nygb, nzgb, nvar)		
				dat  = np.fromfile(fh,dtype='float64',count=int(nxgb*nygb*nzgb*nvar))
				dat  = np.reshape(dat,(nxgb,nygb,nzgb,nvar))
	fh.closed
	
	if not noaxis:
		xyzsize    = np.fromfile(axis+'axes.bin', dtype=np.float64, count=nxgb*nygb*nzgb+6)
		x          = xyzsize[6:6+nxgb]
		y          = xyzsize[6+nxgb:6+nxgb+nygb]
		z          = xyzsize[6+nxgb+nygb:6+nxgb+nygb+nzgb]
		X,Y,Z = np.meshgrid(x,y,z, sparse=False,indexing='ij')
	# X,Y = np.meshgrid(x,y, sparse=False,indexing='ij')
	
	# print('[load restarts] restart file '+fname+' loaded.')
	Z = 1
	if not noaxis:
		return [X,Y,Z,dat]
	else:
		return[dat]	
	# return [X,Y,dat]


def loadrst_cmpreal(n):	
	# print("[grid]    nxgb, nygb, nzgb, nvar: ",nxgb, nygb, nzgb)
	restartdir = '/home/alferezn/Dropbox/Code/compreal_developments/CompReal/wrk/restarts/'

	fname = restartdir + 'restart_' + str(n).zfill(8)	

	with open(fname,"rb") as fh:

				bufsize = np.fromfile(fh,dtype='float64',count=1)
				buf     = np.fromfile(fh,dtype='float64',count=int(bufsize)-1)				
				nxgb, nygb, nzgb, nvar = int(buf[2]),int(buf[3]),int(buf[4]),int(buf[5])
				# print("[restart CompReal] nxgb, nygb, nzgb, nvar: ",nxgb, nygb, nzgb, nvar)		

				dat  = np.fromfile(fh,dtype='float64',count=int(nxgb*nygb*nzgb*nvar))
				dat  = np.reshape(dat,(nxgb,nygb,nzgb,nvar),order='F')
	fh.closed
		
	print('[load restarts from CompReal] restart file '+fname+' loaded.')

	return dat	

def build_q(path,fname,nit):

	[qi1]    = loadrstax(nit,rdir=path,fname=fname+'shell_i1',noaxis=True)
	[qimax]  = loadrstax(nit,rdir=path,fname=fname+'shell_imax',noaxis=True)

	[qj1]    = loadrstax(nit,rdir=path,fname=fname+'shell_j1',noaxis=True)
	[qjmax]  = loadrstax(nit,rdir=path,fname=fname+'shell_jmax',noaxis=True)

	[qcore]  = loadrstax(nit,rdir=path,fname=fname,noaxis=True)

	nxgb,nygb,nzgb,nvar = qcore.shape[0],qcore.shape[1],qcore.shape[2],qcore.shape[3]
	hlo = qi1.shape[0]

	q = np.empty(shape=(nxgb+2*hlo,nygb+2*hlo,nzgb,nvar),dtype = wp)

	q[0:hlo              ,0:nygb+2*hlo,0:nvar] = qi1
	q[hlo+nxgb:2*hlo+nxgb,0:nygb+2*hlo,0:nvar] = qimax

	q[0:nxgb+2*hlo,0:hlo              ,0:nvar] = qj1
	q[0:nxgb+2*hlo,hlo+nygb:2*hlo+nygb,0:nvar] = qjmax

	q[hlo:nxgb+hlo,hlo:hlo+nygb       ,0:nvar] = qcore

	return q

if __name__ == '__main__':
    main()
