














!====================================================================================================
!
! General Boundary Conditions: compute Eqns with Boundary Scheme
!     
!====================================================================================================

subroutine boundarySchemestored_edges_i1_j1_0_1(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

  implicit none

integer,parameter :: wp = kind(0.0D0) ! working precision

!=======================================
!
! Addresses for integers parameters
!
!=======================================



!========================================
!
! Addresses for floating point parameters
!
!========================================


  real(wp), intent(in)    :: param_float(*)

  integer, intent(in) :: neq,neqst
  integer, intent(in) :: nx,ny,nz,hlo
  integer, intent(in) :: ind(1:neq+neqst) 
  integer, intent(in) :: idarray(6)
  integer, intent(in) :: sizeblck(3)
  integer, intent(in) :: nvar_f(3),nvar_e(3)
  
real(wp),intent(inout) :: q(idarray(1):idarray(2),idarray(3):idarray(4),neq),&
                        rhs(idarray(1):idarray(2),idarray(3):idarray(4),neq),&
                        qst(idarray(1):idarray(2),idarray(3):idarray(4),neqst)

real(wp),intent(inout) :: qface_i(idarray(3):idarray(4),nvar_f(1)),&
                       qface_j(idarray(1):idarray(2),nvar_f(2)),&
                       qface_k(1),&
                       qedge_ij(1),&
                       qedge_jk(1),&
                       qedge_ik(1)

! LOCAL VARIABLES
 
  integer :: i,j,k
  integer :: bi,bj,bk
  integer :: biend,bjend,bkend
  integer :: size_bi,size_bj,size_bk
  integer :: indbc(6),idloop(6)



 real(wp) ::  d1_stemp_dx_0_1m2p0p01m2p1k,d1_stemp_dx_0_1m2p0p11m2p1k,d1_stemp_dx_0_1m2p0p21m2p1k &
            ,d1_stemp_dx_0_1m2p01m2p1k &
            ,d1_stemp_dy_0_1m2p01m2p1m1k,d1_stemp_dy_0_1m2p01m2p1p0k,d1_stemp_dy_0_1m2p01m2p1p1k &
            ,d1_stemp_dy_0_1m2p01m2p1k 

  integer :: indvars(neq),indvarsst(neqst)


  indvars   = ind(1:neq)
  indvarsst = ind(1+neq:neq+neqst) 

!f2py intent(in)    :: qst,nx,ny,nz
!f2py intent(inout) :: q,rhs
      
size_bk = sizeblck(3)
size_bj = sizeblck(2)
size_bi = sizeblck(1)

indbc(1)=1
indbc(2)=1

indbc(3)=1
indbc(4)=1

indbc(5)=1
indbc(6)=1


!$OMP DO SCHEDULE(GUIDED,4) COLLAPSE(2)  
  do bk=indbc(5),indbc(6),size_bk  
    do bj=indbc(3),indbc(4),size_bj 
      do bi=indbc(1),indbc(2),size_bi 
    
   
idloop(6) = min( bk+size_bk, indbc(6)+1)-1
idloop(4) = min( bj+size_bj, indbc(4)+1)-1
idloop(2) = min( bi+size_bi, indbc(2)+1)-1

idloop(5) = bk
idloop(3) = bj
idloop(1) = bi 



!***********************************************************
!                                                           
! Start building layers for BC : i1 j1 None ****************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 0 1 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None stemp ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ((dabs(0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_1m2p0p01m2p1k = q(1-2+0+0,1-2+1,indvars(3))

d1_stemp_dx_0_1m2p0p11m2p1k = q(1-2+0+1,1-2+1,indvars(3))

d1_stemp_dx_0_1m2p0p21m2p1k = q(1-2+0+2,1-2+1,indvars(3))

d1_stemp_dx_0_1m2p01m2p1k = -&
          1.5_wp*d1_stemp_dx_0_1m2p0p01m2p1k+&
          2.0_wp*d1_stemp_dx_0_1m2p0p11m2p1k-&
          0.5_wp*d1_stemp_dx_0_1m2p0p21m2p1k

d1_stemp_dx_0_1m2p01m2p1k = d1_stemp_dx_0_1m2p01m2p1k*param_float(1)

d1_stemp_dy_0_1m2p01m2p1m1k = q(1-2+0,1-2+1-1,indvars(2))

d1_stemp_dy_0_1m2p01m2p1p1k = q(1-2+0,1-2+1+1,indvars(2))

d1_stemp_dy_0_1m2p01m2p1k = -&
          0.5_wp*d1_stemp_dy_0_1m2p01m2p1m1k+&
          0.5_wp*d1_stemp_dy_0_1m2p01m2p1p1k

d1_stemp_dy_0_1m2p01m2p1k = d1_stemp_dy_0_1m2p01m2p1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None stemp *****************
!                                                           
!***********************************************************


qst(1-2+0,1-2+1,indvarsst(4)) =  ((dabs(0.5_wp*(qst(1-2+0,1-2+1,indvarsst(11))*(d1_stemp_dy_0_1m2p01m2p1k)-&
                    qst(1-2+0,1-2+1,indvarsst(10))*(d1_stemp_dx_0_1m2p01m2p1k)))))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None SA *******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (stemp+fv2*ReI*nut/(k**2.0_wp*eta**2.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None SA ********************
!                                                           
!***********************************************************


qst(1-2+0,1-2+1,indvarsst(12)) =  (qst(1-2+0,1-2+1,indvarsst(4))+&
                    (1.0_wp-&
                    (q(1-2+0,1-2+1,indvars(5))/1.0_wp)/(1.0_wp+&
                    (q(1-2+0,1-2+1,indvars(5))/1.0_wp)*((q(1-2+0,1-2+1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+0,1-2+1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))))*param_float(1 + 5)*q(1-2+0,1-2+1,indvars(5))/(param_float(9 + 5)**2.0_wp*qst(1-2+0,1-2+1,indvarsst(2))**2.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None SS *******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (max(SA,0.3*stemp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None SS ********************
!                                                           
!***********************************************************


qst(1-2+0,1-2+1,indvarsst(13)) =  (max(qst(1-2+0,1-2+1,indvarsst(12)),0.3*qst(1-2+0,1-2+1,indvarsst(4))))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None r ********
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (min(ReI*(nut/(SS*k**2.0_wp*eta**2.0_wp)),10.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None r *********************
!                                                           
!***********************************************************


qst(1-2+0,1-2+1,indvarsst(14)) =  (min(param_float(1 + 5)*(q(1-2+0,1-2+1,indvars(5))/(qst(1-2+0,1-2+1,indvarsst(13))*param_float(9 + 5)**2.0_wp*qst(1-2+0,1-2+1,indvarsst(2))**2.0_wp)),10.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None g ********
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (r+Cw2*(r**6.0_wp-r))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None g *********************
!                                                           
!***********************************************************


qst(1-2+0,1-2+1,indvarsst(15)) =  (qst(1-2+0,1-2+1,indvarsst(14))+&
                    param_float(11 + 5)*(qst(1-2+0,1-2+1,indvarsst(14))**6.0_wp-&
                    qst(1-2+0,1-2+1,indvarsst(14))))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None fw *******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (g*((1.0_wp+Cw3**6.0_wp)/(g**6.0_wp+Cw3**6.0_wp))**(1.0_wp/6.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None fw ********************
!                                                           
!***********************************************************


qst(1-2+0,1-2+1,indvarsst(16)) =  (qst(1-2+0,1-2+1,indvarsst(15))*((1.0_wp+&
                    param_float(12 + 5)**6.0_wp)/(qst(1-2+0,1-2+1,indvarsst(15))**6.0_wp+&
                    param_float(12 + 5)**6.0_wp))**(1.0_wp/6.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None tau_wall *
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (visc_t*([u]_1y)*deltayI)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None tau_wall **************
!                                                           
!***********************************************************


qst(1-2+0,1-2+1,indvarsst(17)) =  ((1.0_wp)*(1.0_wp+&
                    ((q(1-2+0,1-2+1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+0,1-2+1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+0,1-2+1,indvars(5))/1.0_wp))*param_float(1 + 5)*(d1_stemp_dy_0_1m2p01m2p1k)*qst(1-2+0,1-2+1,indvarsst(11)))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None visc_SA **
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! nut*rho
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None visc_SA ***************
!                                                           
!***********************************************************


qst(1-2+0,1-2+1,indvarsst(18)) =  q(1-2+0,1-2+1,indvars(5))*q(1-2+0,1-2+1,indvars(1))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None visc_turb 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! nut*rho*fv1
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None visc_turb *************
!                                                           
!***********************************************************


qst(1-2+0,1-2+1,indvarsst(19)) =  q(1-2+0,1-2+1,indvars(5))*q(1-2+0,1-2+1,indvars(1))*((q(1-2+0,1-2+1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+0,1-2+1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None Pressure *
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (gamma_m1)*rho*(e)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None Pressure **************
!                                                           
!***********************************************************


qst(1-2+0,1-2+1,indvarsst(20)) =  (param_float(3 + 5))*q(1-2+0,1-2+1,indvars(1))*((q(1-2+0,1-2+1,indvars(4))-&
                    0.5_wp*(q(1-2+0,1-2+1,indvars(2))*q(1-2+0,1-2+1,indvars(2))+&
                    q(1-2+0,1-2+1,indvars(3))*q(1-2+0,1-2+1,indvars(3)))))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None chi_coeff 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! nut/visc*rho
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None chi_coeff *************
!                                                           
!***********************************************************


qst(1-2+0,1-2+1,indvarsst(21)) =  q(1-2+0,1-2+1,indvars(5))/1.0_wp*q(1-2+0,1-2+1,indvars(1))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None Production 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! Cb1*(1.0_wp-ft2)*SS*rho*nut
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None Production ************
!                                                           
!***********************************************************


qst(1-2+0,1-2+1,indvarsst(22)) =  param_float(6 + 5)*(1.0_wp-&
                    (param_float(16 + 5)*exp(-&
                    param_float(17 + 5)*(q(1-2+0,1-2+1,indvars(5))/1.0_wp)**2.0_wp)))*qst(1-2+0,1-2+1,indvarsst(13))*q(1-2+0,1-2+1,indvars(1))*q(1-2+0,1-2+1,indvars(5))


    enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO

 end subroutine boundarySchemestored_edges_i1_j1_0_1


