














!====================================================================================================
!
! General Boundary Conditions: compute Eqns with Boundary Scheme
!     
!====================================================================================================

subroutine boundarySchemestored_faces_imax(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

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



 real(wp) ::  d1_stemp_dx_0_nxp2p0p0jk,d1_stemp_dx_0_nxp2p0m1jk,d1_stemp_dx_0_nxp2p0m2jk &
            ,d1_stemp_dx_0_nxp2p0jk &
            ,d1_stemp_dy_0_nxp2p0jm1k,d1_stemp_dy_0_nxp2p0jp0k,d1_stemp_dy_0_nxp2p0jp1k &
            ,d1_stemp_dy_0_nxp2p0jk 

 real(wp) ::  d1_stemp_dx_0_nxp2m1m1jk,d1_stemp_dx_0_nxp2m1p0jk,d1_stemp_dx_0_nxp2m1p1jk &
            ,d1_stemp_dx_0_nxp2m1jk &
            ,d1_stemp_dy_0_nxp2m1jm1k,d1_stemp_dy_0_nxp2m1jp0k,d1_stemp_dy_0_nxp2m1jp1k &
            ,d1_stemp_dy_0_nxp2m1jk 

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
indbc(4)=ny

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
! Start building layers for BC : imax None None ************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 0 None None ************************************
!                                                           
!***********************************************************


     do j=idloop(3),idloop(4) 


!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None stemp *
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ((dabs(0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_nxp2p0p0jk = q(nx+2+0+0,j,indvars(3))

d1_stemp_dx_0_nxp2p0m1jk = q(nx+2+0-1,j,indvars(3))

d1_stemp_dx_0_nxp2p0m2jk = q(nx+2+0-2,j,indvars(3))

d1_stemp_dx_0_nxp2p0jk = 1.5_wp*d1_stemp_dx_0_nxp2p0p0jk-&
          2.0_wp*d1_stemp_dx_0_nxp2p0m1jk+&
          0.5_wp*d1_stemp_dx_0_nxp2p0m2jk

d1_stemp_dx_0_nxp2p0jk = d1_stemp_dx_0_nxp2p0jk*param_float(1)

d1_stemp_dy_0_nxp2p0jm1k = q(nx+2+0,j-1,indvars(2))

d1_stemp_dy_0_nxp2p0jp1k = q(nx+2+0,j+1,indvars(2))

d1_stemp_dy_0_nxp2p0jk = -&
          0.5_wp*d1_stemp_dy_0_nxp2p0jm1k+&
          0.5_wp*d1_stemp_dy_0_nxp2p0jp1k

d1_stemp_dy_0_nxp2p0jk = d1_stemp_dy_0_nxp2p0jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None stemp **************
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(4)) =  ((dabs(0.5_wp*(qst(nx+2+0,j,indvarsst(11))*(d1_stemp_dy_0_nxp2p0jk)-&
                    qst(nx+2+0,j,indvarsst(10))*(d1_stemp_dx_0_nxp2p0jk)))))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None SA ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (stemp+fv2*ReI*nut/(k**2.0_wp*eta**2.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None SA *****************
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(12)) =  (qst(nx+2+0,j,indvarsst(4))+&
                    (1.0_wp-&
                    (q(nx+2+0,j,indvars(5))/1.0_wp)/(1.0_wp+&
                    (q(nx+2+0,j,indvars(5))/1.0_wp)*((q(nx+2+0,j,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2+0,j,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))))*param_float(1 + 5)*q(nx+2+0,j,indvars(5))/(param_float(9 + 5)**2.0_wp*qst(nx+2+0,j,indvarsst(2))**2.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None SS ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (max(SA,0.3*stemp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None SS *****************
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(13)) =  (max(qst(nx+2+0,j,indvarsst(12)),0.3*qst(nx+2+0,j,indvarsst(4))))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None r *****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (min(ReI*(nut/(SS*k**2.0_wp*eta**2.0_wp)),10.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None r ******************
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(14)) =  (min(param_float(1 + 5)*(q(nx+2+0,j,indvars(5))/(qst(nx+2+0,j,indvarsst(13))*param_float(9 + 5)**2.0_wp*qst(nx+2+0,j,indvarsst(2))**2.0_wp)),10.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None g *****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (r+Cw2*(r**6.0_wp-r))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None g ******************
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(15)) =  (qst(nx+2+0,j,indvarsst(14))+&
                    param_float(11 + 5)*(qst(nx+2+0,j,indvarsst(14))**6.0_wp-&
                    qst(nx+2+0,j,indvarsst(14))))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None fw ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (g*((1.0_wp+Cw3**6.0_wp)/(g**6.0_wp+Cw3**6.0_wp))**(1.0_wp/6.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None fw *****************
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(16)) =  (qst(nx+2+0,j,indvarsst(15))*((1.0_wp+&
                    param_float(12 + 5)**6.0_wp)/(qst(nx+2+0,j,indvarsst(15))**6.0_wp+&
                    param_float(12 + 5)**6.0_wp))**(1.0_wp/6.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None tau_wall 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (visc_t*([u]_1y)*deltayI)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None tau_wall ***********
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(17)) =  ((1.0_wp)*(1.0_wp+&
                    ((q(nx+2+0,j,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2+0,j,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+2+0,j,indvars(5))/1.0_wp))*param_float(1 + 5)*(d1_stemp_dy_0_nxp2p0jk)*qst(nx+2+0,j,indvarsst(11)))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None visc_SA 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! nut*rho
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None visc_SA ************
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(18)) =  q(nx+2+0,j,indvars(5))*q(nx+2+0,j,indvars(1))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None visc_turb 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! nut*rho*fv1
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None visc_turb **********
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(19)) =  q(nx+2+0,j,indvars(5))*q(nx+2+0,j,indvars(1))*((q(nx+2+0,j,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2+0,j,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None Pressure 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (gamma_m1)*rho*(e)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None Pressure ***********
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(20)) =  (param_float(3 + 5))*q(nx+2+0,j,indvars(1))*((q(nx+2+0,j,indvars(4))-&
                    0.5_wp*(q(nx+2+0,j,indvars(2))*q(nx+2+0,j,indvars(2))+&
                    q(nx+2+0,j,indvars(3))*q(nx+2+0,j,indvars(3)))))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None chi_coeff 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! nut/visc*rho
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None chi_coeff **********
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(21)) =  q(nx+2+0,j,indvars(5))/1.0_wp*q(nx+2+0,j,indvars(1))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None Production 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! Cb1*(1.0_wp-ft2)*SS*rho*nut
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None Production *********
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(22)) =  param_float(6 + 5)*(1.0_wp-&
                    (param_float(16 + 5)*exp(-&
                    param_float(17 + 5)*(q(nx+2+0,j,indvars(5))/1.0_wp)**2.0_wp)))*qst(nx+2+0,j,indvarsst(13))*q(nx+2+0,j,indvars(1))*q(nx+2+0,j,indvars(5))

     enddo


!***********************************************************
!                                                           
! Start building layers for BC : imax None None ************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 1 None None ************************************
!                                                           
!***********************************************************


     do j=idloop(3),idloop(4) 


!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None stemp *
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ((dabs(0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_nxp2m1m1jk = q(nx+2-1-1,j,indvars(3))

d1_stemp_dx_0_nxp2m1p1jk = q(nx+2-1+1,j,indvars(3))

d1_stemp_dx_0_nxp2m1jk = -&
          0.5_wp*d1_stemp_dx_0_nxp2m1m1jk+&
          0.5_wp*d1_stemp_dx_0_nxp2m1p1jk

d1_stemp_dx_0_nxp2m1jk = d1_stemp_dx_0_nxp2m1jk*param_float(1)

d1_stemp_dy_0_nxp2m1jm1k = q(nx+2-1,j-1,indvars(2))

d1_stemp_dy_0_nxp2m1jp1k = q(nx+2-1,j+1,indvars(2))

d1_stemp_dy_0_nxp2m1jk = -&
          0.5_wp*d1_stemp_dy_0_nxp2m1jm1k+&
          0.5_wp*d1_stemp_dy_0_nxp2m1jp1k

d1_stemp_dy_0_nxp2m1jk = d1_stemp_dy_0_nxp2m1jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None stemp **************
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(4)) =  ((dabs(0.5_wp*(qst(nx+2-1,j,indvarsst(11))*(d1_stemp_dy_0_nxp2m1jk)-&
                    qst(nx+2-1,j,indvarsst(10))*(d1_stemp_dx_0_nxp2m1jk)))))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None SA ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (stemp+fv2*ReI*nut/(k**2.0_wp*eta**2.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None SA *****************
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(12)) =  (qst(nx+2-1,j,indvarsst(4))+&
                    (1.0_wp-&
                    (q(nx+2-1,j,indvars(5))/1.0_wp)/(1.0_wp+&
                    (q(nx+2-1,j,indvars(5))/1.0_wp)*((q(nx+2-1,j,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2-1,j,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))))*param_float(1 + 5)*q(nx+2-1,j,indvars(5))/(param_float(9 + 5)**2.0_wp*qst(nx+2-1,j,indvarsst(2))**2.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None SS ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (max(SA,0.3*stemp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None SS *****************
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(13)) =  (max(qst(nx+2-1,j,indvarsst(12)),0.3*qst(nx+2-1,j,indvarsst(4))))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None r *****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (min(ReI*(nut/(SS*k**2.0_wp*eta**2.0_wp)),10.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None r ******************
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(14)) =  (min(param_float(1 + 5)*(q(nx+2-1,j,indvars(5))/(qst(nx+2-1,j,indvarsst(13))*param_float(9 + 5)**2.0_wp*qst(nx+2-1,j,indvarsst(2))**2.0_wp)),10.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None g *****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (r+Cw2*(r**6.0_wp-r))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None g ******************
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(15)) =  (qst(nx+2-1,j,indvarsst(14))+&
                    param_float(11 + 5)*(qst(nx+2-1,j,indvarsst(14))**6.0_wp-&
                    qst(nx+2-1,j,indvarsst(14))))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None fw ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (g*((1.0_wp+Cw3**6.0_wp)/(g**6.0_wp+Cw3**6.0_wp))**(1.0_wp/6.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None fw *****************
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(16)) =  (qst(nx+2-1,j,indvarsst(15))*((1.0_wp+&
                    param_float(12 + 5)**6.0_wp)/(qst(nx+2-1,j,indvarsst(15))**6.0_wp+&
                    param_float(12 + 5)**6.0_wp))**(1.0_wp/6.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None tau_wall 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (visc_t*([u]_1y)*deltayI)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None tau_wall ***********
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(17)) =  ((1.0_wp)*(1.0_wp+&
                    ((q(nx+2-1,j,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2-1,j,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+2-1,j,indvars(5))/1.0_wp))*param_float(1 + 5)*(d1_stemp_dy_0_nxp2m1jk)*qst(nx+2-1,j,indvarsst(11)))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None visc_SA 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! nut*rho
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None visc_SA ************
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(18)) =  q(nx+2-1,j,indvars(5))*q(nx+2-1,j,indvars(1))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None visc_turb 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! nut*rho*fv1
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None visc_turb **********
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(19)) =  q(nx+2-1,j,indvars(5))*q(nx+2-1,j,indvars(1))*((q(nx+2-1,j,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2-1,j,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None Pressure 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (gamma_m1)*rho*(e)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None Pressure ***********
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(20)) =  (param_float(3 + 5))*q(nx+2-1,j,indvars(1))*((q(nx+2-1,j,indvars(4))-&
                    0.5_wp*(q(nx+2-1,j,indvars(2))*q(nx+2-1,j,indvars(2))+&
                    q(nx+2-1,j,indvars(3))*q(nx+2-1,j,indvars(3)))))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None chi_coeff 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! nut/visc*rho
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None chi_coeff **********
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(21)) =  q(nx+2-1,j,indvars(5))/1.0_wp*q(nx+2-1,j,indvars(1))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None Production 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! Cb1*(1.0_wp-ft2)*SS*rho*nut
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None Production *********
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(22)) =  param_float(6 + 5)*(1.0_wp-&
                    (param_float(16 + 5)*exp(-&
                    param_float(17 + 5)*(q(nx+2-1,j,indvars(5))/1.0_wp)**2.0_wp)))*qst(nx+2-1,j,indvarsst(13))*q(nx+2-1,j,indvars(1))*q(nx+2-1,j,indvars(5))

     enddo

    enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO

 end subroutine boundarySchemestored_faces_imax


