














!====================================================================================================
!
! General Boundary Conditions: enforce physical BC and compute Eqns with Boundary Scheme
!                              Both steady and unsteady BC can be generated (i.e. on vector q or rhs)
!     
!====================================================================================================

subroutine phyLow_surfbc_faces_j1_6(param_float,hlo,neq,neqst,ind,idloop,idarray,nx,ny,nz,sizeblck,q,qst,q1,q2,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

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
  integer, intent(inout) :: idloop(6) 
  integer, intent(in) :: sizeblck(3)
  integer, intent(in) :: nvar_f(3),nvar_e(3)
    
real(wp),intent(inout) ::  q(1-hlo:nx+hlo,1-hlo:ny+hlo,neq),&
                          q1(1-hlo:nx+hlo,1-hlo:ny+hlo,neq),&
                          q2(1-hlo:nx+hlo,1-hlo:ny+hlo,neq),&
                         rhs(1-hlo:nx+hlo,1-hlo:ny+hlo,neq),&
                        qst(1-hlo:nx+hlo,1-hlo:ny+hlo,neqst)

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
  integer :: indbc(6) 



 real(wp) ::  d1_rhs_rho_dx_0_im11m4p0k,d1_rhs_rho_dx_0_ip01m4p0k,d1_rhs_rho_dx_0_ip11m4p0k &
            ,d1_rhs_rho_dx_0_i1m4p0k &
            ,d1_rhs_rho_dy_0_i1m4p0p0k,d1_rhs_rho_dy_0_i1m4p0p1k,d1_rhs_rho_dy_0_i1m4p0p2k &
            ,d1_rhs_rho_dy_0_i1m4p0k &
            ,d2_rhs_u_dxdx_0_0_im11m4p0k_im1m11m4p0k,d2_rhs_u_dxdx_0_0_im11m4p0k_im1p01m4p0k,d2_rhs_u_dxdx_0_0_im11m4p0k_im1p11m4p0k &
            ,d2_rhs_u_dxdx_0_0_im11m4p0k &
            ,d2_rhs_u_dxdx_0_0_ip11m4p0k_ip1m11m4p0k,d2_rhs_u_dxdx_0_0_ip11m4p0k_ip1p01m4p0k,d2_rhs_u_dxdx_0_0_ip11m4p0k_ip1p11m4p0k &
            ,d2_rhs_u_dxdx_0_0_ip11m4p0k &
            ,d1_rhs_u_dx_1_im11m4p0k,d1_rhs_u_dx_1_ip01m4p0k,d1_rhs_u_dx_1_ip11m4p0k &
            ,d1_rhs_u_dx_1_i1m4p0k &
            ,d1_rhs_u_dx_0_im11m4p0k,d1_rhs_u_dx_0_ip01m4p0k,d1_rhs_u_dx_0_ip11m4p0k &
            ,d1_rhs_u_dx_0_i1m4p0k &
            ,d2_rhs_v_dxdx_0_0_im11m4p0k_im1m11m4p0k,d2_rhs_v_dxdx_0_0_im11m4p0k_im1p01m4p0k,d2_rhs_v_dxdx_0_0_im11m4p0k_im1p11m4p0k &
            ,d2_rhs_v_dxdx_0_0_im11m4p0k &
            ,d2_rhs_v_dxdx_0_0_ip11m4p0k_ip1m11m4p0k,d2_rhs_v_dxdx_0_0_ip11m4p0k_ip1p01m4p0k,d2_rhs_v_dxdx_0_0_ip11m4p0k_ip1p11m4p0k &
            ,d2_rhs_v_dxdx_0_0_ip11m4p0k &
            ,d1_rhs_v_dx_1_im11m4p0k,d1_rhs_v_dx_1_ip01m4p0k,d1_rhs_v_dx_1_ip11m4p0k &
            ,d1_rhs_v_dx_1_i1m4p0k &
            ,d1_rhs_v_dx_0_im11m4p0k,d1_rhs_v_dx_0_ip01m4p0k,d1_rhs_v_dx_0_ip11m4p0k &
            ,d1_rhs_v_dx_0_i1m4p0k &
            ,d2_rhs_et_dxdx_0_2_im11m4p0k_im1m11m4p0k,d2_rhs_et_dxdx_0_2_im11m4p0k_im1p01m4p0k,d2_rhs_et_dxdx_0_2_im11m4p0k_im1p11m4p0k &
            ,d2_rhs_et_dxdx_0_2_im11m4p0k &
            ,d2_rhs_et_dxdx_0_2_ip11m4p0k_ip1m11m4p0k,d2_rhs_et_dxdx_0_2_ip11m4p0k_ip1p01m4p0k,d2_rhs_et_dxdx_0_2_ip11m4p0k_ip1p11m4p0k &
            ,d2_rhs_et_dxdx_0_2_ip11m4p0k &
            ,d1_rhs_et_dx_1_im11m4p0k,d1_rhs_et_dx_1_ip01m4p0k,d1_rhs_et_dx_1_ip11m4p0k &
            ,d1_rhs_et_dx_1_i1m4p0k &
            ,d1_rhs_et_dx_0_im11m4p0k,d1_rhs_et_dx_0_ip01m4p0k,d1_rhs_et_dx_0_ip11m4p0k &
            ,d1_rhs_et_dx_0_i1m4p0k &
            ,d1_rhs_et_dx_2_im11m4p0k,d1_rhs_et_dx_2_ip01m4p0k,d1_rhs_et_dx_2_ip11m4p0k &
            ,d1_rhs_et_dx_2_i1m4p0k &
            ,d1_rhs_et_dy_0_i1m4p0p0k,d1_rhs_et_dy_0_i1m4p0p1k,d1_rhs_et_dy_0_i1m4p0p2k &
            ,d1_rhs_et_dy_0_i1m4p0k &
            ,d2_rhs_et_dydy_4_0_i1m4p0p0k_i1m4p0p0p0k,d2_rhs_et_dydy_4_0_i1m4p0p0k_i1m4p0p0p1k,d2_rhs_et_dydy_4_0_i1m4p0p0k_i1m4p0p0p2k &
            ,d2_rhs_et_dydy_4_0_i1m4p0p0k &
            ,d2_rhs_et_dydy_4_0_i1m4p0p1k_i1m4p0p1m1k,d2_rhs_et_dydy_4_0_i1m4p0p1k_i1m4p0p1p0k,d2_rhs_et_dydy_4_0_i1m4p0p1k_i1m4p0p1p1k &
            ,d2_rhs_et_dydy_4_0_i1m4p0p1k &
            ,d2_rhs_et_dydy_4_0_i1m4p0p2k_i1m4p0p2m1k,d2_rhs_et_dydy_4_0_i1m4p0p2k_i1m4p0p2p0k,d2_rhs_et_dydy_4_0_i1m4p0p2k_i1m4p0p2p1k &
            ,d2_rhs_et_dydy_4_0_i1m4p0p2k &
            ,d1_rhs_et_dy_4_i1m4p0p0k,d1_rhs_et_dy_4_i1m4p0p1k,d1_rhs_et_dy_4_i1m4p0p2k &
            ,d1_rhs_et_dy_4_i1m4p0k &
            ,d1_rhs_nut_dx_0_im11m4p0k,d1_rhs_nut_dx_0_ip01m4p0k,d1_rhs_nut_dx_0_ip11m4p0k &
            ,d1_rhs_nut_dx_0_i1m4p0k &
            ,d1_rhs_nut_dx_1_im11m4p0k,d1_rhs_nut_dx_1_ip01m4p0k,d1_rhs_nut_dx_1_ip11m4p0k &
            ,d1_rhs_nut_dx_1_i1m4p0k &
            ,d2_rhs_nut_dxdx_2_0_im11m4p0k_im1m11m4p0k,d2_rhs_nut_dxdx_2_0_im11m4p0k_im1p01m4p0k,d2_rhs_nut_dxdx_2_0_im11m4p0k_im1p11m4p0k &
            ,d2_rhs_nut_dxdx_2_0_im11m4p0k &
            ,d2_rhs_nut_dxdx_2_0_ip11m4p0k_ip1m11m4p0k,d2_rhs_nut_dxdx_2_0_ip11m4p0k_ip1p01m4p0k,d2_rhs_nut_dxdx_2_0_ip11m4p0k_ip1p11m4p0k &
            ,d2_rhs_nut_dxdx_2_0_ip11m4p0k &
            ,d1_rhs_nut_dx_2_im11m4p0k,d1_rhs_nut_dx_2_ip01m4p0k,d1_rhs_nut_dx_2_ip11m4p0k &
            ,d1_rhs_nut_dx_2_i1m4p0k 

  integer :: indvars(neq),indvarsst(neqst)
  
  indvars   = ind(1:neq)
  indvarsst = ind(1+neq:neq+neqst) 

!f2py intent(in)    :: q1,q2,qst,nx,ny,nz,rhs,nrk
!f2py intent(inout) :: q

      
size_bk = sizeblck(3)
size_bj = sizeblck(2)
size_bi = sizeblck(1)

indbc(1)=1
indbc(2)=nx

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
! Start building layers for BC : None j1 None **************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: None 0 None ************************************
!                                                           
!***********************************************************


 
      do i=idloop(1),idloop(2) 


!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None d(rho)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (deltaxI*([rho*u]_1x))*symm+(rho*([v]_1y)*deltayI)*wall
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_rhs_rho_dx_0_im11m4p0k = q(i-1,1-4+0,indvars(1))*q(i-1,1-4+0,indvars(2))

d1_rhs_rho_dx_0_ip11m4p0k = q(i+1,1-4+0,indvars(1))*q(i+1,1-4+0,indvars(2))

d1_rhs_rho_dx_0_i1m4p0k = -&
          0.5_wp*d1_rhs_rho_dx_0_im11m4p0k+&
          0.5_wp*d1_rhs_rho_dx_0_ip11m4p0k

d1_rhs_rho_dx_0_i1m4p0k = d1_rhs_rho_dx_0_i1m4p0k*param_float(1)

d1_rhs_rho_dy_0_i1m4p0p0k = q(i,1-4+0+0,indvars(3))

d1_rhs_rho_dy_0_i1m4p0p1k = q(i,1-4+0+1,indvars(3))

d1_rhs_rho_dy_0_i1m4p0p2k = q(i,1-4+0+2,indvars(3))

d1_rhs_rho_dy_0_i1m4p0k = -&
          1.5_wp*d1_rhs_rho_dy_0_i1m4p0p0k+&
          2.0_wp*d1_rhs_rho_dy_0_i1m4p0p1k-&
          0.5_wp*d1_rhs_rho_dy_0_i1m4p0p2k

d1_rhs_rho_dy_0_i1m4p0k = d1_rhs_rho_dy_0_i1m4p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None d(rho)/dt **********
!                                                           
!***********************************************************


rhs(i,1-4+0,indvars(1)) =   -  ( (qst(i,1-4+0,indvarsst(10))*(d1_rhs_rho_dx_0_i1m4p0k))*qst(i,1-4+0,indvarsst(5))+&
                    (q(i,1-4+0,indvars(1))*(d1_rhs_rho_dy_0_i1m4p0k)*qst(i,1-4+0,indvarsst(11)))*dabs(1.0_wp-&
                    qst(i,1-4+0,indvarsst(5))) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None d(rho u)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (deltaxI*([-4.0_wp/3.0_wp*visc_t*({u}_1x)*deltaxI]_1x)+deltaxI*([rho*u*u+p]_1x))*symm
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_rhs_u_dxdx_0_0_im11m4p0k_im1m11m4p0k = q(i-1-1,1-4+0,indvars(2))

d2_rhs_u_dxdx_0_0_im11m4p0k_im1p11m4p0k = q(i-1+1,1-4+0,indvars(2))

d2_rhs_u_dxdx_0_0_im11m4p0k = -&
          0.5_wp*d2_rhs_u_dxdx_0_0_im11m4p0k_im1m11m4p0k+&
          0.5_wp*d2_rhs_u_dxdx_0_0_im11m4p0k_im1p11m4p0k

d2_rhs_u_dxdx_0_0_im11m4p0k = d2_rhs_u_dxdx_0_0_im11m4p0k*param_float(1)

d2_rhs_u_dxdx_0_0_ip11m4p0k_ip1m11m4p0k = q(i+1-1,1-4+0,indvars(2))

d2_rhs_u_dxdx_0_0_ip11m4p0k_ip1p11m4p0k = q(i+1+1,1-4+0,indvars(2))

d2_rhs_u_dxdx_0_0_ip11m4p0k = -&
          0.5_wp*d2_rhs_u_dxdx_0_0_ip11m4p0k_ip1m11m4p0k+&
          0.5_wp*d2_rhs_u_dxdx_0_0_ip11m4p0k_ip1p11m4p0k

d2_rhs_u_dxdx_0_0_ip11m4p0k = d2_rhs_u_dxdx_0_0_ip11m4p0k*param_float(1)

d1_rhs_u_dx_1_im11m4p0k = q(i-1,1-4+0,indvars(1))*q(i-1,1-4+0,indvars(2))*q(i-1,1-4+0,indvars(2))+(param_float(3 + 5))*q(i-1,1-4+0,indvars(1))*((q(i-1,1-4+0,indvars(4))-&
                    0.5_wp*(q(i-1,1-4+0,indvars(2))*q(i-1,1-4+0,indvars(2))+&
                    q(i-1,1-4+0,indvars(3))*q(i-1,1-4+0,indvars(3)))))

d1_rhs_u_dx_1_ip11m4p0k = q(i+1,1-4+0,indvars(1))*q(i+1,1-4+0,indvars(2))*q(i+1,1-4+0,indvars(2))+(param_float(3 + 5))*q(i+1,1-4+0,indvars(1))*((q(i+1,1-4+0,indvars(4))-&
                    0.5_wp*(q(i+1,1-4+0,indvars(2))*q(i+1,1-4+0,indvars(2))+&
                    q(i+1,1-4+0,indvars(3))*q(i+1,1-4+0,indvars(3)))))

d1_rhs_u_dx_1_i1m4p0k = -&
          0.5_wp*d1_rhs_u_dx_1_im11m4p0k+&
          0.5_wp*d1_rhs_u_dx_1_ip11m4p0k

d1_rhs_u_dx_1_i1m4p0k = d1_rhs_u_dx_1_i1m4p0k*param_float(1)

d1_rhs_u_dx_0_im11m4p0k = -4.0_wp/3.0_wp*(1.0_wp)*(1.0_wp+&
                    ((q(i-1,1-4+0,indvars(5))/1.0_wp*q(i-1,1-4+0,indvars(1)))**3.0_wp/((q(i-1,1-4+0,indvars(5))/1.0_wp*q(i-1,1-4+0,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i-1,1-4+0,indvars(5))/1.0_wp*q(i-1,1-4+0,indvars(1))))*param_float(1 + 5)*(d2_rhs_u_dxdx_0_0_im11m4p0k)*qst(i-1,1-4+0,indvarsst(10))

d1_rhs_u_dx_0_ip11m4p0k = -4.0_wp/3.0_wp*(1.0_wp)*(1.0_wp+&
                    ((q(i+1,1-4+0,indvars(5))/1.0_wp*q(i+1,1-4+0,indvars(1)))**3.0_wp/((q(i+1,1-4+0,indvars(5))/1.0_wp*q(i+1,1-4+0,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i+1,1-4+0,indvars(5))/1.0_wp*q(i+1,1-4+0,indvars(1))))*param_float(1 + 5)*(d2_rhs_u_dxdx_0_0_ip11m4p0k)*qst(i+1,1-4+0,indvarsst(10))

d1_rhs_u_dx_0_i1m4p0k = -&
          0.5_wp*d1_rhs_u_dx_0_im11m4p0k+&
          0.5_wp*d1_rhs_u_dx_0_ip11m4p0k

d1_rhs_u_dx_0_i1m4p0k = d1_rhs_u_dx_0_i1m4p0k*param_float(1)



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None d(rho u)/dt ********
!                                                           
!***********************************************************


rhs(i,1-4+0,indvars(2)) =   -  ( (qst(i,1-4+0,indvarsst(10))*(d1_rhs_u_dx_0_i1m4p0k)+&
                    qst(i,1-4+0,indvarsst(10))*(d1_rhs_u_dx_1_i1m4p0k))*qst(i,1-4+0,indvarsst(5)) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None d(rho v)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (deltaxI*([-visc_t*({v}_1x)*deltaxI]_1x)+deltaxI*([rho*u*v]_1x))*symm
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_rhs_v_dxdx_0_0_im11m4p0k_im1m11m4p0k = q(i-1-1,1-4+0,indvars(3))

d2_rhs_v_dxdx_0_0_im11m4p0k_im1p11m4p0k = q(i-1+1,1-4+0,indvars(3))

d2_rhs_v_dxdx_0_0_im11m4p0k = -&
          0.5_wp*d2_rhs_v_dxdx_0_0_im11m4p0k_im1m11m4p0k+&
          0.5_wp*d2_rhs_v_dxdx_0_0_im11m4p0k_im1p11m4p0k

d2_rhs_v_dxdx_0_0_im11m4p0k = d2_rhs_v_dxdx_0_0_im11m4p0k*param_float(1)

d2_rhs_v_dxdx_0_0_ip11m4p0k_ip1m11m4p0k = q(i+1-1,1-4+0,indvars(3))

d2_rhs_v_dxdx_0_0_ip11m4p0k_ip1p11m4p0k = q(i+1+1,1-4+0,indvars(3))

d2_rhs_v_dxdx_0_0_ip11m4p0k = -&
          0.5_wp*d2_rhs_v_dxdx_0_0_ip11m4p0k_ip1m11m4p0k+&
          0.5_wp*d2_rhs_v_dxdx_0_0_ip11m4p0k_ip1p11m4p0k

d2_rhs_v_dxdx_0_0_ip11m4p0k = d2_rhs_v_dxdx_0_0_ip11m4p0k*param_float(1)

d1_rhs_v_dx_1_im11m4p0k = q(i-1,1-4+0,indvars(1))*q(i-1,1-4+0,indvars(2))*q(i-1,1-4+0,indvars(3))

d1_rhs_v_dx_1_ip11m4p0k = q(i+1,1-4+0,indvars(1))*q(i+1,1-4+0,indvars(2))*q(i+1,1-4+0,indvars(3))

d1_rhs_v_dx_1_i1m4p0k = -&
          0.5_wp*d1_rhs_v_dx_1_im11m4p0k+&
          0.5_wp*d1_rhs_v_dx_1_ip11m4p0k

d1_rhs_v_dx_1_i1m4p0k = d1_rhs_v_dx_1_i1m4p0k*param_float(1)

d1_rhs_v_dx_0_im11m4p0k = -(1.0_wp)*(1.0_wp+&
                    ((q(i-1,1-4+0,indvars(5))/1.0_wp*q(i-1,1-4+0,indvars(1)))**3.0_wp/((q(i-1,1-4+0,indvars(5))/1.0_wp*q(i-1,1-4+0,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i-1,1-4+0,indvars(5))/1.0_wp*q(i-1,1-4+0,indvars(1))))*param_float(1 + 5)*(d2_rhs_v_dxdx_0_0_im11m4p0k)*qst(i-1,1-4+0,indvarsst(10))

d1_rhs_v_dx_0_ip11m4p0k = -(1.0_wp)*(1.0_wp+&
                    ((q(i+1,1-4+0,indvars(5))/1.0_wp*q(i+1,1-4+0,indvars(1)))**3.0_wp/((q(i+1,1-4+0,indvars(5))/1.0_wp*q(i+1,1-4+0,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i+1,1-4+0,indvars(5))/1.0_wp*q(i+1,1-4+0,indvars(1))))*param_float(1 + 5)*(d2_rhs_v_dxdx_0_0_ip11m4p0k)*qst(i+1,1-4+0,indvarsst(10))

d1_rhs_v_dx_0_i1m4p0k = -&
          0.5_wp*d1_rhs_v_dx_0_im11m4p0k+&
          0.5_wp*d1_rhs_v_dx_0_ip11m4p0k

d1_rhs_v_dx_0_i1m4p0k = d1_rhs_v_dx_0_i1m4p0k*param_float(1)



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None d(rho v)/dt ********
!                                                           
!***********************************************************


rhs(i,1-4+0,indvars(3)) =   -  ( (qst(i,1-4+0,indvarsst(10))*(d1_rhs_v_dx_0_i1m4p0k)+&
                    qst(i,1-4+0,indvarsst(10))*(d1_rhs_v_dx_1_i1m4p0k))*qst(i,1-4+0,indvarsst(5)) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None d(rho et)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (deltaxI*([-4.0_wp/3.0_wp*visc_t*({u}_1x)*deltaxI*u-visc_t*({v}_1x)*deltaxI*v-kappa*({T}_1x)*deltaxI]_1x)+deltaxI*([(rho*et+p)*u]_1x))*symm+(-visc_t*([u]_1y)*([u]_1y)*deltayI**2-4.0_wp/3.0_wp*visc_t*([v]_1y)*([v]_1y)*deltayI**2-kappa*([({T}_1y)*deltayI]_1y)-kappa*([({T}_1x)*deltaxI]_1x)+(rho*et+p)*([v]_1y)*deltayI)*wall
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_rhs_et_dxdx_0_2_im11m4p0k_im1m11m4p0k = ((q(i-1-1,1-4+0,indvars(4))-&
                    0.5_wp*(q(i-1-1,1-4+0,indvars(2))*q(i-1-1,1-4+0,indvars(2))+&
                    q(i-1-1,1-4+0,indvars(3))*q(i-1-1,1-4+0,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dxdx_0_2_im11m4p0k_im1p11m4p0k = ((q(i-1+1,1-4+0,indvars(4))-&
                    0.5_wp*(q(i-1+1,1-4+0,indvars(2))*q(i-1+1,1-4+0,indvars(2))+&
                    q(i-1+1,1-4+0,indvars(3))*q(i-1+1,1-4+0,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dxdx_0_2_im11m4p0k = -&
          0.5_wp*d2_rhs_et_dxdx_0_2_im11m4p0k_im1m11m4p0k+&
          0.5_wp*d2_rhs_et_dxdx_0_2_im11m4p0k_im1p11m4p0k

d2_rhs_et_dxdx_0_2_im11m4p0k = d2_rhs_et_dxdx_0_2_im11m4p0k*param_float(1)

d2_rhs_et_dxdx_0_2_ip11m4p0k_ip1m11m4p0k = ((q(i+1-1,1-4+0,indvars(4))-&
                    0.5_wp*(q(i+1-1,1-4+0,indvars(2))*q(i+1-1,1-4+0,indvars(2))+&
                    q(i+1-1,1-4+0,indvars(3))*q(i+1-1,1-4+0,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dxdx_0_2_ip11m4p0k_ip1p11m4p0k = ((q(i+1+1,1-4+0,indvars(4))-&
                    0.5_wp*(q(i+1+1,1-4+0,indvars(2))*q(i+1+1,1-4+0,indvars(2))+&
                    q(i+1+1,1-4+0,indvars(3))*q(i+1+1,1-4+0,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dxdx_0_2_ip11m4p0k = -&
          0.5_wp*d2_rhs_et_dxdx_0_2_ip11m4p0k_ip1m11m4p0k+&
          0.5_wp*d2_rhs_et_dxdx_0_2_ip11m4p0k_ip1p11m4p0k

d2_rhs_et_dxdx_0_2_ip11m4p0k = d2_rhs_et_dxdx_0_2_ip11m4p0k*param_float(1)

d1_rhs_et_dx_1_im11m4p0k = (q(i-1,1-4+0,indvars(1))*q(i-1,1-4+0,indvars(4))+&
                    (param_float(3 + 5))*q(i-1,1-4+0,indvars(1))*((q(i-1,1-4+0,indvars(4))-&
                    0.5_wp*(q(i-1,1-4+0,indvars(2))*q(i-1,1-4+0,indvars(2))+&
                    q(i-1,1-4+0,indvars(3))*q(i-1,1-4+0,indvars(3))))))*q(i-1,1-4+0,indvars(2))

d1_rhs_et_dx_1_ip11m4p0k = (q(i+1,1-4+0,indvars(1))*q(i+1,1-4+0,indvars(4))+&
                    (param_float(3 + 5))*q(i+1,1-4+0,indvars(1))*((q(i+1,1-4+0,indvars(4))-&
                    0.5_wp*(q(i+1,1-4+0,indvars(2))*q(i+1,1-4+0,indvars(2))+&
                    q(i+1,1-4+0,indvars(3))*q(i+1,1-4+0,indvars(3))))))*q(i+1,1-4+0,indvars(2))

d1_rhs_et_dx_1_i1m4p0k = -&
          0.5_wp*d1_rhs_et_dx_1_im11m4p0k+&
          0.5_wp*d1_rhs_et_dx_1_ip11m4p0k

d1_rhs_et_dx_1_i1m4p0k = d1_rhs_et_dx_1_i1m4p0k*param_float(1)

d1_rhs_et_dx_0_im11m4p0k = -4.0_wp/3.0_wp*(1.0_wp)*(1.0_wp+&
                    ((q(i-1,1-4+0,indvars(5))/1.0_wp*q(i-1,1-4+0,indvars(1)))**3.0_wp/((q(i-1,1-4+0,indvars(5))/1.0_wp*q(i-1,1-4+0,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i-1,1-4+0,indvars(5))/1.0_wp*q(i-1,1-4+0,indvars(1))))*param_float(1 + 5)*(d2_rhs_u_dxdx_0_0_im11m4p0k)*qst(i-1,1-4+0,indvarsst(10))*q(i-1,1-4+0,indvars(2))-&
                    (1.0_wp)*(1.0_wp+&
                    ((q(i-1,1-4+0,indvars(5))/1.0_wp*q(i-1,1-4+0,indvars(1)))**3.0_wp/((q(i-1,1-4+0,indvars(5))/1.0_wp*q(i-1,1-4+0,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i-1,1-4+0,indvars(5))/1.0_wp*q(i-1,1-4+0,indvars(1))))*param_float(1 + 5)*(d2_rhs_v_dxdx_0_0_im11m4p0k)*qst(i-1,1-4+0,indvarsst(10))*q(i-1,1-4+0,indvars(3))-&
                    param_float(2 + 5)*(d2_rhs_et_dxdx_0_2_im11m4p0k)*qst(i-1,1-4+0,indvarsst(10))

d1_rhs_et_dx_0_ip11m4p0k = -4.0_wp/3.0_wp*(1.0_wp)*(1.0_wp+&
                    ((q(i+1,1-4+0,indvars(5))/1.0_wp*q(i+1,1-4+0,indvars(1)))**3.0_wp/((q(i+1,1-4+0,indvars(5))/1.0_wp*q(i+1,1-4+0,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i+1,1-4+0,indvars(5))/1.0_wp*q(i+1,1-4+0,indvars(1))))*param_float(1 + 5)*(d2_rhs_u_dxdx_0_0_ip11m4p0k)*qst(i+1,1-4+0,indvarsst(10))*q(i+1,1-4+0,indvars(2))-&
                    (1.0_wp)*(1.0_wp+&
                    ((q(i+1,1-4+0,indvars(5))/1.0_wp*q(i+1,1-4+0,indvars(1)))**3.0_wp/((q(i+1,1-4+0,indvars(5))/1.0_wp*q(i+1,1-4+0,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i+1,1-4+0,indvars(5))/1.0_wp*q(i+1,1-4+0,indvars(1))))*param_float(1 + 5)*(d2_rhs_v_dxdx_0_0_ip11m4p0k)*qst(i+1,1-4+0,indvarsst(10))*q(i+1,1-4+0,indvars(3))-&
                    param_float(2 + 5)*(d2_rhs_et_dxdx_0_2_ip11m4p0k)*qst(i+1,1-4+0,indvarsst(10))

d1_rhs_et_dx_0_i1m4p0k = -&
          0.5_wp*d1_rhs_et_dx_0_im11m4p0k+&
          0.5_wp*d1_rhs_et_dx_0_ip11m4p0k

d1_rhs_et_dx_0_i1m4p0k = d1_rhs_et_dx_0_i1m4p0k*param_float(1)

d1_rhs_et_dx_2_im11m4p0k = (d2_rhs_et_dxdx_0_2_im11m4p0k)*qst(i-1,1-4+0,indvarsst(10))

d1_rhs_et_dx_2_ip11m4p0k = (d2_rhs_et_dxdx_0_2_ip11m4p0k)*qst(i+1,1-4+0,indvarsst(10))

d1_rhs_et_dx_2_i1m4p0k = -&
          0.5_wp*d1_rhs_et_dx_2_im11m4p0k+&
          0.5_wp*d1_rhs_et_dx_2_ip11m4p0k

d1_rhs_et_dx_2_i1m4p0k = d1_rhs_et_dx_2_i1m4p0k*param_float(1)

d1_rhs_et_dy_0_i1m4p0p0k = q(i,1-4+0+0,indvars(2))

d1_rhs_et_dy_0_i1m4p0p1k = q(i,1-4+0+1,indvars(2))

d1_rhs_et_dy_0_i1m4p0p2k = q(i,1-4+0+2,indvars(2))

d1_rhs_et_dy_0_i1m4p0k = -&
          1.5_wp*d1_rhs_et_dy_0_i1m4p0p0k+&
          2.0_wp*d1_rhs_et_dy_0_i1m4p0p1k-&
          0.5_wp*d1_rhs_et_dy_0_i1m4p0p2k

d1_rhs_et_dy_0_i1m4p0k = d1_rhs_et_dy_0_i1m4p0k*param_float(2)

d2_rhs_et_dydy_4_0_i1m4p0p0k_i1m4p0p0p0k = ((q(i,1-4+0+0+0,indvars(4))-&
                    0.5_wp*(q(i,1-4+0+0+0,indvars(2))*q(i,1-4+0+0+0,indvars(2))+&
                    q(i,1-4+0+0+0,indvars(3))*q(i,1-4+0+0+0,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_4_0_i1m4p0p0k_i1m4p0p0p1k = ((q(i,1-4+0+0+1,indvars(4))-&
                    0.5_wp*(q(i,1-4+0+0+1,indvars(2))*q(i,1-4+0+0+1,indvars(2))+&
                    q(i,1-4+0+0+1,indvars(3))*q(i,1-4+0+0+1,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_4_0_i1m4p0p0k_i1m4p0p0p2k = ((q(i,1-4+0+0+2,indvars(4))-&
                    0.5_wp*(q(i,1-4+0+0+2,indvars(2))*q(i,1-4+0+0+2,indvars(2))+&
                    q(i,1-4+0+0+2,indvars(3))*q(i,1-4+0+0+2,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_4_0_i1m4p0p0k = -&
          1.5_wp*d2_rhs_et_dydy_4_0_i1m4p0p0k_i1m4p0p0p0k+&
          2.0_wp*d2_rhs_et_dydy_4_0_i1m4p0p0k_i1m4p0p0p1k-&
          0.5_wp*d2_rhs_et_dydy_4_0_i1m4p0p0k_i1m4p0p0p2k

d2_rhs_et_dydy_4_0_i1m4p0p0k = d2_rhs_et_dydy_4_0_i1m4p0p0k*param_float(2)

d2_rhs_et_dydy_4_0_i1m4p0p1k_i1m4p0p1m1k = ((q(i,1-4+0+1-1,indvars(4))-&
                    0.5_wp*(q(i,1-4+0+1-1,indvars(2))*q(i,1-4+0+1-1,indvars(2))+&
                    q(i,1-4+0+1-1,indvars(3))*q(i,1-4+0+1-1,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_4_0_i1m4p0p1k_i1m4p0p1p1k = ((q(i,1-4+0+1+1,indvars(4))-&
                    0.5_wp*(q(i,1-4+0+1+1,indvars(2))*q(i,1-4+0+1+1,indvars(2))+&
                    q(i,1-4+0+1+1,indvars(3))*q(i,1-4+0+1+1,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_4_0_i1m4p0p1k = -&
          0.5_wp*d2_rhs_et_dydy_4_0_i1m4p0p1k_i1m4p0p1m1k+&
          0.5_wp*d2_rhs_et_dydy_4_0_i1m4p0p1k_i1m4p0p1p1k

d2_rhs_et_dydy_4_0_i1m4p0p1k = d2_rhs_et_dydy_4_0_i1m4p0p1k*param_float(2)

d2_rhs_et_dydy_4_0_i1m4p0p2k_i1m4p0p2m1k = ((q(i,1-4+0+2-1,indvars(4))-&
                    0.5_wp*(q(i,1-4+0+2-1,indvars(2))*q(i,1-4+0+2-1,indvars(2))+&
                    q(i,1-4+0+2-1,indvars(3))*q(i,1-4+0+2-1,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_4_0_i1m4p0p2k_i1m4p0p2p1k = ((q(i,1-4+0+2+1,indvars(4))-&
                    0.5_wp*(q(i,1-4+0+2+1,indvars(2))*q(i,1-4+0+2+1,indvars(2))+&
                    q(i,1-4+0+2+1,indvars(3))*q(i,1-4+0+2+1,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_4_0_i1m4p0p2k = -&
          0.5_wp*d2_rhs_et_dydy_4_0_i1m4p0p2k_i1m4p0p2m1k+&
          0.5_wp*d2_rhs_et_dydy_4_0_i1m4p0p2k_i1m4p0p2p1k

d2_rhs_et_dydy_4_0_i1m4p0p2k = d2_rhs_et_dydy_4_0_i1m4p0p2k*param_float(2)

d1_rhs_et_dy_4_i1m4p0p0k = (d2_rhs_et_dydy_4_0_i1m4p0p0k)*qst(i,1-4+0+0,indvarsst(11))

d1_rhs_et_dy_4_i1m4p0p1k = (d2_rhs_et_dydy_4_0_i1m4p0p1k)*qst(i,1-4+0+1,indvarsst(11))

d1_rhs_et_dy_4_i1m4p0p2k = (d2_rhs_et_dydy_4_0_i1m4p0p2k)*qst(i,1-4+0+2,indvarsst(11))

d1_rhs_et_dy_4_i1m4p0k = -&
          1.5_wp*d1_rhs_et_dy_4_i1m4p0p0k+&
          2.0_wp*d1_rhs_et_dy_4_i1m4p0p1k-&
          0.5_wp*d1_rhs_et_dy_4_i1m4p0p2k

d1_rhs_et_dy_4_i1m4p0k = d1_rhs_et_dy_4_i1m4p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None d(rho et)/dt *******
!                                                           
!***********************************************************


rhs(i,1-4+0,indvars(4)) =   -  ( (qst(i,1-4+0,indvarsst(10))*(d1_rhs_et_dx_0_i1m4p0k)+&
                    qst(i,1-4+0,indvarsst(10))*(d1_rhs_et_dx_1_i1m4p0k))*qst(i,1-4+0,indvarsst(5))+&
                    (-&
                    (1.0_wp)*(1.0_wp+&
                    ((q(i,1-4+0,indvars(5))/1.0_wp*q(i,1-4+0,indvars(1)))**3.0_wp/((q(i,1-4+0,indvars(5))/1.0_wp*q(i,1-4+0,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i,1-4+0,indvars(5))/1.0_wp*q(i,1-4+0,indvars(1))))*param_float(1 + 5)*(d1_rhs_et_dy_0_i1m4p0k)*(d1_rhs_et_dy_0_i1m4p0k)*qst(i,1-4+0,indvarsst(11))**2-&
                    4.0_wp/3.0_wp*(1.0_wp)*(1.0_wp+&
                    ((q(i,1-4+0,indvars(5))/1.0_wp*q(i,1-4+0,indvars(1)))**3.0_wp/((q(i,1-4+0,indvars(5))/1.0_wp*q(i,1-4+0,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i,1-4+0,indvars(5))/1.0_wp*q(i,1-4+0,indvars(1))))*param_float(1 + 5)*(d1_rhs_rho_dy_0_i1m4p0k)*(d1_rhs_rho_dy_0_i1m4p0k)*qst(i,1-4+0,indvarsst(11))**2-&
                    param_float(2 + 5)*(d1_rhs_et_dy_4_i1m4p0k)-&
                    param_float(2 + 5)*(d1_rhs_et_dx_2_i1m4p0k)+&
                    (q(i,1-4+0,indvars(1))*q(i,1-4+0,indvars(4))+&
                    (param_float(3 + 5))*q(i,1-4+0,indvars(1))*((q(i,1-4+0,indvars(4))-&
                    0.5_wp*(q(i,1-4+0,indvars(2))*q(i,1-4+0,indvars(2))+&
                    q(i,1-4+0,indvars(3))*q(i,1-4+0,indvars(3))))))*(d1_rhs_rho_dy_0_i1m4p0k)*qst(i,1-4+0,indvarsst(11)))*dabs(1.0_wp-&
                    qst(i,1-4+0,indvarsst(5))) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None d(rho nut)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (-ReI*Cb2*sigmaI*((deltaxI)**2*([rho*nut]_1x)*([nut]_1x))-Cb1*(1-ft2)*SS*rho*nut+ReI*(Cw1*fw-Cb1/k**2*ft2)*rho*nut**2/eta**2+deltaxI*([rho*u*nut-sigmaI*(visc+rho*nut)*({nut}_1x)*deltaxI]_1x))*symm
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_rhs_nut_dx_0_im11m4p0k = q(i-1,1-4+0,indvars(1))*q(i-1,1-4+0,indvars(5))

d1_rhs_nut_dx_0_ip11m4p0k = q(i+1,1-4+0,indvars(1))*q(i+1,1-4+0,indvars(5))

d1_rhs_nut_dx_0_i1m4p0k = -&
          0.5_wp*d1_rhs_nut_dx_0_im11m4p0k+&
          0.5_wp*d1_rhs_nut_dx_0_ip11m4p0k

d1_rhs_nut_dx_0_i1m4p0k = d1_rhs_nut_dx_0_i1m4p0k*param_float(1)

d1_rhs_nut_dx_1_im11m4p0k = q(i-1,1-4+0,indvars(5))

d1_rhs_nut_dx_1_ip11m4p0k = q(i+1,1-4+0,indvars(5))

d1_rhs_nut_dx_1_i1m4p0k = -&
          0.5_wp*d1_rhs_nut_dx_1_im11m4p0k+&
          0.5_wp*d1_rhs_nut_dx_1_ip11m4p0k

d1_rhs_nut_dx_1_i1m4p0k = d1_rhs_nut_dx_1_i1m4p0k*param_float(1)

d2_rhs_nut_dxdx_2_0_im11m4p0k_im1m11m4p0k = q(i-1-1,1-4+0,indvars(5))

d2_rhs_nut_dxdx_2_0_im11m4p0k_im1p11m4p0k = q(i-1+1,1-4+0,indvars(5))

d2_rhs_nut_dxdx_2_0_im11m4p0k = -&
          0.5_wp*d2_rhs_nut_dxdx_2_0_im11m4p0k_im1m11m4p0k+&
          0.5_wp*d2_rhs_nut_dxdx_2_0_im11m4p0k_im1p11m4p0k

d2_rhs_nut_dxdx_2_0_im11m4p0k = d2_rhs_nut_dxdx_2_0_im11m4p0k*param_float(1)

d2_rhs_nut_dxdx_2_0_ip11m4p0k_ip1m11m4p0k = q(i+1-1,1-4+0,indvars(5))

d2_rhs_nut_dxdx_2_0_ip11m4p0k_ip1p11m4p0k = q(i+1+1,1-4+0,indvars(5))

d2_rhs_nut_dxdx_2_0_ip11m4p0k = -&
          0.5_wp*d2_rhs_nut_dxdx_2_0_ip11m4p0k_ip1m11m4p0k+&
          0.5_wp*d2_rhs_nut_dxdx_2_0_ip11m4p0k_ip1p11m4p0k

d2_rhs_nut_dxdx_2_0_ip11m4p0k = d2_rhs_nut_dxdx_2_0_ip11m4p0k*param_float(1)

d1_rhs_nut_dx_2_im11m4p0k = q(i-1,1-4+0,indvars(1))*q(i-1,1-4+0,indvars(2))*q(i-1,1-4+0,indvars(5))-&
                    param_float(18 + 5)*(1.0_wp+&
                    q(i-1,1-4+0,indvars(1))*q(i-1,1-4+0,indvars(5)))*(d2_rhs_nut_dxdx_2_0_im11m4p0k)*qst(i-1,1-4+0,indvarsst(10))

d1_rhs_nut_dx_2_ip11m4p0k = q(i+1,1-4+0,indvars(1))*q(i+1,1-4+0,indvars(2))*q(i+1,1-4+0,indvars(5))-&
                    param_float(18 + 5)*(1.0_wp+&
                    q(i+1,1-4+0,indvars(1))*q(i+1,1-4+0,indvars(5)))*(d2_rhs_nut_dxdx_2_0_ip11m4p0k)*qst(i+1,1-4+0,indvarsst(10))

d1_rhs_nut_dx_2_i1m4p0k = -&
          0.5_wp*d1_rhs_nut_dx_2_im11m4p0k+&
          0.5_wp*d1_rhs_nut_dx_2_ip11m4p0k

d1_rhs_nut_dx_2_i1m4p0k = d1_rhs_nut_dx_2_i1m4p0k*param_float(1)



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None d(rho nut)/dt ******
!                                                           
!***********************************************************


rhs(i,1-4+0,indvars(5)) =   -  ( (-&
                    param_float(1 + 5)*param_float(7 + 5)*param_float(18 + 5)*((qst(i,1-4+0,indvarsst(10)))**2*(d1_rhs_nut_dx_0_i1m4p0k)*(d1_rhs_nut_dx_1_i1m4p0k))-&
                    param_float(6 + 5)*(1-&
                    (param_float(16 + 5)*exp(-&
                    param_float(17 + 5)*(q(i,1-4+0,indvars(5))/1.0_wp*q(i,1-4+0,indvars(1)))**2.0_wp)))*qst(i,1-4+0,indvarsst(12))*q(i,1-4+0,indvars(1))*q(i,1-4+0,indvars(5))+&
                    param_float(1 + 5)*(param_float(10 + 5)*(((min(param_float(1 + 5)*(q(i,1-4+0,indvars(5))/(qst(i,1-4+0,indvarsst(12))*param_float(9 + 5)**2.0_wp*qst(i,1-4+0,indvarsst(2))**2.0_wp)),10.0_wp))+&
                    param_float(11 + 5)*((min(param_float(1 + 5)*(q(i,1-4+0,indvars(5))/(qst(i,1-4+0,indvarsst(12))*param_float(9 + 5)**2.0_wp*qst(i,1-4+0,indvarsst(2))**2.0_wp)),10.0_wp))**6.0_wp-&
                    (min(param_float(1 + 5)*(q(i,1-4+0,indvars(5))/(qst(i,1-4+0,indvarsst(12))*param_float(9 + 5)**2.0_wp*qst(i,1-4+0,indvarsst(2))**2.0_wp)),10.0_wp))))*((1.0_wp+&
                    param_float(12 + 5)**6.0_wp)/(((min(param_float(1 + 5)*(q(i,1-4+0,indvars(5))/(qst(i,1-4+0,indvarsst(12))*param_float(9 + 5)**2.0_wp*qst(i,1-4+0,indvarsst(2))**2.0_wp)),10.0_wp))+&
                    param_float(11 + 5)*((min(param_float(1 + 5)*(q(i,1-4+0,indvars(5))/(qst(i,1-4+0,indvarsst(12))*param_float(9 + 5)**2.0_wp*qst(i,1-4+0,indvarsst(2))**2.0_wp)),10.0_wp))**6.0_wp-&
                    (min(param_float(1 + 5)*(q(i,1-4+0,indvars(5))/(qst(i,1-4+0,indvarsst(12))*param_float(9 + 5)**2.0_wp*qst(i,1-4+0,indvarsst(2))**2.0_wp)),10.0_wp))))**6.0_wp+&
                    param_float(12 + 5)**6.0_wp))**(1.0_wp/6.0_wp))-&
                    param_float(6 + 5)/param_float(9 + 5)**2*(param_float(16 + 5)*exp(-&
                    param_float(17 + 5)*(q(i,1-4+0,indvars(5))/1.0_wp*q(i,1-4+0,indvars(1)))**2.0_wp)))*q(i,1-4+0,indvars(1))*q(i,1-4+0,indvars(5))**2/qst(i,1-4+0,indvarsst(2))**2+&
                    qst(i,1-4+0,indvarsst(10))*(d1_rhs_nut_dx_2_i1m4p0k))*qst(i,1-4+0,indvarsst(5)) ) 

   enddo

    enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO

 end subroutine phyLow_surfbc_faces_j1_6


