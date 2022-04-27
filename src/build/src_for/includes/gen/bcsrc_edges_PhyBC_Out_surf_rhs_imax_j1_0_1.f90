














!====================================================================================================
!
! General Boundary Conditions: compute Eqns with Boundary Scheme
!     
!====================================================================================================

subroutine boundaryScheme_edges_PhyBC_Out_surf_rhs_imax_j1_0_1(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

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



 real(wp) ::  d1_rhs_rho_dx_0_nxp4p0p01m4p1k,d1_rhs_rho_dx_0_nxp4p0m11m4p1k,d1_rhs_rho_dx_0_nxp4p0m21m4p1k &
            ,d1_rhs_rho_dx_0_nxp4p01m4p1k &
            ,d1_rhs_rho_dx_1_nxp4p0p01m4p1k,d1_rhs_rho_dx_1_nxp4p0m11m4p1k,d1_rhs_rho_dx_1_nxp4p0m21m4p1k &
            ,d1_rhs_rho_dx_1_nxp4p01m4p1k &
            ,d1_rhs_rho_dx_2_nxp4p0p01m4p1k,d1_rhs_rho_dx_2_nxp4p0m11m4p1k,d1_rhs_rho_dx_2_nxp4p0m21m4p1k &
            ,d1_rhs_rho_dx_2_nxp4p01m4p1k &
            ,d1_rhs_rho_dx_3_nxp4p0p01m4p1k,d1_rhs_rho_dx_3_nxp4p0m11m4p1k,d1_rhs_rho_dx_3_nxp4p0m21m4p1k &
            ,d1_rhs_rho_dx_3_nxp4p01m4p1k &
            ,d1_rhs_rho_dy_0_nxp4p01m4p1m1k,d1_rhs_rho_dy_0_nxp4p01m4p1p0k,d1_rhs_rho_dy_0_nxp4p01m4p1p1k &
            ,d1_rhs_rho_dy_0_nxp4p01m4p1k &
            ,d2_rhs_u_dxdx_6_0_nxp4p0p01m4p1k_nxp4p0p0p01m4p1k,d2_rhs_u_dxdx_6_0_nxp4p0p01m4p1k_nxp4p0p0m11m4p1k,d2_rhs_u_dxdx_6_0_nxp4p0p01m4p1k_nxp4p0p0m21m4p1k &
            ,d2_rhs_u_dxdx_6_0_nxp4p0p01m4p1k &
            ,d2_rhs_u_dxdx_6_0_nxp4p0m11m4p1k_nxp4p0m1m11m4p1k,d2_rhs_u_dxdx_6_0_nxp4p0m11m4p1k_nxp4p0m1p01m4p1k,d2_rhs_u_dxdx_6_0_nxp4p0m11m4p1k_nxp4p0m1p11m4p1k &
            ,d2_rhs_u_dxdx_6_0_nxp4p0m11m4p1k &
            ,d2_rhs_u_dxdx_6_0_nxp4p0m21m4p1k_nxp4p0m2m11m4p1k,d2_rhs_u_dxdx_6_0_nxp4p0m21m4p1k_nxp4p0m2p01m4p1k,d2_rhs_u_dxdx_6_0_nxp4p0m21m4p1k_nxp4p0m2p11m4p1k &
            ,d2_rhs_u_dxdx_6_0_nxp4p0m21m4p1k &
            ,d2_rhs_u_dxdy_6_0_nxp4p0p01m4p1k_nxp4p0p01m4p1m1k,d2_rhs_u_dxdy_6_0_nxp4p0p01m4p1k_nxp4p0p01m4p1p0k,d2_rhs_u_dxdy_6_0_nxp4p0p01m4p1k_nxp4p0p01m4p1p1k &
            ,d2_rhs_u_dxdy_6_0_nxp4p0p01m4p1k &
            ,d2_rhs_u_dxdy_6_0_nxp4p0m11m4p1k_nxp4p0m11m4p1m1k,d2_rhs_u_dxdy_6_0_nxp4p0m11m4p1k_nxp4p0m11m4p1p0k,d2_rhs_u_dxdy_6_0_nxp4p0m11m4p1k_nxp4p0m11m4p1p1k &
            ,d2_rhs_u_dxdy_6_0_nxp4p0m11m4p1k &
            ,d2_rhs_u_dxdy_6_0_nxp4p0m21m4p1k_nxp4p0m21m4p1m1k,d2_rhs_u_dxdy_6_0_nxp4p0m21m4p1k_nxp4p0m21m4p1p0k,d2_rhs_u_dxdy_6_0_nxp4p0m21m4p1k_nxp4p0m21m4p1p1k &
            ,d2_rhs_u_dxdy_6_0_nxp4p0m21m4p1k &
            ,d1_rhs_u_dx_6_nxp4p0p01m4p1k,d1_rhs_u_dx_6_nxp4p0m11m4p1k,d1_rhs_u_dx_6_nxp4p0m21m4p1k &
            ,d1_rhs_u_dx_6_nxp4p01m4p1k &
            ,d1_rhs_u_dy_0_nxp4p01m4p1m1k,d1_rhs_u_dy_0_nxp4p01m4p1p0k,d1_rhs_u_dy_0_nxp4p01m4p1p1k &
            ,d1_rhs_u_dy_0_nxp4p01m4p1k &
            ,d2_rhs_u_dydx_1_0_nxp4p01m4p1m1k_nxp4p0p01m4p1m1k,d2_rhs_u_dydx_1_0_nxp4p01m4p1m1k_nxp4p0m11m4p1m1k,d2_rhs_u_dydx_1_0_nxp4p01m4p1m1k_nxp4p0m21m4p1m1k &
            ,d2_rhs_u_dydx_1_0_nxp4p01m4p1m1k &
            ,d2_rhs_u_dydx_1_0_nxp4p01m4p1p1k_nxp4p0p01m4p1p1k,d2_rhs_u_dydx_1_0_nxp4p01m4p1p1k_nxp4p0m11m4p1p1k,d2_rhs_u_dydx_1_0_nxp4p01m4p1p1k_nxp4p0m21m4p1p1k &
            ,d2_rhs_u_dydx_1_0_nxp4p01m4p1p1k &
            ,d2_rhs_u_dydy_1_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p0k,d2_rhs_u_dydy_1_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p1k,d2_rhs_u_dydy_1_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p2k &
            ,d2_rhs_u_dydy_1_0_nxp4p01m4p1m1k &
            ,d2_rhs_u_dydy_1_0_nxp4p01m4p1p1k_nxp4p01m4p1p1m1k,d2_rhs_u_dydy_1_0_nxp4p01m4p1p1k_nxp4p01m4p1p1p0k,d2_rhs_u_dydy_1_0_nxp4p01m4p1p1k_nxp4p01m4p1p1p1k &
            ,d2_rhs_u_dydy_1_0_nxp4p01m4p1p1k &
            ,d1_rhs_u_dy_1_nxp4p01m4p1m1k,d1_rhs_u_dy_1_nxp4p01m4p1p0k,d1_rhs_u_dy_1_nxp4p01m4p1p1k &
            ,d1_rhs_u_dy_1_nxp4p01m4p1k &
            ,d1_rhs_v_dx_4_nxp4p0p01m4p1k,d1_rhs_v_dx_4_nxp4p0m11m4p1k,d1_rhs_v_dx_4_nxp4p0m21m4p1k &
            ,d1_rhs_v_dx_4_nxp4p01m4p1k &
            ,d2_rhs_v_dxdx_5_0_nxp4p0p01m4p1k_nxp4p0p0p01m4p1k,d2_rhs_v_dxdx_5_0_nxp4p0p01m4p1k_nxp4p0p0m11m4p1k,d2_rhs_v_dxdx_5_0_nxp4p0p01m4p1k_nxp4p0p0m21m4p1k &
            ,d2_rhs_v_dxdx_5_0_nxp4p0p01m4p1k &
            ,d2_rhs_v_dxdx_5_0_nxp4p0m11m4p1k_nxp4p0m1m11m4p1k,d2_rhs_v_dxdx_5_0_nxp4p0m11m4p1k_nxp4p0m1p01m4p1k,d2_rhs_v_dxdx_5_0_nxp4p0m11m4p1k_nxp4p0m1p11m4p1k &
            ,d2_rhs_v_dxdx_5_0_nxp4p0m11m4p1k &
            ,d2_rhs_v_dxdx_5_0_nxp4p0m21m4p1k_nxp4p0m2m11m4p1k,d2_rhs_v_dxdx_5_0_nxp4p0m21m4p1k_nxp4p0m2p01m4p1k,d2_rhs_v_dxdx_5_0_nxp4p0m21m4p1k_nxp4p0m2p11m4p1k &
            ,d2_rhs_v_dxdx_5_0_nxp4p0m21m4p1k &
            ,d2_rhs_v_dxdy_5_0_nxp4p0p01m4p1k_nxp4p0p01m4p1m1k,d2_rhs_v_dxdy_5_0_nxp4p0p01m4p1k_nxp4p0p01m4p1p0k,d2_rhs_v_dxdy_5_0_nxp4p0p01m4p1k_nxp4p0p01m4p1p1k &
            ,d2_rhs_v_dxdy_5_0_nxp4p0p01m4p1k &
            ,d2_rhs_v_dxdy_5_0_nxp4p0m11m4p1k_nxp4p0m11m4p1m1k,d2_rhs_v_dxdy_5_0_nxp4p0m11m4p1k_nxp4p0m11m4p1p0k,d2_rhs_v_dxdy_5_0_nxp4p0m11m4p1k_nxp4p0m11m4p1p1k &
            ,d2_rhs_v_dxdy_5_0_nxp4p0m11m4p1k &
            ,d2_rhs_v_dxdy_5_0_nxp4p0m21m4p1k_nxp4p0m21m4p1m1k,d2_rhs_v_dxdy_5_0_nxp4p0m21m4p1k_nxp4p0m21m4p1p0k,d2_rhs_v_dxdy_5_0_nxp4p0m21m4p1k_nxp4p0m21m4p1p1k &
            ,d2_rhs_v_dxdy_5_0_nxp4p0m21m4p1k &
            ,d1_rhs_v_dx_5_nxp4p0p01m4p1k,d1_rhs_v_dx_5_nxp4p0m11m4p1k,d1_rhs_v_dx_5_nxp4p0m21m4p1k &
            ,d1_rhs_v_dx_5_nxp4p01m4p1k &
            ,d1_rhs_v_dy_0_nxp4p01m4p1m1k,d1_rhs_v_dy_0_nxp4p01m4p1p0k,d1_rhs_v_dy_0_nxp4p01m4p1p1k &
            ,d1_rhs_v_dy_0_nxp4p01m4p1k &
            ,d2_rhs_v_dydx_1_0_nxp4p01m4p1m1k_nxp4p0p01m4p1m1k,d2_rhs_v_dydx_1_0_nxp4p01m4p1m1k_nxp4p0m11m4p1m1k,d2_rhs_v_dydx_1_0_nxp4p01m4p1m1k_nxp4p0m21m4p1m1k &
            ,d2_rhs_v_dydx_1_0_nxp4p01m4p1m1k &
            ,d2_rhs_v_dydx_1_0_nxp4p01m4p1p1k_nxp4p0p01m4p1p1k,d2_rhs_v_dydx_1_0_nxp4p01m4p1p1k_nxp4p0m11m4p1p1k,d2_rhs_v_dydx_1_0_nxp4p01m4p1p1k_nxp4p0m21m4p1p1k &
            ,d2_rhs_v_dydx_1_0_nxp4p01m4p1p1k &
            ,d2_rhs_v_dydy_1_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p0k,d2_rhs_v_dydy_1_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p1k,d2_rhs_v_dydy_1_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p2k &
            ,d2_rhs_v_dydy_1_0_nxp4p01m4p1m1k &
            ,d2_rhs_v_dydy_1_0_nxp4p01m4p1p1k_nxp4p01m4p1p1m1k,d2_rhs_v_dydy_1_0_nxp4p01m4p1p1k_nxp4p01m4p1p1p0k,d2_rhs_v_dydy_1_0_nxp4p01m4p1p1k_nxp4p01m4p1p1p1k &
            ,d2_rhs_v_dydy_1_0_nxp4p01m4p1p1k &
            ,d1_rhs_v_dy_1_nxp4p01m4p1m1k,d1_rhs_v_dy_1_nxp4p01m4p1p0k,d1_rhs_v_dy_1_nxp4p01m4p1p1k &
            ,d1_rhs_v_dy_1_nxp4p01m4p1k &
            ,d2_rhs_et_dxdx_9_0_nxp4p0p01m4p1k_nxp4p0p0p01m4p1k,d2_rhs_et_dxdx_9_0_nxp4p0p01m4p1k_nxp4p0p0m11m4p1k,d2_rhs_et_dxdx_9_0_nxp4p0p01m4p1k_nxp4p0p0m21m4p1k &
            ,d2_rhs_et_dxdx_9_0_nxp4p0p01m4p1k &
            ,d2_rhs_et_dxdx_9_0_nxp4p0m11m4p1k_nxp4p0m1m11m4p1k,d2_rhs_et_dxdx_9_0_nxp4p0m11m4p1k_nxp4p0m1p01m4p1k,d2_rhs_et_dxdx_9_0_nxp4p0m11m4p1k_nxp4p0m1p11m4p1k &
            ,d2_rhs_et_dxdx_9_0_nxp4p0m11m4p1k &
            ,d2_rhs_et_dxdx_9_0_nxp4p0m21m4p1k_nxp4p0m2m11m4p1k,d2_rhs_et_dxdx_9_0_nxp4p0m21m4p1k_nxp4p0m2p01m4p1k,d2_rhs_et_dxdx_9_0_nxp4p0m21m4p1k_nxp4p0m2p11m4p1k &
            ,d2_rhs_et_dxdx_9_0_nxp4p0m21m4p1k &
            ,d1_rhs_et_dx_9_nxp4p0p01m4p1k,d1_rhs_et_dx_9_nxp4p0m11m4p1k,d1_rhs_et_dx_9_nxp4p0m21m4p1k &
            ,d1_rhs_et_dx_9_nxp4p01m4p1k &
            ,d1_rhs_et_dy_0_nxp4p01m4p1m1k,d1_rhs_et_dy_0_nxp4p01m4p1p0k,d1_rhs_et_dy_0_nxp4p01m4p1p1k &
            ,d1_rhs_et_dy_0_nxp4p01m4p1k &
            ,d2_rhs_et_dydy_1_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p0k,d2_rhs_et_dydy_1_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p1k,d2_rhs_et_dydy_1_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p2k &
            ,d2_rhs_et_dydy_1_0_nxp4p01m4p1m1k &
            ,d2_rhs_et_dydy_1_0_nxp4p01m4p1p1k_nxp4p01m4p1p1m1k,d2_rhs_et_dydy_1_0_nxp4p01m4p1p1k_nxp4p01m4p1p1p0k,d2_rhs_et_dydy_1_0_nxp4p01m4p1p1k_nxp4p01m4p1p1p1k &
            ,d2_rhs_et_dydy_1_0_nxp4p01m4p1p1k &
            ,d1_rhs_et_dy_1_nxp4p01m4p1m1k,d1_rhs_et_dy_1_nxp4p01m4p1p0k,d1_rhs_et_dy_1_nxp4p01m4p1p1k &
            ,d1_rhs_et_dy_1_nxp4p01m4p1k &
            ,d2_rhs_nut_dydy_0_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p0k,d2_rhs_nut_dydy_0_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p1k,d2_rhs_nut_dydy_0_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p2k &
            ,d2_rhs_nut_dydy_0_0_nxp4p01m4p1m1k &
            ,d2_rhs_nut_dydy_0_0_nxp4p01m4p1p1k_nxp4p01m4p1p1m1k,d2_rhs_nut_dydy_0_0_nxp4p01m4p1p1k_nxp4p01m4p1p1p0k,d2_rhs_nut_dydy_0_0_nxp4p01m4p1p1k_nxp4p01m4p1p1p1k &
            ,d2_rhs_nut_dydy_0_0_nxp4p01m4p1p1k &
            ,d1_rhs_nut_dy_1_nxp4p01m4p1m1k,d1_rhs_nut_dy_1_nxp4p01m4p1p0k,d1_rhs_nut_dy_1_nxp4p01m4p1p1k &
            ,d1_rhs_nut_dy_1_nxp4p01m4p1k &
            ,d1_rhs_nut_dy_2_nxp4p01m4p1m1k,d1_rhs_nut_dy_2_nxp4p01m4p1p0k,d1_rhs_nut_dy_2_nxp4p01m4p1p1k &
            ,d1_rhs_nut_dy_2_nxp4p01m4p1k &
            ,d1_rhs_nut_dy_0_nxp4p01m4p1m1k,d1_rhs_nut_dy_0_nxp4p01m4p1p0k,d1_rhs_nut_dy_0_nxp4p01m4p1p1k &
            ,d1_rhs_nut_dy_0_nxp4p01m4p1k 

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
! Start building layers for BC : imax j1 None **************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 0 1 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None d(rho)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([rho*u]_1x)+deltayI*([rho*v]_1y)+(1.0_wp/c**2*((u*(-c**2*[rho]_1x+1.0_wp*[p]_1x))*deltaxI+0.5_wp*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI+esse*c*(1-M_imax*M_imax)/L_ref*(p-P0))))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_rhs_rho_dx_0_nxp4p0p01m4p1k = q(nx+4+0+0,1-4+1,indvars(1))*q(nx+4+0+0,1-4+1,indvars(2))

d1_rhs_rho_dx_0_nxp4p0m11m4p1k = q(nx+4+0-1,1-4+1,indvars(1))*q(nx+4+0-1,1-4+1,indvars(2))

d1_rhs_rho_dx_0_nxp4p0m21m4p1k = q(nx+4+0-2,1-4+1,indvars(1))*q(nx+4+0-2,1-4+1,indvars(2))

d1_rhs_rho_dx_0_nxp4p01m4p1k = 1.5_wp*d1_rhs_rho_dx_0_nxp4p0p01m4p1k-&
          2.0_wp*d1_rhs_rho_dx_0_nxp4p0m11m4p1k+&
          0.5_wp*d1_rhs_rho_dx_0_nxp4p0m21m4p1k

d1_rhs_rho_dx_0_nxp4p01m4p1k = d1_rhs_rho_dx_0_nxp4p01m4p1k*param_float(1)

d1_rhs_rho_dx_1_nxp4p0p01m4p1k = q(nx+4+0+0,1-4+1,indvars(1))

d1_rhs_rho_dx_1_nxp4p0m11m4p1k = q(nx+4+0-1,1-4+1,indvars(1))

d1_rhs_rho_dx_1_nxp4p0m21m4p1k = q(nx+4+0-2,1-4+1,indvars(1))

d1_rhs_rho_dx_1_nxp4p01m4p1k = 1.5_wp*d1_rhs_rho_dx_1_nxp4p0p01m4p1k-&
          2.0_wp*d1_rhs_rho_dx_1_nxp4p0m11m4p1k+&
          0.5_wp*d1_rhs_rho_dx_1_nxp4p0m21m4p1k

d1_rhs_rho_dx_1_nxp4p01m4p1k = d1_rhs_rho_dx_1_nxp4p01m4p1k*param_float(1)

d1_rhs_rho_dx_2_nxp4p0p01m4p1k = (param_float(3 + 5))*q(nx+4+0+0,1-4+1,indvars(1))*((q(nx+4+0+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+1,indvars(2))*q(nx+4+0+0,1-4+1,indvars(2))+&
                    q(nx+4+0+0,1-4+1,indvars(3))*q(nx+4+0+0,1-4+1,indvars(3)))))

d1_rhs_rho_dx_2_nxp4p0m11m4p1k = (param_float(3 + 5))*q(nx+4+0-1,1-4+1,indvars(1))*((q(nx+4+0-1,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+1,indvars(2))*q(nx+4+0-1,1-4+1,indvars(2))+&
                    q(nx+4+0-1,1-4+1,indvars(3))*q(nx+4+0-1,1-4+1,indvars(3)))))

d1_rhs_rho_dx_2_nxp4p0m21m4p1k = (param_float(3 + 5))*q(nx+4+0-2,1-4+1,indvars(1))*((q(nx+4+0-2,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+1,indvars(2))*q(nx+4+0-2,1-4+1,indvars(2))+&
                    q(nx+4+0-2,1-4+1,indvars(3))*q(nx+4+0-2,1-4+1,indvars(3)))))

d1_rhs_rho_dx_2_nxp4p01m4p1k = 1.5_wp*d1_rhs_rho_dx_2_nxp4p0p01m4p1k-&
          2.0_wp*d1_rhs_rho_dx_2_nxp4p0m11m4p1k+&
          0.5_wp*d1_rhs_rho_dx_2_nxp4p0m21m4p1k

d1_rhs_rho_dx_2_nxp4p01m4p1k = d1_rhs_rho_dx_2_nxp4p01m4p1k*param_float(1)

d1_rhs_rho_dx_3_nxp4p0p01m4p1k = q(nx+4+0+0,1-4+1,indvars(2))

d1_rhs_rho_dx_3_nxp4p0m11m4p1k = q(nx+4+0-1,1-4+1,indvars(2))

d1_rhs_rho_dx_3_nxp4p0m21m4p1k = q(nx+4+0-2,1-4+1,indvars(2))

d1_rhs_rho_dx_3_nxp4p01m4p1k = 1.5_wp*d1_rhs_rho_dx_3_nxp4p0p01m4p1k-&
          2.0_wp*d1_rhs_rho_dx_3_nxp4p0m11m4p1k+&
          0.5_wp*d1_rhs_rho_dx_3_nxp4p0m21m4p1k

d1_rhs_rho_dx_3_nxp4p01m4p1k = d1_rhs_rho_dx_3_nxp4p01m4p1k*param_float(1)

d1_rhs_rho_dy_0_nxp4p01m4p1m1k = q(nx+4+0,1-4+1-1,indvars(1))*q(nx+4+0,1-4+1-1,indvars(3))

d1_rhs_rho_dy_0_nxp4p01m4p1p1k = q(nx+4+0,1-4+1+1,indvars(1))*q(nx+4+0,1-4+1+1,indvars(3))

d1_rhs_rho_dy_0_nxp4p01m4p1k = -&
          0.5_wp*d1_rhs_rho_dy_0_nxp4p01m4p1m1k+&
          0.5_wp*d1_rhs_rho_dy_0_nxp4p01m4p1p1k

d1_rhs_rho_dy_0_nxp4p01m4p1k = d1_rhs_rho_dy_0_nxp4p01m4p1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None d(rho)/dt *************
!                                                           
!***********************************************************


rhs(nx+4+0,1-4+1,indvars(1)) =   -  ( qst(nx+4+0,1-4+1,indvarsst(10))*(d1_rhs_rho_dx_0_nxp4p01m4p1k)+&
                    qst(nx+4+0,1-4+1,indvarsst(11))*(d1_rhs_rho_dy_0_nxp4p01m4p1k)+&
                    (1.0_wp/(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5**2*((q(nx+4+0,1-4+1,indvars(2))*(-&
                    (param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5**2*d1_rhs_rho_dx_1_nxp4p01m4p1k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p01m4p1k))*qst(nx+4+0,1-4+1,indvarsst(10))+&
                    0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5+&
                    q(nx+4+0,1-4+1,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5*q(nx+4+0,1-4+1,indvars(1))*d1_rhs_rho_dx_3_nxp4p01m4p1k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p01m4p1k))*qst(nx+4+0,1-4+1,indvarsst(10))+&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5*(1-&
                    qface_i(1-4+1,2)*qface_i(1-4+1,2))/param_float(20 + 5)*((param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))-&
                    param_float(26 + 5))))) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None d(rho u)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! u*(1.0_wp/c**2*((u*(-c**2*[rho]_1x+1.0_wp*[p]_1x))*deltaxI+0.5_wp*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI+esse*c*(1-M_imax*M_imax)/L_ref*(p-P0))))+rho*(1.0_wp/(2*rho*c)*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI-esse*c*(1-M_imax*M_imax)/L_ref*(p-P0)))+[rho*u*v]_1y+deltaxI*([-visc_t*(2.0_wp*deltaxI*({u}_1x)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y)))]_1x)+deltayI*([-visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_rhs_u_dxdx_6_0_nxp4p0p01m4p1k_nxp4p0p0p01m4p1k = q(nx+4+0+0+0,1-4+1,indvars(2))

d2_rhs_u_dxdx_6_0_nxp4p0p01m4p1k_nxp4p0p0m11m4p1k = q(nx+4+0+0-1,1-4+1,indvars(2))

d2_rhs_u_dxdx_6_0_nxp4p0p01m4p1k_nxp4p0p0m21m4p1k = q(nx+4+0+0-2,1-4+1,indvars(2))

d2_rhs_u_dxdx_6_0_nxp4p0p01m4p1k = 1.5_wp*d2_rhs_u_dxdx_6_0_nxp4p0p01m4p1k_nxp4p0p0p01m4p1k-&
          2.0_wp*d2_rhs_u_dxdx_6_0_nxp4p0p01m4p1k_nxp4p0p0m11m4p1k+&
          0.5_wp*d2_rhs_u_dxdx_6_0_nxp4p0p01m4p1k_nxp4p0p0m21m4p1k

d2_rhs_u_dxdx_6_0_nxp4p0p01m4p1k = d2_rhs_u_dxdx_6_0_nxp4p0p01m4p1k*param_float(1)

d2_rhs_u_dxdx_6_0_nxp4p0m11m4p1k_nxp4p0m1m11m4p1k = q(nx+4+0-1-1,1-4+1,indvars(2))

d2_rhs_u_dxdx_6_0_nxp4p0m11m4p1k_nxp4p0m1p11m4p1k = q(nx+4+0-1+1,1-4+1,indvars(2))

d2_rhs_u_dxdx_6_0_nxp4p0m11m4p1k = -&
          0.5_wp*d2_rhs_u_dxdx_6_0_nxp4p0m11m4p1k_nxp4p0m1m11m4p1k+&
          0.5_wp*d2_rhs_u_dxdx_6_0_nxp4p0m11m4p1k_nxp4p0m1p11m4p1k

d2_rhs_u_dxdx_6_0_nxp4p0m11m4p1k = d2_rhs_u_dxdx_6_0_nxp4p0m11m4p1k*param_float(1)

d2_rhs_u_dxdx_6_0_nxp4p0m21m4p1k_nxp4p0m2m11m4p1k = q(nx+4+0-2-1,1-4+1,indvars(2))

d2_rhs_u_dxdx_6_0_nxp4p0m21m4p1k_nxp4p0m2p11m4p1k = q(nx+4+0-2+1,1-4+1,indvars(2))

d2_rhs_u_dxdx_6_0_nxp4p0m21m4p1k = -&
          0.5_wp*d2_rhs_u_dxdx_6_0_nxp4p0m21m4p1k_nxp4p0m2m11m4p1k+&
          0.5_wp*d2_rhs_u_dxdx_6_0_nxp4p0m21m4p1k_nxp4p0m2p11m4p1k

d2_rhs_u_dxdx_6_0_nxp4p0m21m4p1k = d2_rhs_u_dxdx_6_0_nxp4p0m21m4p1k*param_float(1)

d2_rhs_u_dxdy_6_0_nxp4p0p01m4p1k_nxp4p0p01m4p1m1k = q(nx+4+0+0,1-4+1-1,indvars(3))

d2_rhs_u_dxdy_6_0_nxp4p0p01m4p1k_nxp4p0p01m4p1p1k = q(nx+4+0+0,1-4+1+1,indvars(3))

d2_rhs_u_dxdy_6_0_nxp4p0p01m4p1k = -&
          0.5_wp*d2_rhs_u_dxdy_6_0_nxp4p0p01m4p1k_nxp4p0p01m4p1m1k+&
          0.5_wp*d2_rhs_u_dxdy_6_0_nxp4p0p01m4p1k_nxp4p0p01m4p1p1k

d2_rhs_u_dxdy_6_0_nxp4p0p01m4p1k = d2_rhs_u_dxdy_6_0_nxp4p0p01m4p1k*param_float(2)

d2_rhs_u_dxdy_6_0_nxp4p0m11m4p1k_nxp4p0m11m4p1m1k = q(nx+4+0-1,1-4+1-1,indvars(3))

d2_rhs_u_dxdy_6_0_nxp4p0m11m4p1k_nxp4p0m11m4p1p1k = q(nx+4+0-1,1-4+1+1,indvars(3))

d2_rhs_u_dxdy_6_0_nxp4p0m11m4p1k = -&
          0.5_wp*d2_rhs_u_dxdy_6_0_nxp4p0m11m4p1k_nxp4p0m11m4p1m1k+&
          0.5_wp*d2_rhs_u_dxdy_6_0_nxp4p0m11m4p1k_nxp4p0m11m4p1p1k

d2_rhs_u_dxdy_6_0_nxp4p0m11m4p1k = d2_rhs_u_dxdy_6_0_nxp4p0m11m4p1k*param_float(2)

d2_rhs_u_dxdy_6_0_nxp4p0m21m4p1k_nxp4p0m21m4p1m1k = q(nx+4+0-2,1-4+1-1,indvars(3))

d2_rhs_u_dxdy_6_0_nxp4p0m21m4p1k_nxp4p0m21m4p1p1k = q(nx+4+0-2,1-4+1+1,indvars(3))

d2_rhs_u_dxdy_6_0_nxp4p0m21m4p1k = -&
          0.5_wp*d2_rhs_u_dxdy_6_0_nxp4p0m21m4p1k_nxp4p0m21m4p1m1k+&
          0.5_wp*d2_rhs_u_dxdy_6_0_nxp4p0m21m4p1k_nxp4p0m21m4p1p1k

d2_rhs_u_dxdy_6_0_nxp4p0m21m4p1k = d2_rhs_u_dxdy_6_0_nxp4p0m21m4p1k*param_float(2)

d1_rhs_u_dx_6_nxp4p0p01m4p1k = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0+0,1-4+1,indvars(5))/1.0_wp*q(nx+4+0+0,1-4+1,indvars(1)))**3.0_wp/((q(nx+4+0+0,1-4+1,indvars(5))/1.0_wp*q(nx+4+0+0,1-4+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0+0,1-4+1,indvars(5))/1.0_wp*q(nx+4+0+0,1-4+1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0+0,1-4+1,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0p01m4p1k)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0+0,1-4+1,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0p01m4p1k)+&
                    qst(nx+4+0+0,1-4+1,indvarsst(11))*(d2_rhs_u_dxdy_6_0_nxp4p0p01m4p1k)))

d1_rhs_u_dx_6_nxp4p0m11m4p1k = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0-1,1-4+1,indvars(5))/1.0_wp*q(nx+4+0-1,1-4+1,indvars(1)))**3.0_wp/((q(nx+4+0-1,1-4+1,indvars(5))/1.0_wp*q(nx+4+0-1,1-4+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0-1,1-4+1,indvars(5))/1.0_wp*q(nx+4+0-1,1-4+1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0-1,1-4+1,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0m11m4p1k)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0-1,1-4+1,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0m11m4p1k)+&
                    qst(nx+4+0-1,1-4+1,indvarsst(11))*(d2_rhs_u_dxdy_6_0_nxp4p0m11m4p1k)))

d1_rhs_u_dx_6_nxp4p0m21m4p1k = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0-2,1-4+1,indvars(5))/1.0_wp*q(nx+4+0-2,1-4+1,indvars(1)))**3.0_wp/((q(nx+4+0-2,1-4+1,indvars(5))/1.0_wp*q(nx+4+0-2,1-4+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0-2,1-4+1,indvars(5))/1.0_wp*q(nx+4+0-2,1-4+1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0-2,1-4+1,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0m21m4p1k)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0-2,1-4+1,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0m21m4p1k)+&
                    qst(nx+4+0-2,1-4+1,indvarsst(11))*(d2_rhs_u_dxdy_6_0_nxp4p0m21m4p1k)))

d1_rhs_u_dx_6_nxp4p01m4p1k = 1.5_wp*d1_rhs_u_dx_6_nxp4p0p01m4p1k-&
          2.0_wp*d1_rhs_u_dx_6_nxp4p0m11m4p1k+&
          0.5_wp*d1_rhs_u_dx_6_nxp4p0m21m4p1k

d1_rhs_u_dx_6_nxp4p01m4p1k = d1_rhs_u_dx_6_nxp4p01m4p1k*param_float(1)

d1_rhs_u_dy_0_nxp4p01m4p1m1k = q(nx+4+0,1-4+1-1,indvars(1))*q(nx+4+0,1-4+1-1,indvars(2))*q(nx+4+0,1-4+1-1,indvars(3))

d1_rhs_u_dy_0_nxp4p01m4p1p1k = q(nx+4+0,1-4+1+1,indvars(1))*q(nx+4+0,1-4+1+1,indvars(2))*q(nx+4+0,1-4+1+1,indvars(3))

d1_rhs_u_dy_0_nxp4p01m4p1k = -&
          0.5_wp*d1_rhs_u_dy_0_nxp4p01m4p1m1k+&
          0.5_wp*d1_rhs_u_dy_0_nxp4p01m4p1p1k

d1_rhs_u_dy_0_nxp4p01m4p1k = d1_rhs_u_dy_0_nxp4p01m4p1k*param_float(2)

d2_rhs_u_dydx_1_0_nxp4p01m4p1m1k_nxp4p0p01m4p1m1k = q(nx+4+0+0,1-4+1-1,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p01m4p1m1k_nxp4p0m11m4p1m1k = q(nx+4+0-1,1-4+1-1,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p01m4p1m1k_nxp4p0m21m4p1m1k = q(nx+4+0-2,1-4+1-1,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p01m4p1m1k = 1.5_wp*d2_rhs_u_dydx_1_0_nxp4p01m4p1m1k_nxp4p0p01m4p1m1k-&
          2.0_wp*d2_rhs_u_dydx_1_0_nxp4p01m4p1m1k_nxp4p0m11m4p1m1k+&
          0.5_wp*d2_rhs_u_dydx_1_0_nxp4p01m4p1m1k_nxp4p0m21m4p1m1k

d2_rhs_u_dydx_1_0_nxp4p01m4p1m1k = d2_rhs_u_dydx_1_0_nxp4p01m4p1m1k*param_float(1)

d2_rhs_u_dydx_1_0_nxp4p01m4p1p1k_nxp4p0p01m4p1p1k = q(nx+4+0+0,1-4+1+1,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p01m4p1p1k_nxp4p0m11m4p1p1k = q(nx+4+0-1,1-4+1+1,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p01m4p1p1k_nxp4p0m21m4p1p1k = q(nx+4+0-2,1-4+1+1,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p01m4p1p1k = 1.5_wp*d2_rhs_u_dydx_1_0_nxp4p01m4p1p1k_nxp4p0p01m4p1p1k-&
          2.0_wp*d2_rhs_u_dydx_1_0_nxp4p01m4p1p1k_nxp4p0m11m4p1p1k+&
          0.5_wp*d2_rhs_u_dydx_1_0_nxp4p01m4p1p1k_nxp4p0m21m4p1p1k

d2_rhs_u_dydx_1_0_nxp4p01m4p1p1k = d2_rhs_u_dydx_1_0_nxp4p01m4p1p1k*param_float(1)

d2_rhs_u_dydy_1_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p0k = q(nx+4+0,1-4+1-1+0,indvars(2))

d2_rhs_u_dydy_1_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p1k = q(nx+4+0,1-4+1-1+1,indvars(2))

d2_rhs_u_dydy_1_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p2k = q(nx+4+0,1-4+1-1+2,indvars(2))

d2_rhs_u_dydy_1_0_nxp4p01m4p1m1k = -&
          1.5_wp*d2_rhs_u_dydy_1_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p0k+&
          2.0_wp*d2_rhs_u_dydy_1_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p1k-&
          0.5_wp*d2_rhs_u_dydy_1_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p2k

d2_rhs_u_dydy_1_0_nxp4p01m4p1m1k = d2_rhs_u_dydy_1_0_nxp4p01m4p1m1k*param_float(2)

d2_rhs_u_dydy_1_0_nxp4p01m4p1p1k_nxp4p01m4p1p1m1k = q(nx+4+0,1-4+1+1-1,indvars(2))

d2_rhs_u_dydy_1_0_nxp4p01m4p1p1k_nxp4p01m4p1p1p1k = q(nx+4+0,1-4+1+1+1,indvars(2))

d2_rhs_u_dydy_1_0_nxp4p01m4p1p1k = -&
          0.5_wp*d2_rhs_u_dydy_1_0_nxp4p01m4p1p1k_nxp4p01m4p1p1m1k+&
          0.5_wp*d2_rhs_u_dydy_1_0_nxp4p01m4p1p1k_nxp4p01m4p1p1p1k

d2_rhs_u_dydy_1_0_nxp4p01m4p1p1k = d2_rhs_u_dydy_1_0_nxp4p01m4p1p1k*param_float(2)

d1_rhs_u_dy_1_nxp4p01m4p1m1k = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0,1-4+1-1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1-1,indvars(1)))**3.0_wp/((q(nx+4+0,1-4+1-1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1-1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0,1-4+1-1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1-1,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0,1-4+1-1,indvarsst(11))*(d2_rhs_u_dydy_1_0_nxp4p01m4p1m1k)+&
                    qst(nx+4+0,1-4+1-1,indvarsst(10))*(d2_rhs_u_dydx_1_0_nxp4p01m4p1m1k))

d1_rhs_u_dy_1_nxp4p01m4p1p1k = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0,1-4+1+1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1+1,indvars(1)))**3.0_wp/((q(nx+4+0,1-4+1+1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0,1-4+1+1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1+1,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0,1-4+1+1,indvarsst(11))*(d2_rhs_u_dydy_1_0_nxp4p01m4p1p1k)+&
                    qst(nx+4+0,1-4+1+1,indvarsst(10))*(d2_rhs_u_dydx_1_0_nxp4p01m4p1p1k))

d1_rhs_u_dy_1_nxp4p01m4p1k = -&
          0.5_wp*d1_rhs_u_dy_1_nxp4p01m4p1m1k+&
          0.5_wp*d1_rhs_u_dy_1_nxp4p01m4p1p1k

d1_rhs_u_dy_1_nxp4p01m4p1k = d1_rhs_u_dy_1_nxp4p01m4p1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None d(rho u)/dt ***********
!                                                           
!***********************************************************


rhs(nx+4+0,1-4+1,indvars(2)) =   -  ( q(nx+4+0,1-4+1,indvars(2))*(1.0_wp/(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5**2*((q(nx+4+0,1-4+1,indvars(2))*(-&
                    (param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5**2*d1_rhs_rho_dx_1_nxp4p01m4p1k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p01m4p1k))*qst(nx+4+0,1-4+1,indvarsst(10))+&
                    0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5+&
                    q(nx+4+0,1-4+1,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5*q(nx+4+0,1-4+1,indvars(1))*d1_rhs_rho_dx_3_nxp4p01m4p1k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p01m4p1k))*qst(nx+4+0,1-4+1,indvarsst(10))+&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5*(1-&
                    qface_i(1-4+1,2)*qface_i(1-4+1,2))/param_float(20 + 5)*((param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))-&
                    param_float(26 + 5)))))+&
                    q(nx+4+0,1-4+1,indvars(1))*(1.0_wp/(2*q(nx+4+0,1-4+1,indvars(1))*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5)*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5+&
                    q(nx+4+0,1-4+1,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5*q(nx+4+0,1-4+1,indvars(1))*d1_rhs_rho_dx_3_nxp4p01m4p1k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p01m4p1k))*qst(nx+4+0,1-4+1,indvarsst(10))-&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5*(1-&
                    qface_i(1-4+1,2)*qface_i(1-4+1,2))/param_float(20 + 5)*((param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))-&
                    param_float(26 + 5))))+&
                    d1_rhs_u_dy_0_nxp4p01m4p1k+&
                    qst(nx+4+0,1-4+1,indvarsst(10))*(d1_rhs_u_dx_6_nxp4p01m4p1k)+&
                    qst(nx+4+0,1-4+1,indvarsst(11))*(d1_rhs_u_dy_1_nxp4p01m4p1k) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None d(rho v)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! v*(1.0_wp/c**2*((u*(-c**2*[rho]_1x+1.0_wp*[p]_1x))*deltaxI+0.5_wp*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI+esse*c*(1-M_imax*M_imax)/L_ref*(p-P0))))+rho*((1.0_wp*u*[v]_1x)*deltaxI)+[rho*v*v]_1y+deltaxI*([-visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x))]_1x)+deltayI*([-visc_t*(2.0_wp*deltayI*({v}_1y)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y)))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_rhs_v_dx_4_nxp4p0p01m4p1k = q(nx+4+0+0,1-4+1,indvars(3))

d1_rhs_v_dx_4_nxp4p0m11m4p1k = q(nx+4+0-1,1-4+1,indvars(3))

d1_rhs_v_dx_4_nxp4p0m21m4p1k = q(nx+4+0-2,1-4+1,indvars(3))

d1_rhs_v_dx_4_nxp4p01m4p1k = 1.5_wp*d1_rhs_v_dx_4_nxp4p0p01m4p1k-&
          2.0_wp*d1_rhs_v_dx_4_nxp4p0m11m4p1k+&
          0.5_wp*d1_rhs_v_dx_4_nxp4p0m21m4p1k

d1_rhs_v_dx_4_nxp4p01m4p1k = d1_rhs_v_dx_4_nxp4p01m4p1k*param_float(1)

d2_rhs_v_dxdx_5_0_nxp4p0p01m4p1k_nxp4p0p0p01m4p1k = q(nx+4+0+0+0,1-4+1,indvars(3))

d2_rhs_v_dxdx_5_0_nxp4p0p01m4p1k_nxp4p0p0m11m4p1k = q(nx+4+0+0-1,1-4+1,indvars(3))

d2_rhs_v_dxdx_5_0_nxp4p0p01m4p1k_nxp4p0p0m21m4p1k = q(nx+4+0+0-2,1-4+1,indvars(3))

d2_rhs_v_dxdx_5_0_nxp4p0p01m4p1k = 1.5_wp*d2_rhs_v_dxdx_5_0_nxp4p0p01m4p1k_nxp4p0p0p01m4p1k-&
          2.0_wp*d2_rhs_v_dxdx_5_0_nxp4p0p01m4p1k_nxp4p0p0m11m4p1k+&
          0.5_wp*d2_rhs_v_dxdx_5_0_nxp4p0p01m4p1k_nxp4p0p0m21m4p1k

d2_rhs_v_dxdx_5_0_nxp4p0p01m4p1k = d2_rhs_v_dxdx_5_0_nxp4p0p01m4p1k*param_float(1)

d2_rhs_v_dxdx_5_0_nxp4p0m11m4p1k_nxp4p0m1m11m4p1k = q(nx+4+0-1-1,1-4+1,indvars(3))

d2_rhs_v_dxdx_5_0_nxp4p0m11m4p1k_nxp4p0m1p11m4p1k = q(nx+4+0-1+1,1-4+1,indvars(3))

d2_rhs_v_dxdx_5_0_nxp4p0m11m4p1k = -&
          0.5_wp*d2_rhs_v_dxdx_5_0_nxp4p0m11m4p1k_nxp4p0m1m11m4p1k+&
          0.5_wp*d2_rhs_v_dxdx_5_0_nxp4p0m11m4p1k_nxp4p0m1p11m4p1k

d2_rhs_v_dxdx_5_0_nxp4p0m11m4p1k = d2_rhs_v_dxdx_5_0_nxp4p0m11m4p1k*param_float(1)

d2_rhs_v_dxdx_5_0_nxp4p0m21m4p1k_nxp4p0m2m11m4p1k = q(nx+4+0-2-1,1-4+1,indvars(3))

d2_rhs_v_dxdx_5_0_nxp4p0m21m4p1k_nxp4p0m2p11m4p1k = q(nx+4+0-2+1,1-4+1,indvars(3))

d2_rhs_v_dxdx_5_0_nxp4p0m21m4p1k = -&
          0.5_wp*d2_rhs_v_dxdx_5_0_nxp4p0m21m4p1k_nxp4p0m2m11m4p1k+&
          0.5_wp*d2_rhs_v_dxdx_5_0_nxp4p0m21m4p1k_nxp4p0m2p11m4p1k

d2_rhs_v_dxdx_5_0_nxp4p0m21m4p1k = d2_rhs_v_dxdx_5_0_nxp4p0m21m4p1k*param_float(1)

d2_rhs_v_dxdy_5_0_nxp4p0p01m4p1k_nxp4p0p01m4p1m1k = q(nx+4+0+0,1-4+1-1,indvars(2))

d2_rhs_v_dxdy_5_0_nxp4p0p01m4p1k_nxp4p0p01m4p1p1k = q(nx+4+0+0,1-4+1+1,indvars(2))

d2_rhs_v_dxdy_5_0_nxp4p0p01m4p1k = -&
          0.5_wp*d2_rhs_v_dxdy_5_0_nxp4p0p01m4p1k_nxp4p0p01m4p1m1k+&
          0.5_wp*d2_rhs_v_dxdy_5_0_nxp4p0p01m4p1k_nxp4p0p01m4p1p1k

d2_rhs_v_dxdy_5_0_nxp4p0p01m4p1k = d2_rhs_v_dxdy_5_0_nxp4p0p01m4p1k*param_float(2)

d2_rhs_v_dxdy_5_0_nxp4p0m11m4p1k_nxp4p0m11m4p1m1k = q(nx+4+0-1,1-4+1-1,indvars(2))

d2_rhs_v_dxdy_5_0_nxp4p0m11m4p1k_nxp4p0m11m4p1p1k = q(nx+4+0-1,1-4+1+1,indvars(2))

d2_rhs_v_dxdy_5_0_nxp4p0m11m4p1k = -&
          0.5_wp*d2_rhs_v_dxdy_5_0_nxp4p0m11m4p1k_nxp4p0m11m4p1m1k+&
          0.5_wp*d2_rhs_v_dxdy_5_0_nxp4p0m11m4p1k_nxp4p0m11m4p1p1k

d2_rhs_v_dxdy_5_0_nxp4p0m11m4p1k = d2_rhs_v_dxdy_5_0_nxp4p0m11m4p1k*param_float(2)

d2_rhs_v_dxdy_5_0_nxp4p0m21m4p1k_nxp4p0m21m4p1m1k = q(nx+4+0-2,1-4+1-1,indvars(2))

d2_rhs_v_dxdy_5_0_nxp4p0m21m4p1k_nxp4p0m21m4p1p1k = q(nx+4+0-2,1-4+1+1,indvars(2))

d2_rhs_v_dxdy_5_0_nxp4p0m21m4p1k = -&
          0.5_wp*d2_rhs_v_dxdy_5_0_nxp4p0m21m4p1k_nxp4p0m21m4p1m1k+&
          0.5_wp*d2_rhs_v_dxdy_5_0_nxp4p0m21m4p1k_nxp4p0m21m4p1p1k

d2_rhs_v_dxdy_5_0_nxp4p0m21m4p1k = d2_rhs_v_dxdy_5_0_nxp4p0m21m4p1k*param_float(2)

d1_rhs_v_dx_5_nxp4p0p01m4p1k = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0+0,1-4+1,indvars(5))/1.0_wp*q(nx+4+0+0,1-4+1,indvars(1)))**3.0_wp/((q(nx+4+0+0,1-4+1,indvars(5))/1.0_wp*q(nx+4+0+0,1-4+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0+0,1-4+1,indvars(5))/1.0_wp*q(nx+4+0+0,1-4+1,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0+0,1-4+1,indvarsst(11))*(d2_rhs_v_dxdy_5_0_nxp4p0p01m4p1k)+&
                    qst(nx+4+0+0,1-4+1,indvarsst(10))*(d2_rhs_v_dxdx_5_0_nxp4p0p01m4p1k))

d1_rhs_v_dx_5_nxp4p0m11m4p1k = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0-1,1-4+1,indvars(5))/1.0_wp*q(nx+4+0-1,1-4+1,indvars(1)))**3.0_wp/((q(nx+4+0-1,1-4+1,indvars(5))/1.0_wp*q(nx+4+0-1,1-4+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0-1,1-4+1,indvars(5))/1.0_wp*q(nx+4+0-1,1-4+1,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0-1,1-4+1,indvarsst(11))*(d2_rhs_v_dxdy_5_0_nxp4p0m11m4p1k)+&
                    qst(nx+4+0-1,1-4+1,indvarsst(10))*(d2_rhs_v_dxdx_5_0_nxp4p0m11m4p1k))

d1_rhs_v_dx_5_nxp4p0m21m4p1k = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0-2,1-4+1,indvars(5))/1.0_wp*q(nx+4+0-2,1-4+1,indvars(1)))**3.0_wp/((q(nx+4+0-2,1-4+1,indvars(5))/1.0_wp*q(nx+4+0-2,1-4+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0-2,1-4+1,indvars(5))/1.0_wp*q(nx+4+0-2,1-4+1,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0-2,1-4+1,indvarsst(11))*(d2_rhs_v_dxdy_5_0_nxp4p0m21m4p1k)+&
                    qst(nx+4+0-2,1-4+1,indvarsst(10))*(d2_rhs_v_dxdx_5_0_nxp4p0m21m4p1k))

d1_rhs_v_dx_5_nxp4p01m4p1k = 1.5_wp*d1_rhs_v_dx_5_nxp4p0p01m4p1k-&
          2.0_wp*d1_rhs_v_dx_5_nxp4p0m11m4p1k+&
          0.5_wp*d1_rhs_v_dx_5_nxp4p0m21m4p1k

d1_rhs_v_dx_5_nxp4p01m4p1k = d1_rhs_v_dx_5_nxp4p01m4p1k*param_float(1)

d1_rhs_v_dy_0_nxp4p01m4p1m1k = q(nx+4+0,1-4+1-1,indvars(1))*q(nx+4+0,1-4+1-1,indvars(3))*q(nx+4+0,1-4+1-1,indvars(3))

d1_rhs_v_dy_0_nxp4p01m4p1p1k = q(nx+4+0,1-4+1+1,indvars(1))*q(nx+4+0,1-4+1+1,indvars(3))*q(nx+4+0,1-4+1+1,indvars(3))

d1_rhs_v_dy_0_nxp4p01m4p1k = -&
          0.5_wp*d1_rhs_v_dy_0_nxp4p01m4p1m1k+&
          0.5_wp*d1_rhs_v_dy_0_nxp4p01m4p1p1k

d1_rhs_v_dy_0_nxp4p01m4p1k = d1_rhs_v_dy_0_nxp4p01m4p1k*param_float(2)

d2_rhs_v_dydx_1_0_nxp4p01m4p1m1k_nxp4p0p01m4p1m1k = q(nx+4+0+0,1-4+1-1,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p01m4p1m1k_nxp4p0m11m4p1m1k = q(nx+4+0-1,1-4+1-1,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p01m4p1m1k_nxp4p0m21m4p1m1k = q(nx+4+0-2,1-4+1-1,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p01m4p1m1k = 1.5_wp*d2_rhs_v_dydx_1_0_nxp4p01m4p1m1k_nxp4p0p01m4p1m1k-&
          2.0_wp*d2_rhs_v_dydx_1_0_nxp4p01m4p1m1k_nxp4p0m11m4p1m1k+&
          0.5_wp*d2_rhs_v_dydx_1_0_nxp4p01m4p1m1k_nxp4p0m21m4p1m1k

d2_rhs_v_dydx_1_0_nxp4p01m4p1m1k = d2_rhs_v_dydx_1_0_nxp4p01m4p1m1k*param_float(1)

d2_rhs_v_dydx_1_0_nxp4p01m4p1p1k_nxp4p0p01m4p1p1k = q(nx+4+0+0,1-4+1+1,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p01m4p1p1k_nxp4p0m11m4p1p1k = q(nx+4+0-1,1-4+1+1,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p01m4p1p1k_nxp4p0m21m4p1p1k = q(nx+4+0-2,1-4+1+1,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p01m4p1p1k = 1.5_wp*d2_rhs_v_dydx_1_0_nxp4p01m4p1p1k_nxp4p0p01m4p1p1k-&
          2.0_wp*d2_rhs_v_dydx_1_0_nxp4p01m4p1p1k_nxp4p0m11m4p1p1k+&
          0.5_wp*d2_rhs_v_dydx_1_0_nxp4p01m4p1p1k_nxp4p0m21m4p1p1k

d2_rhs_v_dydx_1_0_nxp4p01m4p1p1k = d2_rhs_v_dydx_1_0_nxp4p01m4p1p1k*param_float(1)

d2_rhs_v_dydy_1_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p0k = q(nx+4+0,1-4+1-1+0,indvars(3))

d2_rhs_v_dydy_1_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p1k = q(nx+4+0,1-4+1-1+1,indvars(3))

d2_rhs_v_dydy_1_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p2k = q(nx+4+0,1-4+1-1+2,indvars(3))

d2_rhs_v_dydy_1_0_nxp4p01m4p1m1k = -&
          1.5_wp*d2_rhs_v_dydy_1_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p0k+&
          2.0_wp*d2_rhs_v_dydy_1_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p1k-&
          0.5_wp*d2_rhs_v_dydy_1_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p2k

d2_rhs_v_dydy_1_0_nxp4p01m4p1m1k = d2_rhs_v_dydy_1_0_nxp4p01m4p1m1k*param_float(2)

d2_rhs_v_dydy_1_0_nxp4p01m4p1p1k_nxp4p01m4p1p1m1k = q(nx+4+0,1-4+1+1-1,indvars(3))

d2_rhs_v_dydy_1_0_nxp4p01m4p1p1k_nxp4p01m4p1p1p1k = q(nx+4+0,1-4+1+1+1,indvars(3))

d2_rhs_v_dydy_1_0_nxp4p01m4p1p1k = -&
          0.5_wp*d2_rhs_v_dydy_1_0_nxp4p01m4p1p1k_nxp4p01m4p1p1m1k+&
          0.5_wp*d2_rhs_v_dydy_1_0_nxp4p01m4p1p1k_nxp4p01m4p1p1p1k

d2_rhs_v_dydy_1_0_nxp4p01m4p1p1k = d2_rhs_v_dydy_1_0_nxp4p01m4p1p1k*param_float(2)

d1_rhs_v_dy_1_nxp4p01m4p1m1k = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0,1-4+1-1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1-1,indvars(1)))**3.0_wp/((q(nx+4+0,1-4+1-1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1-1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0,1-4+1-1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1-1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0,1-4+1-1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp4p01m4p1m1k)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0,1-4+1-1,indvarsst(10))*(d2_rhs_v_dydx_1_0_nxp4p01m4p1m1k)+&
                    qst(nx+4+0,1-4+1-1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp4p01m4p1m1k)))

d1_rhs_v_dy_1_nxp4p01m4p1p1k = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0,1-4+1+1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1+1,indvars(1)))**3.0_wp/((q(nx+4+0,1-4+1+1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0,1-4+1+1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1+1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0,1-4+1+1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp4p01m4p1p1k)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0,1-4+1+1,indvarsst(10))*(d2_rhs_v_dydx_1_0_nxp4p01m4p1p1k)+&
                    qst(nx+4+0,1-4+1+1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp4p01m4p1p1k)))

d1_rhs_v_dy_1_nxp4p01m4p1k = -&
          0.5_wp*d1_rhs_v_dy_1_nxp4p01m4p1m1k+&
          0.5_wp*d1_rhs_v_dy_1_nxp4p01m4p1p1k

d1_rhs_v_dy_1_nxp4p01m4p1k = d1_rhs_v_dy_1_nxp4p01m4p1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None d(rho v)/dt ***********
!                                                           
!***********************************************************


rhs(nx+4+0,1-4+1,indvars(3)) =   -  ( q(nx+4+0,1-4+1,indvars(3))*(1.0_wp/(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5**2*((q(nx+4+0,1-4+1,indvars(2))*(-&
                    (param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5**2*d1_rhs_rho_dx_1_nxp4p01m4p1k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p01m4p1k))*qst(nx+4+0,1-4+1,indvarsst(10))+&
                    0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5+&
                    q(nx+4+0,1-4+1,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5*q(nx+4+0,1-4+1,indvars(1))*d1_rhs_rho_dx_3_nxp4p01m4p1k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p01m4p1k))*qst(nx+4+0,1-4+1,indvarsst(10))+&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5*(1-&
                    qface_i(1-4+1,2)*qface_i(1-4+1,2))/param_float(20 + 5)*((param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))-&
                    param_float(26 + 5)))))+&
                    q(nx+4+0,1-4+1,indvars(1))*((1.0_wp*q(nx+4+0,1-4+1,indvars(2))*d1_rhs_v_dx_4_nxp4p01m4p1k)*qst(nx+4+0,1-4+1,indvarsst(10)))+&
                    d1_rhs_v_dy_0_nxp4p01m4p1k+&
                    qst(nx+4+0,1-4+1,indvarsst(10))*(d1_rhs_v_dx_5_nxp4p01m4p1k)+&
                    qst(nx+4+0,1-4+1,indvarsst(11))*(d1_rhs_v_dy_1_nxp4p01m4p1k) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None d(rho et)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! 0.5_wp*(u**2+v**2)*(1.0_wp/c**2*((u*(-c**2*[rho]_1x+1.0_wp*[p]_1x))*deltaxI+0.5_wp*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI+esse*c*(1-M_imax*M_imax)/L_ref*(p-P0))))+1.0_wp/gamma_m1*(0.5_wp*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI+esse*c*(1-M_imax*M_imax)/L_ref*(p-P0)))+rho*u*(1.0_wp/(2*rho*c)*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI-esse*c*(1-M_imax*M_imax)/L_ref*(p-P0)))+rho*v*((1.0_wp*u*[v]_1x)*deltaxI)+[(rho*et+p)*v]_1y+deltaxI*([-kappa*deltaxI*({T}_1x)-u*(visc_t*(2.0_wp*deltaxI*({u}_1x)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y))))-v*(visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x)))]_1x)+deltayI*([-kappa*deltayI*({T}_1y)-u*(visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x)))-v*(visc_t*(2.0_wp*deltayI*({v}_1y)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y))))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_rhs_et_dxdx_9_0_nxp4p0p01m4p1k_nxp4p0p0p01m4p1k = ((q(nx+4+0+0+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0+0,1-4+1,indvars(2))*q(nx+4+0+0+0,1-4+1,indvars(2))+&
                    q(nx+4+0+0+0,1-4+1,indvars(3))*q(nx+4+0+0+0,1-4+1,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dxdx_9_0_nxp4p0p01m4p1k_nxp4p0p0m11m4p1k = ((q(nx+4+0+0-1,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0-1,1-4+1,indvars(2))*q(nx+4+0+0-1,1-4+1,indvars(2))+&
                    q(nx+4+0+0-1,1-4+1,indvars(3))*q(nx+4+0+0-1,1-4+1,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dxdx_9_0_nxp4p0p01m4p1k_nxp4p0p0m21m4p1k = ((q(nx+4+0+0-2,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0-2,1-4+1,indvars(2))*q(nx+4+0+0-2,1-4+1,indvars(2))+&
                    q(nx+4+0+0-2,1-4+1,indvars(3))*q(nx+4+0+0-2,1-4+1,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dxdx_9_0_nxp4p0p01m4p1k = 1.5_wp*d2_rhs_et_dxdx_9_0_nxp4p0p01m4p1k_nxp4p0p0p01m4p1k-&
          2.0_wp*d2_rhs_et_dxdx_9_0_nxp4p0p01m4p1k_nxp4p0p0m11m4p1k+&
          0.5_wp*d2_rhs_et_dxdx_9_0_nxp4p0p01m4p1k_nxp4p0p0m21m4p1k

d2_rhs_et_dxdx_9_0_nxp4p0p01m4p1k = d2_rhs_et_dxdx_9_0_nxp4p0p01m4p1k*param_float(1)

d2_rhs_et_dxdx_9_0_nxp4p0m11m4p1k_nxp4p0m1m11m4p1k = ((q(nx+4+0-1-1,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1-1,1-4+1,indvars(2))*q(nx+4+0-1-1,1-4+1,indvars(2))+&
                    q(nx+4+0-1-1,1-4+1,indvars(3))*q(nx+4+0-1-1,1-4+1,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dxdx_9_0_nxp4p0m11m4p1k_nxp4p0m1p11m4p1k = ((q(nx+4+0-1+1,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1+1,1-4+1,indvars(2))*q(nx+4+0-1+1,1-4+1,indvars(2))+&
                    q(nx+4+0-1+1,1-4+1,indvars(3))*q(nx+4+0-1+1,1-4+1,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dxdx_9_0_nxp4p0m11m4p1k = -&
          0.5_wp*d2_rhs_et_dxdx_9_0_nxp4p0m11m4p1k_nxp4p0m1m11m4p1k+&
          0.5_wp*d2_rhs_et_dxdx_9_0_nxp4p0m11m4p1k_nxp4p0m1p11m4p1k

d2_rhs_et_dxdx_9_0_nxp4p0m11m4p1k = d2_rhs_et_dxdx_9_0_nxp4p0m11m4p1k*param_float(1)

d2_rhs_et_dxdx_9_0_nxp4p0m21m4p1k_nxp4p0m2m11m4p1k = ((q(nx+4+0-2-1,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2-1,1-4+1,indvars(2))*q(nx+4+0-2-1,1-4+1,indvars(2))+&
                    q(nx+4+0-2-1,1-4+1,indvars(3))*q(nx+4+0-2-1,1-4+1,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dxdx_9_0_nxp4p0m21m4p1k_nxp4p0m2p11m4p1k = ((q(nx+4+0-2+1,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2+1,1-4+1,indvars(2))*q(nx+4+0-2+1,1-4+1,indvars(2))+&
                    q(nx+4+0-2+1,1-4+1,indvars(3))*q(nx+4+0-2+1,1-4+1,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dxdx_9_0_nxp4p0m21m4p1k = -&
          0.5_wp*d2_rhs_et_dxdx_9_0_nxp4p0m21m4p1k_nxp4p0m2m11m4p1k+&
          0.5_wp*d2_rhs_et_dxdx_9_0_nxp4p0m21m4p1k_nxp4p0m2p11m4p1k

d2_rhs_et_dxdx_9_0_nxp4p0m21m4p1k = d2_rhs_et_dxdx_9_0_nxp4p0m21m4p1k*param_float(1)

d1_rhs_et_dx_9_nxp4p0p01m4p1k = -param_float(2 + 5)*qst(nx+4+0+0,1-4+1,indvarsst(10))*(d2_rhs_et_dxdx_9_0_nxp4p0p01m4p1k)-&
                    q(nx+4+0+0,1-4+1,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0+0,1-4+1,indvars(5))/1.0_wp*q(nx+4+0+0,1-4+1,indvars(1)))**3.0_wp/((q(nx+4+0+0,1-4+1,indvars(5))/1.0_wp*q(nx+4+0+0,1-4+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0+0,1-4+1,indvars(5))/1.0_wp*q(nx+4+0+0,1-4+1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0+0,1-4+1,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0p01m4p1k)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0+0,1-4+1,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0p01m4p1k)+&
                    qst(nx+4+0+0,1-4+1,indvarsst(11))*(d2_rhs_u_dxdy_6_0_nxp4p0p01m4p1k))))-&
                    q(nx+4+0+0,1-4+1,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0+0,1-4+1,indvars(5))/1.0_wp*q(nx+4+0+0,1-4+1,indvars(1)))**3.0_wp/((q(nx+4+0+0,1-4+1,indvars(5))/1.0_wp*q(nx+4+0+0,1-4+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0+0,1-4+1,indvars(5))/1.0_wp*q(nx+4+0+0,1-4+1,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0+0,1-4+1,indvarsst(11))*(d2_rhs_v_dxdy_5_0_nxp4p0p01m4p1k)+&
                    qst(nx+4+0+0,1-4+1,indvarsst(10))*(d2_rhs_v_dxdx_5_0_nxp4p0p01m4p1k)))

d1_rhs_et_dx_9_nxp4p0m11m4p1k = -param_float(2 + 5)*qst(nx+4+0-1,1-4+1,indvarsst(10))*(d2_rhs_et_dxdx_9_0_nxp4p0m11m4p1k)-&
                    q(nx+4+0-1,1-4+1,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0-1,1-4+1,indvars(5))/1.0_wp*q(nx+4+0-1,1-4+1,indvars(1)))**3.0_wp/((q(nx+4+0-1,1-4+1,indvars(5))/1.0_wp*q(nx+4+0-1,1-4+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0-1,1-4+1,indvars(5))/1.0_wp*q(nx+4+0-1,1-4+1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0-1,1-4+1,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0m11m4p1k)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0-1,1-4+1,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0m11m4p1k)+&
                    qst(nx+4+0-1,1-4+1,indvarsst(11))*(d2_rhs_u_dxdy_6_0_nxp4p0m11m4p1k))))-&
                    q(nx+4+0-1,1-4+1,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0-1,1-4+1,indvars(5))/1.0_wp*q(nx+4+0-1,1-4+1,indvars(1)))**3.0_wp/((q(nx+4+0-1,1-4+1,indvars(5))/1.0_wp*q(nx+4+0-1,1-4+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0-1,1-4+1,indvars(5))/1.0_wp*q(nx+4+0-1,1-4+1,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0-1,1-4+1,indvarsst(11))*(d2_rhs_v_dxdy_5_0_nxp4p0m11m4p1k)+&
                    qst(nx+4+0-1,1-4+1,indvarsst(10))*(d2_rhs_v_dxdx_5_0_nxp4p0m11m4p1k)))

d1_rhs_et_dx_9_nxp4p0m21m4p1k = -param_float(2 + 5)*qst(nx+4+0-2,1-4+1,indvarsst(10))*(d2_rhs_et_dxdx_9_0_nxp4p0m21m4p1k)-&
                    q(nx+4+0-2,1-4+1,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0-2,1-4+1,indvars(5))/1.0_wp*q(nx+4+0-2,1-4+1,indvars(1)))**3.0_wp/((q(nx+4+0-2,1-4+1,indvars(5))/1.0_wp*q(nx+4+0-2,1-4+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0-2,1-4+1,indvars(5))/1.0_wp*q(nx+4+0-2,1-4+1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0-2,1-4+1,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0m21m4p1k)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0-2,1-4+1,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0m21m4p1k)+&
                    qst(nx+4+0-2,1-4+1,indvarsst(11))*(d2_rhs_u_dxdy_6_0_nxp4p0m21m4p1k))))-&
                    q(nx+4+0-2,1-4+1,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0-2,1-4+1,indvars(5))/1.0_wp*q(nx+4+0-2,1-4+1,indvars(1)))**3.0_wp/((q(nx+4+0-2,1-4+1,indvars(5))/1.0_wp*q(nx+4+0-2,1-4+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0-2,1-4+1,indvars(5))/1.0_wp*q(nx+4+0-2,1-4+1,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0-2,1-4+1,indvarsst(11))*(d2_rhs_v_dxdy_5_0_nxp4p0m21m4p1k)+&
                    qst(nx+4+0-2,1-4+1,indvarsst(10))*(d2_rhs_v_dxdx_5_0_nxp4p0m21m4p1k)))

d1_rhs_et_dx_9_nxp4p01m4p1k = 1.5_wp*d1_rhs_et_dx_9_nxp4p0p01m4p1k-&
          2.0_wp*d1_rhs_et_dx_9_nxp4p0m11m4p1k+&
          0.5_wp*d1_rhs_et_dx_9_nxp4p0m21m4p1k

d1_rhs_et_dx_9_nxp4p01m4p1k = d1_rhs_et_dx_9_nxp4p01m4p1k*param_float(1)

d1_rhs_et_dy_0_nxp4p01m4p1m1k = (q(nx+4+0,1-4+1-1,indvars(1))*q(nx+4+0,1-4+1-1,indvars(4))+&
                    (param_float(3 + 5))*q(nx+4+0,1-4+1-1,indvars(1))*((q(nx+4+0,1-4+1-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1-1,indvars(2))*q(nx+4+0,1-4+1-1,indvars(2))+&
                    q(nx+4+0,1-4+1-1,indvars(3))*q(nx+4+0,1-4+1-1,indvars(3))))))*q(nx+4+0,1-4+1-1,indvars(3))

d1_rhs_et_dy_0_nxp4p01m4p1p1k = (q(nx+4+0,1-4+1+1,indvars(1))*q(nx+4+0,1-4+1+1,indvars(4))+&
                    (param_float(3 + 5))*q(nx+4+0,1-4+1+1,indvars(1))*((q(nx+4+0,1-4+1+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1+1,indvars(2))*q(nx+4+0,1-4+1+1,indvars(2))+&
                    q(nx+4+0,1-4+1+1,indvars(3))*q(nx+4+0,1-4+1+1,indvars(3))))))*q(nx+4+0,1-4+1+1,indvars(3))

d1_rhs_et_dy_0_nxp4p01m4p1k = -&
          0.5_wp*d1_rhs_et_dy_0_nxp4p01m4p1m1k+&
          0.5_wp*d1_rhs_et_dy_0_nxp4p01m4p1p1k

d1_rhs_et_dy_0_nxp4p01m4p1k = d1_rhs_et_dy_0_nxp4p01m4p1k*param_float(2)

d2_rhs_et_dydy_1_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p0k = ((q(nx+4+0,1-4+1-1+0,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1-1+0,indvars(2))*q(nx+4+0,1-4+1-1+0,indvars(2))+&
                    q(nx+4+0,1-4+1-1+0,indvars(3))*q(nx+4+0,1-4+1-1+0,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_1_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p1k = ((q(nx+4+0,1-4+1-1+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1-1+1,indvars(2))*q(nx+4+0,1-4+1-1+1,indvars(2))+&
                    q(nx+4+0,1-4+1-1+1,indvars(3))*q(nx+4+0,1-4+1-1+1,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_1_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p2k = ((q(nx+4+0,1-4+1-1+2,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1-1+2,indvars(2))*q(nx+4+0,1-4+1-1+2,indvars(2))+&
                    q(nx+4+0,1-4+1-1+2,indvars(3))*q(nx+4+0,1-4+1-1+2,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_1_0_nxp4p01m4p1m1k = -&
          1.5_wp*d2_rhs_et_dydy_1_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p0k+&
          2.0_wp*d2_rhs_et_dydy_1_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p1k-&
          0.5_wp*d2_rhs_et_dydy_1_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p2k

d2_rhs_et_dydy_1_0_nxp4p01m4p1m1k = d2_rhs_et_dydy_1_0_nxp4p01m4p1m1k*param_float(2)

d2_rhs_et_dydy_1_0_nxp4p01m4p1p1k_nxp4p01m4p1p1m1k = ((q(nx+4+0,1-4+1+1-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1+1-1,indvars(2))*q(nx+4+0,1-4+1+1-1,indvars(2))+&
                    q(nx+4+0,1-4+1+1-1,indvars(3))*q(nx+4+0,1-4+1+1-1,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_1_0_nxp4p01m4p1p1k_nxp4p01m4p1p1p1k = ((q(nx+4+0,1-4+1+1+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1+1+1,indvars(2))*q(nx+4+0,1-4+1+1+1,indvars(2))+&
                    q(nx+4+0,1-4+1+1+1,indvars(3))*q(nx+4+0,1-4+1+1+1,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_1_0_nxp4p01m4p1p1k = -&
          0.5_wp*d2_rhs_et_dydy_1_0_nxp4p01m4p1p1k_nxp4p01m4p1p1m1k+&
          0.5_wp*d2_rhs_et_dydy_1_0_nxp4p01m4p1p1k_nxp4p01m4p1p1p1k

d2_rhs_et_dydy_1_0_nxp4p01m4p1p1k = d2_rhs_et_dydy_1_0_nxp4p01m4p1p1k*param_float(2)

d1_rhs_et_dy_1_nxp4p01m4p1m1k = -param_float(2 + 5)*qst(nx+4+0,1-4+1-1,indvarsst(11))*(d2_rhs_et_dydy_1_0_nxp4p01m4p1m1k)-&
                    q(nx+4+0,1-4+1-1,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0,1-4+1-1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1-1,indvars(1)))**3.0_wp/((q(nx+4+0,1-4+1-1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1-1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0,1-4+1-1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1-1,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0,1-4+1-1,indvarsst(11))*(d2_rhs_u_dydy_1_0_nxp4p01m4p1m1k)+&
                    qst(nx+4+0,1-4+1-1,indvarsst(10))*(d2_rhs_u_dydx_1_0_nxp4p01m4p1m1k)))-&
                    q(nx+4+0,1-4+1-1,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0,1-4+1-1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1-1,indvars(1)))**3.0_wp/((q(nx+4+0,1-4+1-1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1-1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0,1-4+1-1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1-1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0,1-4+1-1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp4p01m4p1m1k)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0,1-4+1-1,indvarsst(10))*(d2_rhs_v_dydx_1_0_nxp4p01m4p1m1k)+&
                    qst(nx+4+0,1-4+1-1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp4p01m4p1m1k))))

d1_rhs_et_dy_1_nxp4p01m4p1p1k = -param_float(2 + 5)*qst(nx+4+0,1-4+1+1,indvarsst(11))*(d2_rhs_et_dydy_1_0_nxp4p01m4p1p1k)-&
                    q(nx+4+0,1-4+1+1,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0,1-4+1+1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1+1,indvars(1)))**3.0_wp/((q(nx+4+0,1-4+1+1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0,1-4+1+1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1+1,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0,1-4+1+1,indvarsst(11))*(d2_rhs_u_dydy_1_0_nxp4p01m4p1p1k)+&
                    qst(nx+4+0,1-4+1+1,indvarsst(10))*(d2_rhs_u_dydx_1_0_nxp4p01m4p1p1k)))-&
                    q(nx+4+0,1-4+1+1,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0,1-4+1+1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1+1,indvars(1)))**3.0_wp/((q(nx+4+0,1-4+1+1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0,1-4+1+1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1+1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0,1-4+1+1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp4p01m4p1p1k)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0,1-4+1+1,indvarsst(10))*(d2_rhs_v_dydx_1_0_nxp4p01m4p1p1k)+&
                    qst(nx+4+0,1-4+1+1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp4p01m4p1p1k))))

d1_rhs_et_dy_1_nxp4p01m4p1k = -&
          0.5_wp*d1_rhs_et_dy_1_nxp4p01m4p1m1k+&
          0.5_wp*d1_rhs_et_dy_1_nxp4p01m4p1p1k

d1_rhs_et_dy_1_nxp4p01m4p1k = d1_rhs_et_dy_1_nxp4p01m4p1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None d(rho et)/dt **********
!                                                           
!***********************************************************


rhs(nx+4+0,1-4+1,indvars(4)) =   -  ( 0.5_wp*(q(nx+4+0,1-4+1,indvars(2))**2+&
                    q(nx+4+0,1-4+1,indvars(3))**2)*(1.0_wp/(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5**2*((q(nx+4+0,1-4+1,indvars(2))*(-&
                    (param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5**2*d1_rhs_rho_dx_1_nxp4p01m4p1k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p01m4p1k))*qst(nx+4+0,1-4+1,indvarsst(10))+&
                    0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5+&
                    q(nx+4+0,1-4+1,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5*q(nx+4+0,1-4+1,indvars(1))*d1_rhs_rho_dx_3_nxp4p01m4p1k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p01m4p1k))*qst(nx+4+0,1-4+1,indvarsst(10))+&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5*(1-&
                    qface_i(1-4+1,2)*qface_i(1-4+1,2))/param_float(20 + 5)*((param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))-&
                    param_float(26 + 5)))))+&
                    1.0_wp/param_float(3 + 5)*(0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5+&
                    q(nx+4+0,1-4+1,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5*q(nx+4+0,1-4+1,indvars(1))*d1_rhs_rho_dx_3_nxp4p01m4p1k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p01m4p1k))*qst(nx+4+0,1-4+1,indvarsst(10))+&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5*(1-&
                    qface_i(1-4+1,2)*qface_i(1-4+1,2))/param_float(20 + 5)*((param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))-&
                    param_float(26 + 5))))+&
                    q(nx+4+0,1-4+1,indvars(1))*q(nx+4+0,1-4+1,indvars(2))*(1.0_wp/(2*q(nx+4+0,1-4+1,indvars(1))*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5)*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5+&
                    q(nx+4+0,1-4+1,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5*q(nx+4+0,1-4+1,indvars(1))*d1_rhs_rho_dx_3_nxp4p01m4p1k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p01m4p1k))*qst(nx+4+0,1-4+1,indvarsst(10))-&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))/q(nx+4+0,1-4+1,indvars(1)))**0.5*(1-&
                    qface_i(1-4+1,2)*qface_i(1-4+1,2))/param_float(20 + 5)*((param_float(3 + 5))*q(nx+4+0,1-4+1,indvars(1))*((q(nx+4+0,1-4+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+1,indvars(2))*q(nx+4+0,1-4+1,indvars(2))+&
                    q(nx+4+0,1-4+1,indvars(3))*q(nx+4+0,1-4+1,indvars(3)))))-&
                    param_float(26 + 5))))+&
                    q(nx+4+0,1-4+1,indvars(1))*q(nx+4+0,1-4+1,indvars(3))*((1.0_wp*q(nx+4+0,1-4+1,indvars(2))*d1_rhs_v_dx_4_nxp4p01m4p1k)*qst(nx+4+0,1-4+1,indvarsst(10)))+&
                    d1_rhs_et_dy_0_nxp4p01m4p1k+&
                    qst(nx+4+0,1-4+1,indvarsst(10))*(d1_rhs_et_dx_9_nxp4p01m4p1k)+&
                    qst(nx+4+0,1-4+1,indvarsst(11))*(d1_rhs_et_dy_1_nxp4p01m4p1k) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None d(rho nut)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ([rho*u*nut-visc_t*sigmaI*deltayI*({nut}_1y)]_1y)*deltayI+-ReI*Cb2*sigmaI*((deltayI)**2*([rho*nut]_1y)*([nut]_1y))-Cb1*(1-ft2)*SS*rho*nut+ReI*(Cw1*fw-Cb1/k**2*ft2)*rho*nut**2/eta**2
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_rhs_nut_dydy_0_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p0k = q(nx+4+0,1-4+1-1+0,indvars(5))

d2_rhs_nut_dydy_0_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p1k = q(nx+4+0,1-4+1-1+1,indvars(5))

d2_rhs_nut_dydy_0_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p2k = q(nx+4+0,1-4+1-1+2,indvars(5))

d2_rhs_nut_dydy_0_0_nxp4p01m4p1m1k = -&
          1.5_wp*d2_rhs_nut_dydy_0_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p0k+&
          2.0_wp*d2_rhs_nut_dydy_0_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p1k-&
          0.5_wp*d2_rhs_nut_dydy_0_0_nxp4p01m4p1m1k_nxp4p01m4p1m1p2k

d2_rhs_nut_dydy_0_0_nxp4p01m4p1m1k = d2_rhs_nut_dydy_0_0_nxp4p01m4p1m1k*param_float(2)

d2_rhs_nut_dydy_0_0_nxp4p01m4p1p1k_nxp4p01m4p1p1m1k = q(nx+4+0,1-4+1+1-1,indvars(5))

d2_rhs_nut_dydy_0_0_nxp4p01m4p1p1k_nxp4p01m4p1p1p1k = q(nx+4+0,1-4+1+1+1,indvars(5))

d2_rhs_nut_dydy_0_0_nxp4p01m4p1p1k = -&
          0.5_wp*d2_rhs_nut_dydy_0_0_nxp4p01m4p1p1k_nxp4p01m4p1p1m1k+&
          0.5_wp*d2_rhs_nut_dydy_0_0_nxp4p01m4p1p1k_nxp4p01m4p1p1p1k

d2_rhs_nut_dydy_0_0_nxp4p01m4p1p1k = d2_rhs_nut_dydy_0_0_nxp4p01m4p1p1k*param_float(2)

d1_rhs_nut_dy_1_nxp4p01m4p1m1k = q(nx+4+0,1-4+1-1,indvars(1))*q(nx+4+0,1-4+1-1,indvars(5))

d1_rhs_nut_dy_1_nxp4p01m4p1p1k = q(nx+4+0,1-4+1+1,indvars(1))*q(nx+4+0,1-4+1+1,indvars(5))

d1_rhs_nut_dy_1_nxp4p01m4p1k = -&
          0.5_wp*d1_rhs_nut_dy_1_nxp4p01m4p1m1k+&
          0.5_wp*d1_rhs_nut_dy_1_nxp4p01m4p1p1k

d1_rhs_nut_dy_1_nxp4p01m4p1k = d1_rhs_nut_dy_1_nxp4p01m4p1k*param_float(2)

d1_rhs_nut_dy_2_nxp4p01m4p1m1k = q(nx+4+0,1-4+1-1,indvars(5))

d1_rhs_nut_dy_2_nxp4p01m4p1p1k = q(nx+4+0,1-4+1+1,indvars(5))

d1_rhs_nut_dy_2_nxp4p01m4p1k = -&
          0.5_wp*d1_rhs_nut_dy_2_nxp4p01m4p1m1k+&
          0.5_wp*d1_rhs_nut_dy_2_nxp4p01m4p1p1k

d1_rhs_nut_dy_2_nxp4p01m4p1k = d1_rhs_nut_dy_2_nxp4p01m4p1k*param_float(2)

d1_rhs_nut_dy_0_nxp4p01m4p1m1k = q(nx+4+0,1-4+1-1,indvars(1))*q(nx+4+0,1-4+1-1,indvars(2))*q(nx+4+0,1-4+1-1,indvars(5))-&
                    (1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0,1-4+1-1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1-1,indvars(1)))**3.0_wp/((q(nx+4+0,1-4+1-1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1-1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0,1-4+1-1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1-1,indvars(1))))*param_float(1 + 5)*param_float(18 + 5)*qst(nx+4+0,1-4+1-1,indvarsst(11))*(d2_rhs_nut_dydy_0_0_nxp4p01m4p1m1k)

d1_rhs_nut_dy_0_nxp4p01m4p1p1k = q(nx+4+0,1-4+1+1,indvars(1))*q(nx+4+0,1-4+1+1,indvars(2))*q(nx+4+0,1-4+1+1,indvars(5))-&
                    (1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0,1-4+1+1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1+1,indvars(1)))**3.0_wp/((q(nx+4+0,1-4+1+1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0,1-4+1+1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1+1,indvars(1))))*param_float(1 + 5)*param_float(18 + 5)*qst(nx+4+0,1-4+1+1,indvarsst(11))*(d2_rhs_nut_dydy_0_0_nxp4p01m4p1p1k)

d1_rhs_nut_dy_0_nxp4p01m4p1k = -&
          0.5_wp*d1_rhs_nut_dy_0_nxp4p01m4p1m1k+&
          0.5_wp*d1_rhs_nut_dy_0_nxp4p01m4p1p1k

d1_rhs_nut_dy_0_nxp4p01m4p1k = d1_rhs_nut_dy_0_nxp4p01m4p1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None d(rho nut)/dt *********
!                                                           
!***********************************************************


rhs(nx+4+0,1-4+1,indvars(5)) =   -  ( (d1_rhs_nut_dy_0_nxp4p01m4p1k)*qst(nx+4+0,1-4+1,indvarsst(11))+&
                    -&
                    param_float(1 + 5)*param_float(7 + 5)*param_float(18 + 5)*((qst(nx+4+0,1-4+1,indvarsst(11)))**2*(d1_rhs_nut_dy_1_nxp4p01m4p1k)*(d1_rhs_nut_dy_2_nxp4p01m4p1k))-&
                    param_float(6 + 5)*(1-&
                    (param_float(16 + 5)*exp(-&
                    param_float(17 + 5)*(q(nx+4+0,1-4+1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1,indvars(1)))**2.0_wp)))*qst(nx+4+0,1-4+1,indvarsst(12))*q(nx+4+0,1-4+1,indvars(1))*q(nx+4+0,1-4+1,indvars(5))+&
                    param_float(1 + 5)*(param_float(10 + 5)*(((min(param_float(1 + 5)*(q(nx+4+0,1-4+1,indvars(5))/(qst(nx+4+0,1-4+1,indvarsst(12))*param_float(9 + 5)**2.0_wp*qst(nx+4+0,1-4+1,indvarsst(2))**2.0_wp)),10.0_wp))+&
                    param_float(11 + 5)*((min(param_float(1 + 5)*(q(nx+4+0,1-4+1,indvars(5))/(qst(nx+4+0,1-4+1,indvarsst(12))*param_float(9 + 5)**2.0_wp*qst(nx+4+0,1-4+1,indvarsst(2))**2.0_wp)),10.0_wp))**6.0_wp-&
                    (min(param_float(1 + 5)*(q(nx+4+0,1-4+1,indvars(5))/(qst(nx+4+0,1-4+1,indvarsst(12))*param_float(9 + 5)**2.0_wp*qst(nx+4+0,1-4+1,indvarsst(2))**2.0_wp)),10.0_wp))))*((1.0_wp+&
                    param_float(12 + 5)**6.0_wp)/(((min(param_float(1 + 5)*(q(nx+4+0,1-4+1,indvars(5))/(qst(nx+4+0,1-4+1,indvarsst(12))*param_float(9 + 5)**2.0_wp*qst(nx+4+0,1-4+1,indvarsst(2))**2.0_wp)),10.0_wp))+&
                    param_float(11 + 5)*((min(param_float(1 + 5)*(q(nx+4+0,1-4+1,indvars(5))/(qst(nx+4+0,1-4+1,indvarsst(12))*param_float(9 + 5)**2.0_wp*qst(nx+4+0,1-4+1,indvarsst(2))**2.0_wp)),10.0_wp))**6.0_wp-&
                    (min(param_float(1 + 5)*(q(nx+4+0,1-4+1,indvars(5))/(qst(nx+4+0,1-4+1,indvarsst(12))*param_float(9 + 5)**2.0_wp*qst(nx+4+0,1-4+1,indvarsst(2))**2.0_wp)),10.0_wp))))**6.0_wp+&
                    param_float(12 + 5)**6.0_wp))**(1.0_wp/6.0_wp))-&
                    param_float(6 + 5)/param_float(9 + 5)**2*(param_float(16 + 5)*exp(-&
                    param_float(17 + 5)*(q(nx+4+0,1-4+1,indvars(5))/1.0_wp*q(nx+4+0,1-4+1,indvars(1)))**2.0_wp)))*q(nx+4+0,1-4+1,indvars(1))*q(nx+4+0,1-4+1,indvars(5))**2/qst(nx+4+0,1-4+1,indvarsst(2))**2 ) 


    enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO

 end subroutine boundaryScheme_edges_PhyBC_Out_surf_rhs_imax_j1_0_1


