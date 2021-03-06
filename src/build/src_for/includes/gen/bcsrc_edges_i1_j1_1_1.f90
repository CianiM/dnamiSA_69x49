














!====================================================================================================
!
! General Boundary Conditions: compute Eqns with Boundary Scheme
!     
!====================================================================================================

subroutine boundaryScheme_edges_i1_j1_1_1(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

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



 real(wp) ::  d1_conv_rho_dx_0_1m2p1m11m2p1k,d1_conv_rho_dx_0_1m2p1p01m2p1k,d1_conv_rho_dx_0_1m2p1p11m2p1k &
            ,d1_conv_rho_dx_0_1m2p11m2p1k &
            ,d1_conv_rho_dy_0_1m2p11m2p1m1k,d1_conv_rho_dy_0_1m2p11m2p1p0k,d1_conv_rho_dy_0_1m2p11m2p1p1k &
            ,d1_conv_rho_dy_0_1m2p11m2p1k &
            ,d1_conv_rhou_dx_0_1m2p1m11m2p1k,d1_conv_rhou_dx_0_1m2p1p01m2p1k,d1_conv_rhou_dx_0_1m2p1p11m2p1k &
            ,d1_conv_rhou_dx_0_1m2p11m2p1k &
            ,d1_conv_rhou_dy_0_1m2p11m2p1m1k,d1_conv_rhou_dy_0_1m2p11m2p1p0k,d1_conv_rhou_dy_0_1m2p11m2p1p1k &
            ,d1_conv_rhou_dy_0_1m2p11m2p1k &
            ,d1_conv_rhov_dx_0_1m2p1m11m2p1k,d1_conv_rhov_dx_0_1m2p1p01m2p1k,d1_conv_rhov_dx_0_1m2p1p11m2p1k &
            ,d1_conv_rhov_dx_0_1m2p11m2p1k &
            ,d1_conv_rhov_dy_0_1m2p11m2p1m1k,d1_conv_rhov_dy_0_1m2p11m2p1p0k,d1_conv_rhov_dy_0_1m2p11m2p1p1k &
            ,d1_conv_rhov_dy_0_1m2p11m2p1k &
            ,d1_conv_rhoet_dx_0_1m2p1m11m2p1k,d1_conv_rhoet_dx_0_1m2p1p01m2p1k,d1_conv_rhoet_dx_0_1m2p1p11m2p1k &
            ,d1_conv_rhoet_dx_0_1m2p11m2p1k &
            ,d1_conv_rhoet_dy_0_1m2p11m2p1m1k,d1_conv_rhoet_dy_0_1m2p11m2p1p0k,d1_conv_rhoet_dy_0_1m2p11m2p1p1k &
            ,d1_conv_rhoet_dy_0_1m2p11m2p1k &
            ,d1_conv_rhonut_dx_0_1m2p1m11m2p1k,d1_conv_rhonut_dx_0_1m2p1p01m2p1k,d1_conv_rhonut_dx_0_1m2p1p11m2p1k &
            ,d1_conv_rhonut_dx_0_1m2p11m2p1k &
            ,d1_conv_rhonut_dy_0_1m2p11m2p1m1k,d1_conv_rhonut_dy_0_1m2p11m2p1p0k,d1_conv_rhonut_dy_0_1m2p11m2p1p1k &
            ,d1_conv_rhonut_dy_0_1m2p11m2p1k 

 real(wp) ::  d2_dif_rhou_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p01m2p1k,d2_dif_rhou_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p11m2p1k,d2_dif_rhou_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p21m2p1k &
            ,d2_dif_rhou_dxdx_0_0_1m2p1m11m2p1k &
            ,d2_dif_rhou_dxdx_0_0_1m2p1p11m2p1k_1m2p1p1m11m2p1k,d2_dif_rhou_dxdx_0_0_1m2p1p11m2p1k_1m2p1p1p01m2p1k,d2_dif_rhou_dxdx_0_0_1m2p1p11m2p1k_1m2p1p1p11m2p1k &
            ,d2_dif_rhou_dxdx_0_0_1m2p1p11m2p1k &
            ,d2_dif_rhou_dxdy_0_0_1m2p1m11m2p1k_1m2p1m11m2p1m1k,d2_dif_rhou_dxdy_0_0_1m2p1m11m2p1k_1m2p1m11m2p1p0k,d2_dif_rhou_dxdy_0_0_1m2p1m11m2p1k_1m2p1m11m2p1p1k &
            ,d2_dif_rhou_dxdy_0_0_1m2p1m11m2p1k &
            ,d2_dif_rhou_dxdy_0_0_1m2p1p11m2p1k_1m2p1p11m2p1m1k,d2_dif_rhou_dxdy_0_0_1m2p1p11m2p1k_1m2p1p11m2p1p0k,d2_dif_rhou_dxdy_0_0_1m2p1p11m2p1k_1m2p1p11m2p1p1k &
            ,d2_dif_rhou_dxdy_0_0_1m2p1p11m2p1k &
            ,d1_dif_rhou_dx_0_1m2p1m11m2p1k,d1_dif_rhou_dx_0_1m2p1p01m2p1k,d1_dif_rhou_dx_0_1m2p1p11m2p1k &
            ,d1_dif_rhou_dx_0_1m2p11m2p1k &
            ,d2_dif_rhou_dydx_0_0_1m2p11m2p1m1k_1m2p1m11m2p1m1k,d2_dif_rhou_dydx_0_0_1m2p11m2p1m1k_1m2p1p01m2p1m1k,d2_dif_rhou_dydx_0_0_1m2p11m2p1m1k_1m2p1p11m2p1m1k &
            ,d2_dif_rhou_dydx_0_0_1m2p11m2p1m1k &
            ,d2_dif_rhou_dydx_0_0_1m2p11m2p1p1k_1m2p1m11m2p1p1k,d2_dif_rhou_dydx_0_0_1m2p11m2p1p1k_1m2p1p01m2p1p1k,d2_dif_rhou_dydx_0_0_1m2p11m2p1p1k_1m2p1p11m2p1p1k &
            ,d2_dif_rhou_dydx_0_0_1m2p11m2p1p1k &
            ,d2_dif_rhou_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p0k,d2_dif_rhou_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p1k,d2_dif_rhou_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p2k &
            ,d2_dif_rhou_dydy_0_0_1m2p11m2p1m1k &
            ,d2_dif_rhou_dydy_0_0_1m2p11m2p1p1k_1m2p11m2p1p1m1k,d2_dif_rhou_dydy_0_0_1m2p11m2p1p1k_1m2p11m2p1p1p0k,d2_dif_rhou_dydy_0_0_1m2p11m2p1p1k_1m2p11m2p1p1p1k &
            ,d2_dif_rhou_dydy_0_0_1m2p11m2p1p1k &
            ,d1_dif_rhou_dy_0_1m2p11m2p1m1k,d1_dif_rhou_dy_0_1m2p11m2p1p0k,d1_dif_rhou_dy_0_1m2p11m2p1p1k &
            ,d1_dif_rhou_dy_0_1m2p11m2p1k &
            ,d2_dif_rhov_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p01m2p1k,d2_dif_rhov_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p11m2p1k,d2_dif_rhov_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p21m2p1k &
            ,d2_dif_rhov_dxdx_0_0_1m2p1m11m2p1k &
            ,d2_dif_rhov_dxdx_0_0_1m2p1p11m2p1k_1m2p1p1m11m2p1k,d2_dif_rhov_dxdx_0_0_1m2p1p11m2p1k_1m2p1p1p01m2p1k,d2_dif_rhov_dxdx_0_0_1m2p1p11m2p1k_1m2p1p1p11m2p1k &
            ,d2_dif_rhov_dxdx_0_0_1m2p1p11m2p1k &
            ,d2_dif_rhov_dxdy_0_0_1m2p1m11m2p1k_1m2p1m11m2p1m1k,d2_dif_rhov_dxdy_0_0_1m2p1m11m2p1k_1m2p1m11m2p1p0k,d2_dif_rhov_dxdy_0_0_1m2p1m11m2p1k_1m2p1m11m2p1p1k &
            ,d2_dif_rhov_dxdy_0_0_1m2p1m11m2p1k &
            ,d2_dif_rhov_dxdy_0_0_1m2p1p11m2p1k_1m2p1p11m2p1m1k,d2_dif_rhov_dxdy_0_0_1m2p1p11m2p1k_1m2p1p11m2p1p0k,d2_dif_rhov_dxdy_0_0_1m2p1p11m2p1k_1m2p1p11m2p1p1k &
            ,d2_dif_rhov_dxdy_0_0_1m2p1p11m2p1k &
            ,d1_dif_rhov_dx_0_1m2p1m11m2p1k,d1_dif_rhov_dx_0_1m2p1p01m2p1k,d1_dif_rhov_dx_0_1m2p1p11m2p1k &
            ,d1_dif_rhov_dx_0_1m2p11m2p1k &
            ,d2_dif_rhov_dydx_0_0_1m2p11m2p1m1k_1m2p1m11m2p1m1k,d2_dif_rhov_dydx_0_0_1m2p11m2p1m1k_1m2p1p01m2p1m1k,d2_dif_rhov_dydx_0_0_1m2p11m2p1m1k_1m2p1p11m2p1m1k &
            ,d2_dif_rhov_dydx_0_0_1m2p11m2p1m1k &
            ,d2_dif_rhov_dydx_0_0_1m2p11m2p1p1k_1m2p1m11m2p1p1k,d2_dif_rhov_dydx_0_0_1m2p11m2p1p1k_1m2p1p01m2p1p1k,d2_dif_rhov_dydx_0_0_1m2p11m2p1p1k_1m2p1p11m2p1p1k &
            ,d2_dif_rhov_dydx_0_0_1m2p11m2p1p1k &
            ,d2_dif_rhov_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p0k,d2_dif_rhov_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p1k,d2_dif_rhov_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p2k &
            ,d2_dif_rhov_dydy_0_0_1m2p11m2p1m1k &
            ,d2_dif_rhov_dydy_0_0_1m2p11m2p1p1k_1m2p11m2p1p1m1k,d2_dif_rhov_dydy_0_0_1m2p11m2p1p1k_1m2p11m2p1p1p0k,d2_dif_rhov_dydy_0_0_1m2p11m2p1p1k_1m2p11m2p1p1p1k &
            ,d2_dif_rhov_dydy_0_0_1m2p11m2p1p1k &
            ,d1_dif_rhov_dy_0_1m2p11m2p1m1k,d1_dif_rhov_dy_0_1m2p11m2p1p0k,d1_dif_rhov_dy_0_1m2p11m2p1p1k &
            ,d1_dif_rhov_dy_0_1m2p11m2p1k &
            ,d2_dif_rhoet_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p01m2p1k,d2_dif_rhoet_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p11m2p1k,d2_dif_rhoet_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p21m2p1k &
            ,d2_dif_rhoet_dxdx_0_0_1m2p1m11m2p1k &
            ,d2_dif_rhoet_dxdx_0_0_1m2p1p11m2p1k_1m2p1p1m11m2p1k,d2_dif_rhoet_dxdx_0_0_1m2p1p11m2p1k_1m2p1p1p01m2p1k,d2_dif_rhoet_dxdx_0_0_1m2p1p11m2p1k_1m2p1p1p11m2p1k &
            ,d2_dif_rhoet_dxdx_0_0_1m2p1p11m2p1k &
            ,d1_dif_rhoet_dx_0_1m2p1m11m2p1k,d1_dif_rhoet_dx_0_1m2p1p01m2p1k,d1_dif_rhoet_dx_0_1m2p1p11m2p1k &
            ,d1_dif_rhoet_dx_0_1m2p11m2p1k &
            ,d2_dif_rhoet_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p0k,d2_dif_rhoet_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p1k,d2_dif_rhoet_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p2k &
            ,d2_dif_rhoet_dydy_0_0_1m2p11m2p1m1k &
            ,d2_dif_rhoet_dydy_0_0_1m2p11m2p1p1k_1m2p11m2p1p1m1k,d2_dif_rhoet_dydy_0_0_1m2p11m2p1p1k_1m2p11m2p1p1p0k,d2_dif_rhoet_dydy_0_0_1m2p11m2p1p1k_1m2p11m2p1p1p1k &
            ,d2_dif_rhoet_dydy_0_0_1m2p11m2p1p1k &
            ,d1_dif_rhoet_dy_0_1m2p11m2p1m1k,d1_dif_rhoet_dy_0_1m2p11m2p1p0k,d1_dif_rhoet_dy_0_1m2p11m2p1p1k &
            ,d1_dif_rhoet_dy_0_1m2p11m2p1k &
            ,d2_dif_rhonut_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p01m2p1k,d2_dif_rhonut_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p11m2p1k,d2_dif_rhonut_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p21m2p1k &
            ,d2_dif_rhonut_dxdx_0_0_1m2p1m11m2p1k &
            ,d2_dif_rhonut_dxdx_0_0_1m2p1p11m2p1k_1m2p1p1m11m2p1k,d2_dif_rhonut_dxdx_0_0_1m2p1p11m2p1k_1m2p1p1p01m2p1k,d2_dif_rhonut_dxdx_0_0_1m2p1p11m2p1k_1m2p1p1p11m2p1k &
            ,d2_dif_rhonut_dxdx_0_0_1m2p1p11m2p1k &
            ,d1_dif_rhonut_dx_1_1m2p1m11m2p1k,d1_dif_rhonut_dx_1_1m2p1p01m2p1k,d1_dif_rhonut_dx_1_1m2p1p11m2p1k &
            ,d1_dif_rhonut_dx_1_1m2p11m2p1k &
            ,d1_dif_rhonut_dx_2_1m2p1m11m2p1k,d1_dif_rhonut_dx_2_1m2p1p01m2p1k,d1_dif_rhonut_dx_2_1m2p1p11m2p1k &
            ,d1_dif_rhonut_dx_2_1m2p11m2p1k &
            ,d1_dif_rhonut_dx_0_1m2p1m11m2p1k,d1_dif_rhonut_dx_0_1m2p1p01m2p1k,d1_dif_rhonut_dx_0_1m2p1p11m2p1k &
            ,d1_dif_rhonut_dx_0_1m2p11m2p1k &
            ,d2_dif_rhonut_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p0k,d2_dif_rhonut_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p1k,d2_dif_rhonut_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p2k &
            ,d2_dif_rhonut_dydy_0_0_1m2p11m2p1m1k &
            ,d2_dif_rhonut_dydy_0_0_1m2p11m2p1p1k_1m2p11m2p1p1m1k,d2_dif_rhonut_dydy_0_0_1m2p11m2p1p1k_1m2p11m2p1p1p0k,d2_dif_rhonut_dydy_0_0_1m2p11m2p1p1k_1m2p11m2p1p1p1k &
            ,d2_dif_rhonut_dydy_0_0_1m2p11m2p1p1k &
            ,d1_dif_rhonut_dy_1_1m2p11m2p1m1k,d1_dif_rhonut_dy_1_1m2p11m2p1p0k,d1_dif_rhonut_dy_1_1m2p11m2p1p1k &
            ,d1_dif_rhonut_dy_1_1m2p11m2p1k &
            ,d1_dif_rhonut_dy_2_1m2p11m2p1m1k,d1_dif_rhonut_dy_2_1m2p11m2p1p0k,d1_dif_rhonut_dy_2_1m2p11m2p1p1k &
            ,d1_dif_rhonut_dy_2_1m2p11m2p1k &
            ,d1_dif_rhonut_dy_0_1m2p11m2p1m1k,d1_dif_rhonut_dy_0_1m2p11m2p1p0k,d1_dif_rhonut_dy_0_1m2p11m2p1p1k &
            ,d1_dif_rhonut_dy_0_1m2p11m2p1k 

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
! BC layer: 1 1 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None d(rho)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([rho*u]_1x)+deltayI*([rho*v]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_conv_rho_dx_0_1m2p1m11m2p1k = q(1-2+1-1,1-2+1,indvars(1))*q(1-2+1-1,1-2+1,indvars(2))

d1_conv_rho_dx_0_1m2p1p11m2p1k = q(1-2+1+1,1-2+1,indvars(1))*q(1-2+1+1,1-2+1,indvars(2))

d1_conv_rho_dx_0_1m2p11m2p1k = -&
          0.5_wp*d1_conv_rho_dx_0_1m2p1m11m2p1k+&
          0.5_wp*d1_conv_rho_dx_0_1m2p1p11m2p1k

d1_conv_rho_dx_0_1m2p11m2p1k = d1_conv_rho_dx_0_1m2p11m2p1k*param_float(1)

d1_conv_rho_dy_0_1m2p11m2p1m1k = q(1-2+1,1-2+1-1,indvars(1))*q(1-2+1,1-2+1-1,indvars(3))

d1_conv_rho_dy_0_1m2p11m2p1p1k = q(1-2+1,1-2+1+1,indvars(1))*q(1-2+1,1-2+1+1,indvars(3))

d1_conv_rho_dy_0_1m2p11m2p1k = -&
          0.5_wp*d1_conv_rho_dy_0_1m2p11m2p1m1k+&
          0.5_wp*d1_conv_rho_dy_0_1m2p11m2p1p1k

d1_conv_rho_dy_0_1m2p11m2p1k = d1_conv_rho_dy_0_1m2p11m2p1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None d(rho)/dt *************
!                                                           
!***********************************************************


rhs(1-2+1,1-2+1,indvars(1)) =   -  ( qst(1-2+1,1-2+1,indvarsst(10))*(d1_conv_rho_dx_0_1m2p11m2p1k)+&
                    qst(1-2+1,1-2+1,indvarsst(11))*(d1_conv_rho_dy_0_1m2p11m2p1k) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None d(rho u)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([rho*u*u+p]_1x)+deltayI*([rho*v*u]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_conv_rhou_dx_0_1m2p1m11m2p1k = q(1-2+1-1,1-2+1,indvars(1))*q(1-2+1-1,1-2+1,indvars(2))*q(1-2+1-1,1-2+1,indvars(2))+(param_float(3 + 5))*q(1-2+1-1,1-2+1,indvars(1))*((q(1-2+1-1,1-2+1,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,1-2+1,indvars(2))*q(1-2+1-1,1-2+1,indvars(2))+&
                    q(1-2+1-1,1-2+1,indvars(3))*q(1-2+1-1,1-2+1,indvars(3)))))

d1_conv_rhou_dx_0_1m2p1p11m2p1k = q(1-2+1+1,1-2+1,indvars(1))*q(1-2+1+1,1-2+1,indvars(2))*q(1-2+1+1,1-2+1,indvars(2))+(param_float(3 + 5))*q(1-2+1+1,1-2+1,indvars(1))*((q(1-2+1+1,1-2+1,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,1-2+1,indvars(2))*q(1-2+1+1,1-2+1,indvars(2))+&
                    q(1-2+1+1,1-2+1,indvars(3))*q(1-2+1+1,1-2+1,indvars(3)))))

d1_conv_rhou_dx_0_1m2p11m2p1k = -&
          0.5_wp*d1_conv_rhou_dx_0_1m2p1m11m2p1k+&
          0.5_wp*d1_conv_rhou_dx_0_1m2p1p11m2p1k

d1_conv_rhou_dx_0_1m2p11m2p1k = d1_conv_rhou_dx_0_1m2p11m2p1k*param_float(1)

d1_conv_rhou_dy_0_1m2p11m2p1m1k = q(1-2+1,1-2+1-1,indvars(1))*q(1-2+1,1-2+1-1,indvars(3))*q(1-2+1,1-2+1-1,indvars(2))

d1_conv_rhou_dy_0_1m2p11m2p1p1k = q(1-2+1,1-2+1+1,indvars(1))*q(1-2+1,1-2+1+1,indvars(3))*q(1-2+1,1-2+1+1,indvars(2))

d1_conv_rhou_dy_0_1m2p11m2p1k = -&
          0.5_wp*d1_conv_rhou_dy_0_1m2p11m2p1m1k+&
          0.5_wp*d1_conv_rhou_dy_0_1m2p11m2p1p1k

d1_conv_rhou_dy_0_1m2p11m2p1k = d1_conv_rhou_dy_0_1m2p11m2p1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None d(rho u)/dt ***********
!                                                           
!***********************************************************


rhs(1-2+1,1-2+1,indvars(2)) =   -  ( qst(1-2+1,1-2+1,indvarsst(10))*(d1_conv_rhou_dx_0_1m2p11m2p1k)+&
                    qst(1-2+1,1-2+1,indvarsst(11))*(d1_conv_rhou_dy_0_1m2p11m2p1k) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None d(rho v)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([rho*u*v]_1x)+deltayI*([rho*v*v+p]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_conv_rhov_dx_0_1m2p1m11m2p1k = q(1-2+1-1,1-2+1,indvars(1))*q(1-2+1-1,1-2+1,indvars(2))*q(1-2+1-1,1-2+1,indvars(3))

d1_conv_rhov_dx_0_1m2p1p11m2p1k = q(1-2+1+1,1-2+1,indvars(1))*q(1-2+1+1,1-2+1,indvars(2))*q(1-2+1+1,1-2+1,indvars(3))

d1_conv_rhov_dx_0_1m2p11m2p1k = -&
          0.5_wp*d1_conv_rhov_dx_0_1m2p1m11m2p1k+&
          0.5_wp*d1_conv_rhov_dx_0_1m2p1p11m2p1k

d1_conv_rhov_dx_0_1m2p11m2p1k = d1_conv_rhov_dx_0_1m2p11m2p1k*param_float(1)

d1_conv_rhov_dy_0_1m2p11m2p1m1k = q(1-2+1,1-2+1-1,indvars(1))*q(1-2+1,1-2+1-1,indvars(3))*q(1-2+1,1-2+1-1,indvars(3))+(param_float(3 + 5))*q(1-2+1,1-2+1-1,indvars(1))*((q(1-2+1,1-2+1-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,1-2+1-1,indvars(2))*q(1-2+1,1-2+1-1,indvars(2))+&
                    q(1-2+1,1-2+1-1,indvars(3))*q(1-2+1,1-2+1-1,indvars(3)))))

d1_conv_rhov_dy_0_1m2p11m2p1p1k = q(1-2+1,1-2+1+1,indvars(1))*q(1-2+1,1-2+1+1,indvars(3))*q(1-2+1,1-2+1+1,indvars(3))+(param_float(3 + 5))*q(1-2+1,1-2+1+1,indvars(1))*((q(1-2+1,1-2+1+1,indvars(4))-&
                    0.5_wp*(q(1-2+1,1-2+1+1,indvars(2))*q(1-2+1,1-2+1+1,indvars(2))+&
                    q(1-2+1,1-2+1+1,indvars(3))*q(1-2+1,1-2+1+1,indvars(3)))))

d1_conv_rhov_dy_0_1m2p11m2p1k = -&
          0.5_wp*d1_conv_rhov_dy_0_1m2p11m2p1m1k+&
          0.5_wp*d1_conv_rhov_dy_0_1m2p11m2p1p1k

d1_conv_rhov_dy_0_1m2p11m2p1k = d1_conv_rhov_dy_0_1m2p11m2p1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None d(rho v)/dt ***********
!                                                           
!***********************************************************


rhs(1-2+1,1-2+1,indvars(3)) =   -  ( qst(1-2+1,1-2+1,indvarsst(10))*(d1_conv_rhov_dx_0_1m2p11m2p1k)+&
                    qst(1-2+1,1-2+1,indvarsst(11))*(d1_conv_rhov_dy_0_1m2p11m2p1k) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None d(rho et)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([(rho*et+p)*u]_1x)+deltayI*([(rho*et+p)*v]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_conv_rhoet_dx_0_1m2p1m11m2p1k = (q(1-2+1-1,1-2+1,indvars(1))*q(1-2+1-1,1-2+1,indvars(4))+&
                    (param_float(3 + 5))*q(1-2+1-1,1-2+1,indvars(1))*((q(1-2+1-1,1-2+1,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,1-2+1,indvars(2))*q(1-2+1-1,1-2+1,indvars(2))+&
                    q(1-2+1-1,1-2+1,indvars(3))*q(1-2+1-1,1-2+1,indvars(3))))))*q(1-2+1-1,1-2+1,indvars(2))

d1_conv_rhoet_dx_0_1m2p1p11m2p1k = (q(1-2+1+1,1-2+1,indvars(1))*q(1-2+1+1,1-2+1,indvars(4))+&
                    (param_float(3 + 5))*q(1-2+1+1,1-2+1,indvars(1))*((q(1-2+1+1,1-2+1,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,1-2+1,indvars(2))*q(1-2+1+1,1-2+1,indvars(2))+&
                    q(1-2+1+1,1-2+1,indvars(3))*q(1-2+1+1,1-2+1,indvars(3))))))*q(1-2+1+1,1-2+1,indvars(2))

d1_conv_rhoet_dx_0_1m2p11m2p1k = -&
          0.5_wp*d1_conv_rhoet_dx_0_1m2p1m11m2p1k+&
          0.5_wp*d1_conv_rhoet_dx_0_1m2p1p11m2p1k

d1_conv_rhoet_dx_0_1m2p11m2p1k = d1_conv_rhoet_dx_0_1m2p11m2p1k*param_float(1)

d1_conv_rhoet_dy_0_1m2p11m2p1m1k = (q(1-2+1,1-2+1-1,indvars(1))*q(1-2+1,1-2+1-1,indvars(4))+&
                    (param_float(3 + 5))*q(1-2+1,1-2+1-1,indvars(1))*((q(1-2+1,1-2+1-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,1-2+1-1,indvars(2))*q(1-2+1,1-2+1-1,indvars(2))+&
                    q(1-2+1,1-2+1-1,indvars(3))*q(1-2+1,1-2+1-1,indvars(3))))))*q(1-2+1,1-2+1-1,indvars(3))

d1_conv_rhoet_dy_0_1m2p11m2p1p1k = (q(1-2+1,1-2+1+1,indvars(1))*q(1-2+1,1-2+1+1,indvars(4))+&
                    (param_float(3 + 5))*q(1-2+1,1-2+1+1,indvars(1))*((q(1-2+1,1-2+1+1,indvars(4))-&
                    0.5_wp*(q(1-2+1,1-2+1+1,indvars(2))*q(1-2+1,1-2+1+1,indvars(2))+&
                    q(1-2+1,1-2+1+1,indvars(3))*q(1-2+1,1-2+1+1,indvars(3))))))*q(1-2+1,1-2+1+1,indvars(3))

d1_conv_rhoet_dy_0_1m2p11m2p1k = -&
          0.5_wp*d1_conv_rhoet_dy_0_1m2p11m2p1m1k+&
          0.5_wp*d1_conv_rhoet_dy_0_1m2p11m2p1p1k

d1_conv_rhoet_dy_0_1m2p11m2p1k = d1_conv_rhoet_dy_0_1m2p11m2p1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None d(rho et)/dt **********
!                                                           
!***********************************************************


rhs(1-2+1,1-2+1,indvars(4)) =   -  ( qst(1-2+1,1-2+1,indvarsst(10))*(d1_conv_rhoet_dx_0_1m2p11m2p1k)+&
                    qst(1-2+1,1-2+1,indvarsst(11))*(d1_conv_rhoet_dy_0_1m2p11m2p1k) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None d(rho nut)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([(rho*nut)*u]_1x)+deltayI*([(rho*nut)*v]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_conv_rhonut_dx_0_1m2p1m11m2p1k = (q(1-2+1-1,1-2+1,indvars(1))*q(1-2+1-1,1-2+1,indvars(5)))*q(1-2+1-1,1-2+1,indvars(2))

d1_conv_rhonut_dx_0_1m2p1p11m2p1k = (q(1-2+1+1,1-2+1,indvars(1))*q(1-2+1+1,1-2+1,indvars(5)))*q(1-2+1+1,1-2+1,indvars(2))

d1_conv_rhonut_dx_0_1m2p11m2p1k = -&
          0.5_wp*d1_conv_rhonut_dx_0_1m2p1m11m2p1k+&
          0.5_wp*d1_conv_rhonut_dx_0_1m2p1p11m2p1k

d1_conv_rhonut_dx_0_1m2p11m2p1k = d1_conv_rhonut_dx_0_1m2p11m2p1k*param_float(1)

d1_conv_rhonut_dy_0_1m2p11m2p1m1k = (q(1-2+1,1-2+1-1,indvars(1))*q(1-2+1,1-2+1-1,indvars(5)))*q(1-2+1,1-2+1-1,indvars(3))

d1_conv_rhonut_dy_0_1m2p11m2p1p1k = (q(1-2+1,1-2+1+1,indvars(1))*q(1-2+1,1-2+1+1,indvars(5)))*q(1-2+1,1-2+1+1,indvars(3))

d1_conv_rhonut_dy_0_1m2p11m2p1k = -&
          0.5_wp*d1_conv_rhonut_dy_0_1m2p11m2p1m1k+&
          0.5_wp*d1_conv_rhonut_dy_0_1m2p11m2p1p1k

d1_conv_rhonut_dy_0_1m2p11m2p1k = d1_conv_rhonut_dy_0_1m2p11m2p1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None d(rho nut)/dt *********
!                                                           
!***********************************************************


rhs(1-2+1,1-2+1,indvars(5)) =   -  ( qst(1-2+1,1-2+1,indvarsst(10))*(d1_conv_rhonut_dx_0_1m2p11m2p1k)+&
                    qst(1-2+1,1-2+1,indvarsst(11))*(d1_conv_rhonut_dy_0_1m2p11m2p1k) ) 



!***********************************************************
!                                                           
! Start building layers for BC : i1 j1 None ****************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 1 1 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None d(rho u)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([-visc_t*(2.0_wp*deltaxI*({u}_1x)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y)))]_1x)+deltayI*([-visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_dif_rhou_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p01m2p1k = q(1-2+1-1+0,1-2+1,indvars(2))

d2_dif_rhou_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p11m2p1k = q(1-2+1-1+1,1-2+1,indvars(2))

d2_dif_rhou_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p21m2p1k = q(1-2+1-1+2,1-2+1,indvars(2))

d2_dif_rhou_dxdx_0_0_1m2p1m11m2p1k = -&
          1.5_wp*d2_dif_rhou_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p01m2p1k+&
          2.0_wp*d2_dif_rhou_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p11m2p1k-&
          0.5_wp*d2_dif_rhou_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p21m2p1k

d2_dif_rhou_dxdx_0_0_1m2p1m11m2p1k = d2_dif_rhou_dxdx_0_0_1m2p1m11m2p1k*param_float(1)

d2_dif_rhou_dxdx_0_0_1m2p1p11m2p1k_1m2p1p1m11m2p1k = q(1-2+1+1-1,1-2+1,indvars(2))

d2_dif_rhou_dxdx_0_0_1m2p1p11m2p1k_1m2p1p1p11m2p1k = q(1-2+1+1+1,1-2+1,indvars(2))

d2_dif_rhou_dxdx_0_0_1m2p1p11m2p1k = -&
          0.5_wp*d2_dif_rhou_dxdx_0_0_1m2p1p11m2p1k_1m2p1p1m11m2p1k+&
          0.5_wp*d2_dif_rhou_dxdx_0_0_1m2p1p11m2p1k_1m2p1p1p11m2p1k

d2_dif_rhou_dxdx_0_0_1m2p1p11m2p1k = d2_dif_rhou_dxdx_0_0_1m2p1p11m2p1k*param_float(1)

d2_dif_rhou_dxdy_0_0_1m2p1m11m2p1k_1m2p1m11m2p1m1k = q(1-2+1-1,1-2+1-1,indvars(3))

d2_dif_rhou_dxdy_0_0_1m2p1m11m2p1k_1m2p1m11m2p1p1k = q(1-2+1-1,1-2+1+1,indvars(3))

d2_dif_rhou_dxdy_0_0_1m2p1m11m2p1k = -&
          0.5_wp*d2_dif_rhou_dxdy_0_0_1m2p1m11m2p1k_1m2p1m11m2p1m1k+&
          0.5_wp*d2_dif_rhou_dxdy_0_0_1m2p1m11m2p1k_1m2p1m11m2p1p1k

d2_dif_rhou_dxdy_0_0_1m2p1m11m2p1k = d2_dif_rhou_dxdy_0_0_1m2p1m11m2p1k*param_float(2)

d2_dif_rhou_dxdy_0_0_1m2p1p11m2p1k_1m2p1p11m2p1m1k = q(1-2+1+1,1-2+1-1,indvars(3))

d2_dif_rhou_dxdy_0_0_1m2p1p11m2p1k_1m2p1p11m2p1p1k = q(1-2+1+1,1-2+1+1,indvars(3))

d2_dif_rhou_dxdy_0_0_1m2p1p11m2p1k = -&
          0.5_wp*d2_dif_rhou_dxdy_0_0_1m2p1p11m2p1k_1m2p1p11m2p1m1k+&
          0.5_wp*d2_dif_rhou_dxdy_0_0_1m2p1p11m2p1k_1m2p1p11m2p1p1k

d2_dif_rhou_dxdy_0_0_1m2p1p11m2p1k = d2_dif_rhou_dxdy_0_0_1m2p1p11m2p1k*param_float(2)

d1_dif_rhou_dx_0_1m2p1m11m2p1k = -(1.0_wp)*(1.0_wp+&
                    ((q(1-2+1-1,1-2+1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1-1,1-2+1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1-1,1-2+1,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(1-2+1-1,1-2+1,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_1m2p1m11m2p1k)-&
                    2.0_wp/3.0_wp*(qst(1-2+1-1,1-2+1,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_1m2p1m11m2p1k)+&
                    qst(1-2+1-1,1-2+1,indvarsst(11))*(d2_dif_rhou_dxdy_0_0_1m2p1m11m2p1k)))

d1_dif_rhou_dx_0_1m2p1p11m2p1k = -(1.0_wp)*(1.0_wp+&
                    ((q(1-2+1+1,1-2+1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1+1,1-2+1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1+1,1-2+1,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(1-2+1+1,1-2+1,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_1m2p1p11m2p1k)-&
                    2.0_wp/3.0_wp*(qst(1-2+1+1,1-2+1,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_1m2p1p11m2p1k)+&
                    qst(1-2+1+1,1-2+1,indvarsst(11))*(d2_dif_rhou_dxdy_0_0_1m2p1p11m2p1k)))

d1_dif_rhou_dx_0_1m2p11m2p1k = -&
          0.5_wp*d1_dif_rhou_dx_0_1m2p1m11m2p1k+&
          0.5_wp*d1_dif_rhou_dx_0_1m2p1p11m2p1k

d1_dif_rhou_dx_0_1m2p11m2p1k = d1_dif_rhou_dx_0_1m2p11m2p1k*param_float(1)

d2_dif_rhou_dydx_0_0_1m2p11m2p1m1k_1m2p1m11m2p1m1k = q(1-2+1-1,1-2+1-1,indvars(3))

d2_dif_rhou_dydx_0_0_1m2p11m2p1m1k_1m2p1p11m2p1m1k = q(1-2+1+1,1-2+1-1,indvars(3))

d2_dif_rhou_dydx_0_0_1m2p11m2p1m1k = -&
          0.5_wp*d2_dif_rhou_dydx_0_0_1m2p11m2p1m1k_1m2p1m11m2p1m1k+&
          0.5_wp*d2_dif_rhou_dydx_0_0_1m2p11m2p1m1k_1m2p1p11m2p1m1k

d2_dif_rhou_dydx_0_0_1m2p11m2p1m1k = d2_dif_rhou_dydx_0_0_1m2p11m2p1m1k*param_float(1)

d2_dif_rhou_dydx_0_0_1m2p11m2p1p1k_1m2p1m11m2p1p1k = q(1-2+1-1,1-2+1+1,indvars(3))

d2_dif_rhou_dydx_0_0_1m2p11m2p1p1k_1m2p1p11m2p1p1k = q(1-2+1+1,1-2+1+1,indvars(3))

d2_dif_rhou_dydx_0_0_1m2p11m2p1p1k = -&
          0.5_wp*d2_dif_rhou_dydx_0_0_1m2p11m2p1p1k_1m2p1m11m2p1p1k+&
          0.5_wp*d2_dif_rhou_dydx_0_0_1m2p11m2p1p1k_1m2p1p11m2p1p1k

d2_dif_rhou_dydx_0_0_1m2p11m2p1p1k = d2_dif_rhou_dydx_0_0_1m2p11m2p1p1k*param_float(1)

d2_dif_rhou_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p0k = q(1-2+1,1-2+1-1+0,indvars(2))

d2_dif_rhou_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p1k = q(1-2+1,1-2+1-1+1,indvars(2))

d2_dif_rhou_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p2k = q(1-2+1,1-2+1-1+2,indvars(2))

d2_dif_rhou_dydy_0_0_1m2p11m2p1m1k = -&
          1.5_wp*d2_dif_rhou_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p0k+&
          2.0_wp*d2_dif_rhou_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p1k-&
          0.5_wp*d2_dif_rhou_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p2k

d2_dif_rhou_dydy_0_0_1m2p11m2p1m1k = d2_dif_rhou_dydy_0_0_1m2p11m2p1m1k*param_float(2)

d2_dif_rhou_dydy_0_0_1m2p11m2p1p1k_1m2p11m2p1p1m1k = q(1-2+1,1-2+1+1-1,indvars(2))

d2_dif_rhou_dydy_0_0_1m2p11m2p1p1k_1m2p11m2p1p1p1k = q(1-2+1,1-2+1+1+1,indvars(2))

d2_dif_rhou_dydy_0_0_1m2p11m2p1p1k = -&
          0.5_wp*d2_dif_rhou_dydy_0_0_1m2p11m2p1p1k_1m2p11m2p1p1m1k+&
          0.5_wp*d2_dif_rhou_dydy_0_0_1m2p11m2p1p1k_1m2p11m2p1p1p1k

d2_dif_rhou_dydy_0_0_1m2p11m2p1p1k = d2_dif_rhou_dydy_0_0_1m2p11m2p1p1k*param_float(2)

d1_dif_rhou_dy_0_1m2p11m2p1m1k = -(1.0_wp)*(1.0_wp+&
                    ((q(1-2+1,1-2+1-1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1,1-2+1-1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1,1-2+1-1,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(1-2+1,1-2+1-1,indvarsst(11))*(d2_dif_rhou_dydy_0_0_1m2p11m2p1m1k)+&
                    qst(1-2+1,1-2+1-1,indvarsst(10))*(d2_dif_rhou_dydx_0_0_1m2p11m2p1m1k))

d1_dif_rhou_dy_0_1m2p11m2p1p1k = -(1.0_wp)*(1.0_wp+&
                    ((q(1-2+1,1-2+1+1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1,1-2+1+1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1,1-2+1+1,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(1-2+1,1-2+1+1,indvarsst(11))*(d2_dif_rhou_dydy_0_0_1m2p11m2p1p1k)+&
                    qst(1-2+1,1-2+1+1,indvarsst(10))*(d2_dif_rhou_dydx_0_0_1m2p11m2p1p1k))

d1_dif_rhou_dy_0_1m2p11m2p1k = -&
          0.5_wp*d1_dif_rhou_dy_0_1m2p11m2p1m1k+&
          0.5_wp*d1_dif_rhou_dy_0_1m2p11m2p1p1k

d1_dif_rhou_dy_0_1m2p11m2p1k = d1_dif_rhou_dy_0_1m2p11m2p1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None d(rho u)/dt ***********
!                                                           
!***********************************************************


rhs(1-2+1,1-2+1,indvars(2)) = rhs(1-2+1,1-2+1,indvars(2))  -  ( qst(1-2+1,1-2+1,indvarsst(10))*(d1_dif_rhou_dx_0_1m2p11m2p1k)+&
                    qst(1-2+1,1-2+1,indvarsst(11))*(d1_dif_rhou_dy_0_1m2p11m2p1k) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None d(rho v)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([-visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x))]_1x)+deltayI*([-visc_t*(2.0_wp*deltayI*({v}_1y)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y)))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_dif_rhov_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p01m2p1k = q(1-2+1-1+0,1-2+1,indvars(3))

d2_dif_rhov_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p11m2p1k = q(1-2+1-1+1,1-2+1,indvars(3))

d2_dif_rhov_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p21m2p1k = q(1-2+1-1+2,1-2+1,indvars(3))

d2_dif_rhov_dxdx_0_0_1m2p1m11m2p1k = -&
          1.5_wp*d2_dif_rhov_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p01m2p1k+&
          2.0_wp*d2_dif_rhov_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p11m2p1k-&
          0.5_wp*d2_dif_rhov_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p21m2p1k

d2_dif_rhov_dxdx_0_0_1m2p1m11m2p1k = d2_dif_rhov_dxdx_0_0_1m2p1m11m2p1k*param_float(1)

d2_dif_rhov_dxdx_0_0_1m2p1p11m2p1k_1m2p1p1m11m2p1k = q(1-2+1+1-1,1-2+1,indvars(3))

d2_dif_rhov_dxdx_0_0_1m2p1p11m2p1k_1m2p1p1p11m2p1k = q(1-2+1+1+1,1-2+1,indvars(3))

d2_dif_rhov_dxdx_0_0_1m2p1p11m2p1k = -&
          0.5_wp*d2_dif_rhov_dxdx_0_0_1m2p1p11m2p1k_1m2p1p1m11m2p1k+&
          0.5_wp*d2_dif_rhov_dxdx_0_0_1m2p1p11m2p1k_1m2p1p1p11m2p1k

d2_dif_rhov_dxdx_0_0_1m2p1p11m2p1k = d2_dif_rhov_dxdx_0_0_1m2p1p11m2p1k*param_float(1)

d2_dif_rhov_dxdy_0_0_1m2p1m11m2p1k_1m2p1m11m2p1m1k = q(1-2+1-1,1-2+1-1,indvars(2))

d2_dif_rhov_dxdy_0_0_1m2p1m11m2p1k_1m2p1m11m2p1p1k = q(1-2+1-1,1-2+1+1,indvars(2))

d2_dif_rhov_dxdy_0_0_1m2p1m11m2p1k = -&
          0.5_wp*d2_dif_rhov_dxdy_0_0_1m2p1m11m2p1k_1m2p1m11m2p1m1k+&
          0.5_wp*d2_dif_rhov_dxdy_0_0_1m2p1m11m2p1k_1m2p1m11m2p1p1k

d2_dif_rhov_dxdy_0_0_1m2p1m11m2p1k = d2_dif_rhov_dxdy_0_0_1m2p1m11m2p1k*param_float(2)

d2_dif_rhov_dxdy_0_0_1m2p1p11m2p1k_1m2p1p11m2p1m1k = q(1-2+1+1,1-2+1-1,indvars(2))

d2_dif_rhov_dxdy_0_0_1m2p1p11m2p1k_1m2p1p11m2p1p1k = q(1-2+1+1,1-2+1+1,indvars(2))

d2_dif_rhov_dxdy_0_0_1m2p1p11m2p1k = -&
          0.5_wp*d2_dif_rhov_dxdy_0_0_1m2p1p11m2p1k_1m2p1p11m2p1m1k+&
          0.5_wp*d2_dif_rhov_dxdy_0_0_1m2p1p11m2p1k_1m2p1p11m2p1p1k

d2_dif_rhov_dxdy_0_0_1m2p1p11m2p1k = d2_dif_rhov_dxdy_0_0_1m2p1p11m2p1k*param_float(2)

d1_dif_rhov_dx_0_1m2p1m11m2p1k = -(1.0_wp)*(1.0_wp+&
                    ((q(1-2+1-1,1-2+1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1-1,1-2+1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1-1,1-2+1,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(1-2+1-1,1-2+1,indvarsst(11))*(d2_dif_rhov_dxdy_0_0_1m2p1m11m2p1k)+&
                    qst(1-2+1-1,1-2+1,indvarsst(10))*(d2_dif_rhov_dxdx_0_0_1m2p1m11m2p1k))

d1_dif_rhov_dx_0_1m2p1p11m2p1k = -(1.0_wp)*(1.0_wp+&
                    ((q(1-2+1+1,1-2+1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1+1,1-2+1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1+1,1-2+1,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(1-2+1+1,1-2+1,indvarsst(11))*(d2_dif_rhov_dxdy_0_0_1m2p1p11m2p1k)+&
                    qst(1-2+1+1,1-2+1,indvarsst(10))*(d2_dif_rhov_dxdx_0_0_1m2p1p11m2p1k))

d1_dif_rhov_dx_0_1m2p11m2p1k = -&
          0.5_wp*d1_dif_rhov_dx_0_1m2p1m11m2p1k+&
          0.5_wp*d1_dif_rhov_dx_0_1m2p1p11m2p1k

d1_dif_rhov_dx_0_1m2p11m2p1k = d1_dif_rhov_dx_0_1m2p11m2p1k*param_float(1)

d2_dif_rhov_dydx_0_0_1m2p11m2p1m1k_1m2p1m11m2p1m1k = q(1-2+1-1,1-2+1-1,indvars(2))

d2_dif_rhov_dydx_0_0_1m2p11m2p1m1k_1m2p1p11m2p1m1k = q(1-2+1+1,1-2+1-1,indvars(2))

d2_dif_rhov_dydx_0_0_1m2p11m2p1m1k = -&
          0.5_wp*d2_dif_rhov_dydx_0_0_1m2p11m2p1m1k_1m2p1m11m2p1m1k+&
          0.5_wp*d2_dif_rhov_dydx_0_0_1m2p11m2p1m1k_1m2p1p11m2p1m1k

d2_dif_rhov_dydx_0_0_1m2p11m2p1m1k = d2_dif_rhov_dydx_0_0_1m2p11m2p1m1k*param_float(1)

d2_dif_rhov_dydx_0_0_1m2p11m2p1p1k_1m2p1m11m2p1p1k = q(1-2+1-1,1-2+1+1,indvars(2))

d2_dif_rhov_dydx_0_0_1m2p11m2p1p1k_1m2p1p11m2p1p1k = q(1-2+1+1,1-2+1+1,indvars(2))

d2_dif_rhov_dydx_0_0_1m2p11m2p1p1k = -&
          0.5_wp*d2_dif_rhov_dydx_0_0_1m2p11m2p1p1k_1m2p1m11m2p1p1k+&
          0.5_wp*d2_dif_rhov_dydx_0_0_1m2p11m2p1p1k_1m2p1p11m2p1p1k

d2_dif_rhov_dydx_0_0_1m2p11m2p1p1k = d2_dif_rhov_dydx_0_0_1m2p11m2p1p1k*param_float(1)

d2_dif_rhov_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p0k = q(1-2+1,1-2+1-1+0,indvars(3))

d2_dif_rhov_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p1k = q(1-2+1,1-2+1-1+1,indvars(3))

d2_dif_rhov_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p2k = q(1-2+1,1-2+1-1+2,indvars(3))

d2_dif_rhov_dydy_0_0_1m2p11m2p1m1k = -&
          1.5_wp*d2_dif_rhov_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p0k+&
          2.0_wp*d2_dif_rhov_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p1k-&
          0.5_wp*d2_dif_rhov_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p2k

d2_dif_rhov_dydy_0_0_1m2p11m2p1m1k = d2_dif_rhov_dydy_0_0_1m2p11m2p1m1k*param_float(2)

d2_dif_rhov_dydy_0_0_1m2p11m2p1p1k_1m2p11m2p1p1m1k = q(1-2+1,1-2+1+1-1,indvars(3))

d2_dif_rhov_dydy_0_0_1m2p11m2p1p1k_1m2p11m2p1p1p1k = q(1-2+1,1-2+1+1+1,indvars(3))

d2_dif_rhov_dydy_0_0_1m2p11m2p1p1k = -&
          0.5_wp*d2_dif_rhov_dydy_0_0_1m2p11m2p1p1k_1m2p11m2p1p1m1k+&
          0.5_wp*d2_dif_rhov_dydy_0_0_1m2p11m2p1p1k_1m2p11m2p1p1p1k

d2_dif_rhov_dydy_0_0_1m2p11m2p1p1k = d2_dif_rhov_dydy_0_0_1m2p11m2p1p1k*param_float(2)

d1_dif_rhov_dy_0_1m2p11m2p1m1k = -(1.0_wp)*(1.0_wp+&
                    ((q(1-2+1,1-2+1-1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1,1-2+1-1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1,1-2+1-1,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(1-2+1,1-2+1-1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_1m2p11m2p1m1k)-&
                    2.0_wp/3.0_wp*(qst(1-2+1,1-2+1-1,indvarsst(10))*(d2_dif_rhov_dydx_0_0_1m2p11m2p1m1k)+&
                    qst(1-2+1,1-2+1-1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_1m2p11m2p1m1k)))

d1_dif_rhov_dy_0_1m2p11m2p1p1k = -(1.0_wp)*(1.0_wp+&
                    ((q(1-2+1,1-2+1+1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1,1-2+1+1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1,1-2+1+1,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(1-2+1,1-2+1+1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_1m2p11m2p1p1k)-&
                    2.0_wp/3.0_wp*(qst(1-2+1,1-2+1+1,indvarsst(10))*(d2_dif_rhov_dydx_0_0_1m2p11m2p1p1k)+&
                    qst(1-2+1,1-2+1+1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_1m2p11m2p1p1k)))

d1_dif_rhov_dy_0_1m2p11m2p1k = -&
          0.5_wp*d1_dif_rhov_dy_0_1m2p11m2p1m1k+&
          0.5_wp*d1_dif_rhov_dy_0_1m2p11m2p1p1k

d1_dif_rhov_dy_0_1m2p11m2p1k = d1_dif_rhov_dy_0_1m2p11m2p1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None d(rho v)/dt ***********
!                                                           
!***********************************************************


rhs(1-2+1,1-2+1,indvars(3)) = rhs(1-2+1,1-2+1,indvars(3))  -  ( qst(1-2+1,1-2+1,indvarsst(10))*(d1_dif_rhov_dx_0_1m2p11m2p1k)+&
                    qst(1-2+1,1-2+1,indvarsst(11))*(d1_dif_rhov_dy_0_1m2p11m2p1k) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None d(rho et)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([-kappa*deltaxI*({T}_1x)-u*(visc_t*(2.0_wp*deltaxI*({u}_1x)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y))))-v*(visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x)))]_1x)+deltayI*([-kappa*deltayI*({T}_1y)-u*(visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x)))-v*(visc_t*(2.0_wp*deltayI*({v}_1y)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y))))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_dif_rhoet_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p01m2p1k = ((q(1-2+1-1+0,1-2+1,indvars(4))-&
                    0.5_wp*(q(1-2+1-1+0,1-2+1,indvars(2))*q(1-2+1-1+0,1-2+1,indvars(2))+&
                    q(1-2+1-1+0,1-2+1,indvars(3))*q(1-2+1-1+0,1-2+1,indvars(3)))))/(param_float(4 + 5))

d2_dif_rhoet_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p11m2p1k = ((q(1-2+1-1+1,1-2+1,indvars(4))-&
                    0.5_wp*(q(1-2+1-1+1,1-2+1,indvars(2))*q(1-2+1-1+1,1-2+1,indvars(2))+&
                    q(1-2+1-1+1,1-2+1,indvars(3))*q(1-2+1-1+1,1-2+1,indvars(3)))))/(param_float(4 + 5))

d2_dif_rhoet_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p21m2p1k = ((q(1-2+1-1+2,1-2+1,indvars(4))-&
                    0.5_wp*(q(1-2+1-1+2,1-2+1,indvars(2))*q(1-2+1-1+2,1-2+1,indvars(2))+&
                    q(1-2+1-1+2,1-2+1,indvars(3))*q(1-2+1-1+2,1-2+1,indvars(3)))))/(param_float(4 + 5))

d2_dif_rhoet_dxdx_0_0_1m2p1m11m2p1k = -&
          1.5_wp*d2_dif_rhoet_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p01m2p1k+&
          2.0_wp*d2_dif_rhoet_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p11m2p1k-&
          0.5_wp*d2_dif_rhoet_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p21m2p1k

d2_dif_rhoet_dxdx_0_0_1m2p1m11m2p1k = d2_dif_rhoet_dxdx_0_0_1m2p1m11m2p1k*param_float(1)

d2_dif_rhoet_dxdx_0_0_1m2p1p11m2p1k_1m2p1p1m11m2p1k = ((q(1-2+1+1-1,1-2+1,indvars(4))-&
                    0.5_wp*(q(1-2+1+1-1,1-2+1,indvars(2))*q(1-2+1+1-1,1-2+1,indvars(2))+&
                    q(1-2+1+1-1,1-2+1,indvars(3))*q(1-2+1+1-1,1-2+1,indvars(3)))))/(param_float(4 + 5))

d2_dif_rhoet_dxdx_0_0_1m2p1p11m2p1k_1m2p1p1p11m2p1k = ((q(1-2+1+1+1,1-2+1,indvars(4))-&
                    0.5_wp*(q(1-2+1+1+1,1-2+1,indvars(2))*q(1-2+1+1+1,1-2+1,indvars(2))+&
                    q(1-2+1+1+1,1-2+1,indvars(3))*q(1-2+1+1+1,1-2+1,indvars(3)))))/(param_float(4 + 5))

d2_dif_rhoet_dxdx_0_0_1m2p1p11m2p1k = -&
          0.5_wp*d2_dif_rhoet_dxdx_0_0_1m2p1p11m2p1k_1m2p1p1m11m2p1k+&
          0.5_wp*d2_dif_rhoet_dxdx_0_0_1m2p1p11m2p1k_1m2p1p1p11m2p1k

d2_dif_rhoet_dxdx_0_0_1m2p1p11m2p1k = d2_dif_rhoet_dxdx_0_0_1m2p1p11m2p1k*param_float(1)

d1_dif_rhoet_dx_0_1m2p1m11m2p1k = -param_float(2 + 5)*qst(1-2+1-1,1-2+1,indvarsst(10))*(d2_dif_rhoet_dxdx_0_0_1m2p1m11m2p1k)-&
                    q(1-2+1-1,1-2+1,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(1-2+1-1,1-2+1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1-1,1-2+1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1-1,1-2+1,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(1-2+1-1,1-2+1,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_1m2p1m11m2p1k)-&
                    2.0_wp/3.0_wp*(qst(1-2+1-1,1-2+1,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_1m2p1m11m2p1k)+&
                    qst(1-2+1-1,1-2+1,indvarsst(11))*(d2_dif_rhou_dxdy_0_0_1m2p1m11m2p1k))))-&
                    q(1-2+1-1,1-2+1,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(1-2+1-1,1-2+1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1-1,1-2+1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1-1,1-2+1,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(1-2+1-1,1-2+1,indvarsst(11))*(d2_dif_rhov_dxdy_0_0_1m2p1m11m2p1k)+&
                    qst(1-2+1-1,1-2+1,indvarsst(10))*(d2_dif_rhov_dxdx_0_0_1m2p1m11m2p1k)))

d1_dif_rhoet_dx_0_1m2p1p11m2p1k = -param_float(2 + 5)*qst(1-2+1+1,1-2+1,indvarsst(10))*(d2_dif_rhoet_dxdx_0_0_1m2p1p11m2p1k)-&
                    q(1-2+1+1,1-2+1,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(1-2+1+1,1-2+1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1+1,1-2+1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1+1,1-2+1,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(1-2+1+1,1-2+1,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_1m2p1p11m2p1k)-&
                    2.0_wp/3.0_wp*(qst(1-2+1+1,1-2+1,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_1m2p1p11m2p1k)+&
                    qst(1-2+1+1,1-2+1,indvarsst(11))*(d2_dif_rhou_dxdy_0_0_1m2p1p11m2p1k))))-&
                    q(1-2+1+1,1-2+1,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(1-2+1+1,1-2+1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1+1,1-2+1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1+1,1-2+1,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(1-2+1+1,1-2+1,indvarsst(11))*(d2_dif_rhov_dxdy_0_0_1m2p1p11m2p1k)+&
                    qst(1-2+1+1,1-2+1,indvarsst(10))*(d2_dif_rhov_dxdx_0_0_1m2p1p11m2p1k)))

d1_dif_rhoet_dx_0_1m2p11m2p1k = -&
          0.5_wp*d1_dif_rhoet_dx_0_1m2p1m11m2p1k+&
          0.5_wp*d1_dif_rhoet_dx_0_1m2p1p11m2p1k

d1_dif_rhoet_dx_0_1m2p11m2p1k = d1_dif_rhoet_dx_0_1m2p11m2p1k*param_float(1)

d2_dif_rhoet_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p0k = ((q(1-2+1,1-2+1-1+0,indvars(4))-&
                    0.5_wp*(q(1-2+1,1-2+1-1+0,indvars(2))*q(1-2+1,1-2+1-1+0,indvars(2))+&
                    q(1-2+1,1-2+1-1+0,indvars(3))*q(1-2+1,1-2+1-1+0,indvars(3)))))/(param_float(4 + 5))

d2_dif_rhoet_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p1k = ((q(1-2+1,1-2+1-1+1,indvars(4))-&
                    0.5_wp*(q(1-2+1,1-2+1-1+1,indvars(2))*q(1-2+1,1-2+1-1+1,indvars(2))+&
                    q(1-2+1,1-2+1-1+1,indvars(3))*q(1-2+1,1-2+1-1+1,indvars(3)))))/(param_float(4 + 5))

d2_dif_rhoet_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p2k = ((q(1-2+1,1-2+1-1+2,indvars(4))-&
                    0.5_wp*(q(1-2+1,1-2+1-1+2,indvars(2))*q(1-2+1,1-2+1-1+2,indvars(2))+&
                    q(1-2+1,1-2+1-1+2,indvars(3))*q(1-2+1,1-2+1-1+2,indvars(3)))))/(param_float(4 + 5))

d2_dif_rhoet_dydy_0_0_1m2p11m2p1m1k = -&
          1.5_wp*d2_dif_rhoet_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p0k+&
          2.0_wp*d2_dif_rhoet_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p1k-&
          0.5_wp*d2_dif_rhoet_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p2k

d2_dif_rhoet_dydy_0_0_1m2p11m2p1m1k = d2_dif_rhoet_dydy_0_0_1m2p11m2p1m1k*param_float(2)

d2_dif_rhoet_dydy_0_0_1m2p11m2p1p1k_1m2p11m2p1p1m1k = ((q(1-2+1,1-2+1+1-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,1-2+1+1-1,indvars(2))*q(1-2+1,1-2+1+1-1,indvars(2))+&
                    q(1-2+1,1-2+1+1-1,indvars(3))*q(1-2+1,1-2+1+1-1,indvars(3)))))/(param_float(4 + 5))

d2_dif_rhoet_dydy_0_0_1m2p11m2p1p1k_1m2p11m2p1p1p1k = ((q(1-2+1,1-2+1+1+1,indvars(4))-&
                    0.5_wp*(q(1-2+1,1-2+1+1+1,indvars(2))*q(1-2+1,1-2+1+1+1,indvars(2))+&
                    q(1-2+1,1-2+1+1+1,indvars(3))*q(1-2+1,1-2+1+1+1,indvars(3)))))/(param_float(4 + 5))

d2_dif_rhoet_dydy_0_0_1m2p11m2p1p1k = -&
          0.5_wp*d2_dif_rhoet_dydy_0_0_1m2p11m2p1p1k_1m2p11m2p1p1m1k+&
          0.5_wp*d2_dif_rhoet_dydy_0_0_1m2p11m2p1p1k_1m2p11m2p1p1p1k

d2_dif_rhoet_dydy_0_0_1m2p11m2p1p1k = d2_dif_rhoet_dydy_0_0_1m2p11m2p1p1k*param_float(2)

d1_dif_rhoet_dy_0_1m2p11m2p1m1k = -param_float(2 + 5)*qst(1-2+1,1-2+1-1,indvarsst(11))*(d2_dif_rhoet_dydy_0_0_1m2p11m2p1m1k)-&
                    q(1-2+1,1-2+1-1,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(1-2+1,1-2+1-1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1,1-2+1-1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1,1-2+1-1,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(1-2+1,1-2+1-1,indvarsst(11))*(d2_dif_rhou_dydy_0_0_1m2p11m2p1m1k)+&
                    qst(1-2+1,1-2+1-1,indvarsst(10))*(d2_dif_rhou_dydx_0_0_1m2p11m2p1m1k)))-&
                    q(1-2+1,1-2+1-1,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(1-2+1,1-2+1-1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1,1-2+1-1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1,1-2+1-1,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(1-2+1,1-2+1-1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_1m2p11m2p1m1k)-&
                    2.0_wp/3.0_wp*(qst(1-2+1,1-2+1-1,indvarsst(10))*(d2_dif_rhov_dydx_0_0_1m2p11m2p1m1k)+&
                    qst(1-2+1,1-2+1-1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_1m2p11m2p1m1k))))

d1_dif_rhoet_dy_0_1m2p11m2p1p1k = -param_float(2 + 5)*qst(1-2+1,1-2+1+1,indvarsst(11))*(d2_dif_rhoet_dydy_0_0_1m2p11m2p1p1k)-&
                    q(1-2+1,1-2+1+1,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(1-2+1,1-2+1+1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1,1-2+1+1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1,1-2+1+1,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(1-2+1,1-2+1+1,indvarsst(11))*(d2_dif_rhou_dydy_0_0_1m2p11m2p1p1k)+&
                    qst(1-2+1,1-2+1+1,indvarsst(10))*(d2_dif_rhou_dydx_0_0_1m2p11m2p1p1k)))-&
                    q(1-2+1,1-2+1+1,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(1-2+1,1-2+1+1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1,1-2+1+1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1,1-2+1+1,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(1-2+1,1-2+1+1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_1m2p11m2p1p1k)-&
                    2.0_wp/3.0_wp*(qst(1-2+1,1-2+1+1,indvarsst(10))*(d2_dif_rhov_dydx_0_0_1m2p11m2p1p1k)+&
                    qst(1-2+1,1-2+1+1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_1m2p11m2p1p1k))))

d1_dif_rhoet_dy_0_1m2p11m2p1k = -&
          0.5_wp*d1_dif_rhoet_dy_0_1m2p11m2p1m1k+&
          0.5_wp*d1_dif_rhoet_dy_0_1m2p11m2p1p1k

d1_dif_rhoet_dy_0_1m2p11m2p1k = d1_dif_rhoet_dy_0_1m2p11m2p1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None d(rho et)/dt **********
!                                                           
!***********************************************************


rhs(1-2+1,1-2+1,indvars(4)) = rhs(1-2+1,1-2+1,indvars(4))  -  ( qst(1-2+1,1-2+1,indvarsst(10))*(d1_dif_rhoet_dx_0_1m2p11m2p1k)+&
                    qst(1-2+1,1-2+1,indvarsst(11))*(d1_dif_rhoet_dy_0_1m2p11m2p1k) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None d(rho nut)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([-(ReI*(1.0_wp+chi)*sigmaI*deltaxI*({nut}_1x))]_1x)+deltayI*([-(ReI*(1.0_wp+chi)*sigmaI*deltayI*({nut}_1y))]_1y)-ReI*Cb2*sigmaI*((deltaxI)**2*([rho*nut]_1x)*([nut]_1x)+(deltayI)**2*([rho*nut]_1y)*([nut]_1y))-Cb1*(1.0_wp-ft2)*SS*rho*nut+ReI*(Cw1*fw-Cb1/k**2.0_wp*ft2)*rho*(nut/eta)**2.0_wp
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_dif_rhonut_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p01m2p1k = q(1-2+1-1+0,1-2+1,indvars(5))

d2_dif_rhonut_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p11m2p1k = q(1-2+1-1+1,1-2+1,indvars(5))

d2_dif_rhonut_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p21m2p1k = q(1-2+1-1+2,1-2+1,indvars(5))

d2_dif_rhonut_dxdx_0_0_1m2p1m11m2p1k = -&
          1.5_wp*d2_dif_rhonut_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p01m2p1k+&
          2.0_wp*d2_dif_rhonut_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p11m2p1k-&
          0.5_wp*d2_dif_rhonut_dxdx_0_0_1m2p1m11m2p1k_1m2p1m1p21m2p1k

d2_dif_rhonut_dxdx_0_0_1m2p1m11m2p1k = d2_dif_rhonut_dxdx_0_0_1m2p1m11m2p1k*param_float(1)

d2_dif_rhonut_dxdx_0_0_1m2p1p11m2p1k_1m2p1p1m11m2p1k = q(1-2+1+1-1,1-2+1,indvars(5))

d2_dif_rhonut_dxdx_0_0_1m2p1p11m2p1k_1m2p1p1p11m2p1k = q(1-2+1+1+1,1-2+1,indvars(5))

d2_dif_rhonut_dxdx_0_0_1m2p1p11m2p1k = -&
          0.5_wp*d2_dif_rhonut_dxdx_0_0_1m2p1p11m2p1k_1m2p1p1m11m2p1k+&
          0.5_wp*d2_dif_rhonut_dxdx_0_0_1m2p1p11m2p1k_1m2p1p1p11m2p1k

d2_dif_rhonut_dxdx_0_0_1m2p1p11m2p1k = d2_dif_rhonut_dxdx_0_0_1m2p1p11m2p1k*param_float(1)

d1_dif_rhonut_dx_1_1m2p1m11m2p1k = q(1-2+1-1,1-2+1,indvars(1))*q(1-2+1-1,1-2+1,indvars(5))

d1_dif_rhonut_dx_1_1m2p1p11m2p1k = q(1-2+1+1,1-2+1,indvars(1))*q(1-2+1+1,1-2+1,indvars(5))

d1_dif_rhonut_dx_1_1m2p11m2p1k = -&
          0.5_wp*d1_dif_rhonut_dx_1_1m2p1m11m2p1k+&
          0.5_wp*d1_dif_rhonut_dx_1_1m2p1p11m2p1k

d1_dif_rhonut_dx_1_1m2p11m2p1k = d1_dif_rhonut_dx_1_1m2p11m2p1k*param_float(1)

d1_dif_rhonut_dx_2_1m2p1m11m2p1k = q(1-2+1-1,1-2+1,indvars(5))

d1_dif_rhonut_dx_2_1m2p1p11m2p1k = q(1-2+1+1,1-2+1,indvars(5))

d1_dif_rhonut_dx_2_1m2p11m2p1k = -&
          0.5_wp*d1_dif_rhonut_dx_2_1m2p1m11m2p1k+&
          0.5_wp*d1_dif_rhonut_dx_2_1m2p1p11m2p1k

d1_dif_rhonut_dx_2_1m2p11m2p1k = d1_dif_rhonut_dx_2_1m2p11m2p1k*param_float(1)

d1_dif_rhonut_dx_0_1m2p1m11m2p1k = -(param_float(1 + 5)*(1.0_wp+&
                    (q(1-2+1-1,1-2+1,indvars(5))/1.0_wp))*param_float(18 + 5)*qst(1-2+1-1,1-2+1,indvarsst(10))*(d2_dif_rhonut_dxdx_0_0_1m2p1m11m2p1k))

d1_dif_rhonut_dx_0_1m2p1p11m2p1k = -(param_float(1 + 5)*(1.0_wp+&
                    (q(1-2+1+1,1-2+1,indvars(5))/1.0_wp))*param_float(18 + 5)*qst(1-2+1+1,1-2+1,indvarsst(10))*(d2_dif_rhonut_dxdx_0_0_1m2p1p11m2p1k))

d1_dif_rhonut_dx_0_1m2p11m2p1k = -&
          0.5_wp*d1_dif_rhonut_dx_0_1m2p1m11m2p1k+&
          0.5_wp*d1_dif_rhonut_dx_0_1m2p1p11m2p1k

d1_dif_rhonut_dx_0_1m2p11m2p1k = d1_dif_rhonut_dx_0_1m2p11m2p1k*param_float(1)

d2_dif_rhonut_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p0k = q(1-2+1,1-2+1-1+0,indvars(5))

d2_dif_rhonut_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p1k = q(1-2+1,1-2+1-1+1,indvars(5))

d2_dif_rhonut_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p2k = q(1-2+1,1-2+1-1+2,indvars(5))

d2_dif_rhonut_dydy_0_0_1m2p11m2p1m1k = -&
          1.5_wp*d2_dif_rhonut_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p0k+&
          2.0_wp*d2_dif_rhonut_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p1k-&
          0.5_wp*d2_dif_rhonut_dydy_0_0_1m2p11m2p1m1k_1m2p11m2p1m1p2k

d2_dif_rhonut_dydy_0_0_1m2p11m2p1m1k = d2_dif_rhonut_dydy_0_0_1m2p11m2p1m1k*param_float(2)

d2_dif_rhonut_dydy_0_0_1m2p11m2p1p1k_1m2p11m2p1p1m1k = q(1-2+1,1-2+1+1-1,indvars(5))

d2_dif_rhonut_dydy_0_0_1m2p11m2p1p1k_1m2p11m2p1p1p1k = q(1-2+1,1-2+1+1+1,indvars(5))

d2_dif_rhonut_dydy_0_0_1m2p11m2p1p1k = -&
          0.5_wp*d2_dif_rhonut_dydy_0_0_1m2p11m2p1p1k_1m2p11m2p1p1m1k+&
          0.5_wp*d2_dif_rhonut_dydy_0_0_1m2p11m2p1p1k_1m2p11m2p1p1p1k

d2_dif_rhonut_dydy_0_0_1m2p11m2p1p1k = d2_dif_rhonut_dydy_0_0_1m2p11m2p1p1k*param_float(2)

d1_dif_rhonut_dy_1_1m2p11m2p1m1k = q(1-2+1,1-2+1-1,indvars(1))*q(1-2+1,1-2+1-1,indvars(5))

d1_dif_rhonut_dy_1_1m2p11m2p1p1k = q(1-2+1,1-2+1+1,indvars(1))*q(1-2+1,1-2+1+1,indvars(5))

d1_dif_rhonut_dy_1_1m2p11m2p1k = -&
          0.5_wp*d1_dif_rhonut_dy_1_1m2p11m2p1m1k+&
          0.5_wp*d1_dif_rhonut_dy_1_1m2p11m2p1p1k

d1_dif_rhonut_dy_1_1m2p11m2p1k = d1_dif_rhonut_dy_1_1m2p11m2p1k*param_float(2)

d1_dif_rhonut_dy_2_1m2p11m2p1m1k = q(1-2+1,1-2+1-1,indvars(5))

d1_dif_rhonut_dy_2_1m2p11m2p1p1k = q(1-2+1,1-2+1+1,indvars(5))

d1_dif_rhonut_dy_2_1m2p11m2p1k = -&
          0.5_wp*d1_dif_rhonut_dy_2_1m2p11m2p1m1k+&
          0.5_wp*d1_dif_rhonut_dy_2_1m2p11m2p1p1k

d1_dif_rhonut_dy_2_1m2p11m2p1k = d1_dif_rhonut_dy_2_1m2p11m2p1k*param_float(2)

d1_dif_rhonut_dy_0_1m2p11m2p1m1k = -(param_float(1 + 5)*(1.0_wp+&
                    (q(1-2+1,1-2+1-1,indvars(5))/1.0_wp))*param_float(18 + 5)*qst(1-2+1,1-2+1-1,indvarsst(11))*(d2_dif_rhonut_dydy_0_0_1m2p11m2p1m1k))

d1_dif_rhonut_dy_0_1m2p11m2p1p1k = -(param_float(1 + 5)*(1.0_wp+&
                    (q(1-2+1,1-2+1+1,indvars(5))/1.0_wp))*param_float(18 + 5)*qst(1-2+1,1-2+1+1,indvarsst(11))*(d2_dif_rhonut_dydy_0_0_1m2p11m2p1p1k))

d1_dif_rhonut_dy_0_1m2p11m2p1k = -&
          0.5_wp*d1_dif_rhonut_dy_0_1m2p11m2p1m1k+&
          0.5_wp*d1_dif_rhonut_dy_0_1m2p11m2p1p1k

d1_dif_rhonut_dy_0_1m2p11m2p1k = d1_dif_rhonut_dy_0_1m2p11m2p1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None d(rho nut)/dt *********
!                                                           
!***********************************************************


rhs(1-2+1,1-2+1,indvars(5)) = rhs(1-2+1,1-2+1,indvars(5))  -  ( qst(1-2+1,1-2+1,indvarsst(10))*(d1_dif_rhonut_dx_0_1m2p11m2p1k)+&
                    qst(1-2+1,1-2+1,indvarsst(11))*(d1_dif_rhonut_dy_0_1m2p11m2p1k)-&
                    param_float(1 + 5)*param_float(7 + 5)*param_float(18 + 5)*((qst(1-2+1,1-2+1,indvarsst(10)))**2*(d1_dif_rhonut_dx_1_1m2p11m2p1k)*(d1_dif_rhonut_dx_2_1m2p11m2p1k)+&
                    (qst(1-2+1,1-2+1,indvarsst(11)))**2*(d1_dif_rhonut_dy_1_1m2p11m2p1k)*(d1_dif_rhonut_dy_2_1m2p11m2p1k))-&
                    param_float(6 + 5)*(1.0_wp-&
                    (param_float(16 + 5)*exp(-&
                    param_float(17 + 5)*(q(1-2+1,1-2+1,indvars(5))/1.0_wp)**2.0_wp)))*qst(1-2+1,1-2+1,indvarsst(13))*q(1-2+1,1-2+1,indvars(1))*q(1-2+1,1-2+1,indvars(5))+&
                    param_float(1 + 5)*(param_float(10 + 5)*qst(1-2+1,1-2+1,indvarsst(16))-&
                    param_float(6 + 5)/param_float(9 + 5)**2.0_wp*(param_float(16 + 5)*exp(-&
                    param_float(17 + 5)*(q(1-2+1,1-2+1,indvars(5))/1.0_wp)**2.0_wp)))*q(1-2+1,1-2+1,indvars(1))*(q(1-2+1,1-2+1,indvars(5))/qst(1-2+1,1-2+1,indvarsst(2)))**2.0_wp ) 


    enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO

 end subroutine boundaryScheme_edges_i1_j1_1_1


