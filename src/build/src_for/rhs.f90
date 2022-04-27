












subroutine cmprhs(param_float,ind,idloop,idarray,neq,neqst,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

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


integer, intent(in) :: neq,neqst
integer, intent(in), dimension(6)   :: idloop, idarray
integer, intent(in), dimension(neq+neqst) :: ind
integer, intent(in), dimension(3) :: sizeblck
integer, intent(in) :: nvar_f(3),nvar_e(3)

real(wp), intent(in)    :: param_float(*)

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



 real(wp) ::  d1_conv_rho_dx_0_im2jk,d1_conv_rho_dx_0_im1jk,d1_conv_rho_dx_0_ip0jk,d1_conv_rho_dx_0_ip1jk,d1_conv_rho_dx_0_ip2jk &
            ,d1_conv_rho_dx_0_ijk &
            ,d1_conv_rho_dy_0_ijm2k,d1_conv_rho_dy_0_ijm1k,d1_conv_rho_dy_0_ijp0k,d1_conv_rho_dy_0_ijp1k,d1_conv_rho_dy_0_ijp2k &
            ,d1_conv_rho_dy_0_ijk &
            ,d1_conv_rhou_dx_0_im2jk,d1_conv_rhou_dx_0_im1jk,d1_conv_rhou_dx_0_ip0jk,d1_conv_rhou_dx_0_ip1jk,d1_conv_rhou_dx_0_ip2jk &
            ,d1_conv_rhou_dx_0_ijk &
            ,d1_conv_rhou_dy_0_ijm2k,d1_conv_rhou_dy_0_ijm1k,d1_conv_rhou_dy_0_ijp0k,d1_conv_rhou_dy_0_ijp1k,d1_conv_rhou_dy_0_ijp2k &
            ,d1_conv_rhou_dy_0_ijk &
            ,d1_conv_rhov_dx_0_im2jk,d1_conv_rhov_dx_0_im1jk,d1_conv_rhov_dx_0_ip0jk,d1_conv_rhov_dx_0_ip1jk,d1_conv_rhov_dx_0_ip2jk &
            ,d1_conv_rhov_dx_0_ijk &
            ,d1_conv_rhov_dy_0_ijm2k,d1_conv_rhov_dy_0_ijm1k,d1_conv_rhov_dy_0_ijp0k,d1_conv_rhov_dy_0_ijp1k,d1_conv_rhov_dy_0_ijp2k &
            ,d1_conv_rhov_dy_0_ijk &
            ,d1_conv_et_dx_0_im2jk,d1_conv_et_dx_0_im1jk,d1_conv_et_dx_0_ip0jk,d1_conv_et_dx_0_ip1jk,d1_conv_et_dx_0_ip2jk &
            ,d1_conv_et_dx_0_ijk &
            ,d1_conv_et_dy_0_ijm2k,d1_conv_et_dy_0_ijm1k,d1_conv_et_dy_0_ijp0k,d1_conv_et_dy_0_ijp1k,d1_conv_et_dy_0_ijp2k &
            ,d1_conv_et_dy_0_ijk &
            ,d2_conv_nut_dxdx_0_0_im2jk_im2m2jk,d2_conv_nut_dxdx_0_0_im2jk_im2m1jk,d2_conv_nut_dxdx_0_0_im2jk_im2p0jk,d2_conv_nut_dxdx_0_0_im2jk_im2p1jk,d2_conv_nut_dxdx_0_0_im2jk_im2p2jk &
            ,d2_conv_nut_dxdx_0_0_im2jk &
            ,d2_conv_nut_dxdx_0_0_im1jk_im1m2jk,d2_conv_nut_dxdx_0_0_im1jk_im1m1jk,d2_conv_nut_dxdx_0_0_im1jk_im1p0jk,d2_conv_nut_dxdx_0_0_im1jk_im1p1jk,d2_conv_nut_dxdx_0_0_im1jk_im1p2jk &
            ,d2_conv_nut_dxdx_0_0_im1jk &
            ,d2_conv_nut_dxdx_0_0_ip1jk_ip1m2jk,d2_conv_nut_dxdx_0_0_ip1jk_ip1m1jk,d2_conv_nut_dxdx_0_0_ip1jk_ip1p0jk,d2_conv_nut_dxdx_0_0_ip1jk_ip1p1jk,d2_conv_nut_dxdx_0_0_ip1jk_ip1p2jk &
            ,d2_conv_nut_dxdx_0_0_ip1jk &
            ,d2_conv_nut_dxdx_0_0_ip2jk_ip2m2jk,d2_conv_nut_dxdx_0_0_ip2jk_ip2m1jk,d2_conv_nut_dxdx_0_0_ip2jk_ip2p0jk,d2_conv_nut_dxdx_0_0_ip2jk_ip2p1jk,d2_conv_nut_dxdx_0_0_ip2jk_ip2p2jk &
            ,d2_conv_nut_dxdx_0_0_ip2jk &
            ,d1_conv_nut_dx_0_im2jk,d1_conv_nut_dx_0_im1jk,d1_conv_nut_dx_0_ip0jk,d1_conv_nut_dx_0_ip1jk,d1_conv_nut_dx_0_ip2jk &
            ,d1_conv_nut_dx_0_ijk &
            ,d2_conv_nut_dydy_0_0_ijm2k_ijm2m2k,d2_conv_nut_dydy_0_0_ijm2k_ijm2m1k,d2_conv_nut_dydy_0_0_ijm2k_ijm2p0k,d2_conv_nut_dydy_0_0_ijm2k_ijm2p1k,d2_conv_nut_dydy_0_0_ijm2k_ijm2p2k &
            ,d2_conv_nut_dydy_0_0_ijm2k &
            ,d2_conv_nut_dydy_0_0_ijm1k_ijm1m2k,d2_conv_nut_dydy_0_0_ijm1k_ijm1m1k,d2_conv_nut_dydy_0_0_ijm1k_ijm1p0k,d2_conv_nut_dydy_0_0_ijm1k_ijm1p1k,d2_conv_nut_dydy_0_0_ijm1k_ijm1p2k &
            ,d2_conv_nut_dydy_0_0_ijm1k &
            ,d2_conv_nut_dydy_0_0_ijp1k_ijp1m2k,d2_conv_nut_dydy_0_0_ijp1k_ijp1m1k,d2_conv_nut_dydy_0_0_ijp1k_ijp1p0k,d2_conv_nut_dydy_0_0_ijp1k_ijp1p1k,d2_conv_nut_dydy_0_0_ijp1k_ijp1p2k &
            ,d2_conv_nut_dydy_0_0_ijp1k &
            ,d2_conv_nut_dydy_0_0_ijp2k_ijp2m2k,d2_conv_nut_dydy_0_0_ijp2k_ijp2m1k,d2_conv_nut_dydy_0_0_ijp2k_ijp2p0k,d2_conv_nut_dydy_0_0_ijp2k_ijp2p1k,d2_conv_nut_dydy_0_0_ijp2k_ijp2p2k &
            ,d2_conv_nut_dydy_0_0_ijp2k &
            ,d1_conv_nut_dy_0_ijm2k,d1_conv_nut_dy_0_ijm1k,d1_conv_nut_dy_0_ijp0k,d1_conv_nut_dy_0_ijp1k,d1_conv_nut_dy_0_ijp2k &
            ,d1_conv_nut_dy_0_ijk 

 real(wp) ::  d2_dif_rhou_dxdx_0_0_im1jk_im1m1jk,d2_dif_rhou_dxdx_0_0_im1jk_im1p0jk,d2_dif_rhou_dxdx_0_0_im1jk_im1p1jk &
            ,d2_dif_rhou_dxdx_0_0_im1jk &
            ,d2_dif_rhou_dxdx_0_0_ip1jk_ip1m1jk,d2_dif_rhou_dxdx_0_0_ip1jk_ip1p0jk,d2_dif_rhou_dxdx_0_0_ip1jk_ip1p1jk &
            ,d2_dif_rhou_dxdx_0_0_ip1jk &
            ,d2_dif_rhou_dxdy_0_0_im1jk_im1jm1k,d2_dif_rhou_dxdy_0_0_im1jk_im1jp0k,d2_dif_rhou_dxdy_0_0_im1jk_im1jp1k &
            ,d2_dif_rhou_dxdy_0_0_im1jk &
            ,d2_dif_rhou_dxdy_0_0_ip1jk_ip1jm1k,d2_dif_rhou_dxdy_0_0_ip1jk_ip1jp0k,d2_dif_rhou_dxdy_0_0_ip1jk_ip1jp1k &
            ,d2_dif_rhou_dxdy_0_0_ip1jk &
            ,d1_dif_rhou_dx_0_im1jk,d1_dif_rhou_dx_0_ip0jk,d1_dif_rhou_dx_0_ip1jk &
            ,d1_dif_rhou_dx_0_ijk &
            ,d2_dif_rhou_dydx_0_0_ijm1k_im1jm1k,d2_dif_rhou_dydx_0_0_ijm1k_ip0jm1k,d2_dif_rhou_dydx_0_0_ijm1k_ip1jm1k &
            ,d2_dif_rhou_dydx_0_0_ijm1k &
            ,d2_dif_rhou_dydx_0_0_ijp1k_im1jp1k,d2_dif_rhou_dydx_0_0_ijp1k_ip0jp1k,d2_dif_rhou_dydx_0_0_ijp1k_ip1jp1k &
            ,d2_dif_rhou_dydx_0_0_ijp1k &
            ,d2_dif_rhou_dydy_0_0_ijm1k_ijm1m1k,d2_dif_rhou_dydy_0_0_ijm1k_ijm1p0k,d2_dif_rhou_dydy_0_0_ijm1k_ijm1p1k &
            ,d2_dif_rhou_dydy_0_0_ijm1k &
            ,d2_dif_rhou_dydy_0_0_ijp1k_ijp1m1k,d2_dif_rhou_dydy_0_0_ijp1k_ijp1p0k,d2_dif_rhou_dydy_0_0_ijp1k_ijp1p1k &
            ,d2_dif_rhou_dydy_0_0_ijp1k &
            ,d1_dif_rhou_dy_0_ijm1k,d1_dif_rhou_dy_0_ijp0k,d1_dif_rhou_dy_0_ijp1k &
            ,d1_dif_rhou_dy_0_ijk &
            ,d2_dif_rhov_dxdx_0_0_im1jk_im1m1jk,d2_dif_rhov_dxdx_0_0_im1jk_im1p0jk,d2_dif_rhov_dxdx_0_0_im1jk_im1p1jk &
            ,d2_dif_rhov_dxdx_0_0_im1jk &
            ,d2_dif_rhov_dxdx_0_0_ip1jk_ip1m1jk,d2_dif_rhov_dxdx_0_0_ip1jk_ip1p0jk,d2_dif_rhov_dxdx_0_0_ip1jk_ip1p1jk &
            ,d2_dif_rhov_dxdx_0_0_ip1jk &
            ,d2_dif_rhov_dxdy_0_0_im1jk_im1jm1k,d2_dif_rhov_dxdy_0_0_im1jk_im1jp0k,d2_dif_rhov_dxdy_0_0_im1jk_im1jp1k &
            ,d2_dif_rhov_dxdy_0_0_im1jk &
            ,d2_dif_rhov_dxdy_0_0_ip1jk_ip1jm1k,d2_dif_rhov_dxdy_0_0_ip1jk_ip1jp0k,d2_dif_rhov_dxdy_0_0_ip1jk_ip1jp1k &
            ,d2_dif_rhov_dxdy_0_0_ip1jk &
            ,d1_dif_rhov_dx_0_im1jk,d1_dif_rhov_dx_0_ip0jk,d1_dif_rhov_dx_0_ip1jk &
            ,d1_dif_rhov_dx_0_ijk &
            ,d2_dif_rhov_dydx_0_0_ijm1k_im1jm1k,d2_dif_rhov_dydx_0_0_ijm1k_ip0jm1k,d2_dif_rhov_dydx_0_0_ijm1k_ip1jm1k &
            ,d2_dif_rhov_dydx_0_0_ijm1k &
            ,d2_dif_rhov_dydx_0_0_ijp1k_im1jp1k,d2_dif_rhov_dydx_0_0_ijp1k_ip0jp1k,d2_dif_rhov_dydx_0_0_ijp1k_ip1jp1k &
            ,d2_dif_rhov_dydx_0_0_ijp1k &
            ,d2_dif_rhov_dydy_0_0_ijm1k_ijm1m1k,d2_dif_rhov_dydy_0_0_ijm1k_ijm1p0k,d2_dif_rhov_dydy_0_0_ijm1k_ijm1p1k &
            ,d2_dif_rhov_dydy_0_0_ijm1k &
            ,d2_dif_rhov_dydy_0_0_ijp1k_ijp1m1k,d2_dif_rhov_dydy_0_0_ijp1k_ijp1p0k,d2_dif_rhov_dydy_0_0_ijp1k_ijp1p1k &
            ,d2_dif_rhov_dydy_0_0_ijp1k &
            ,d1_dif_rhov_dy_0_ijm1k,d1_dif_rhov_dy_0_ijp0k,d1_dif_rhov_dy_0_ijp1k &
            ,d1_dif_rhov_dy_0_ijk &
            ,d2_dif_et_dxdx_0_0_im1jk_im1m1jk,d2_dif_et_dxdx_0_0_im1jk_im1p0jk,d2_dif_et_dxdx_0_0_im1jk_im1p1jk &
            ,d2_dif_et_dxdx_0_0_im1jk &
            ,d2_dif_et_dxdx_0_0_ip1jk_ip1m1jk,d2_dif_et_dxdx_0_0_ip1jk_ip1p0jk,d2_dif_et_dxdx_0_0_ip1jk_ip1p1jk &
            ,d2_dif_et_dxdx_0_0_ip1jk &
            ,d1_dif_et_dx_0_im1jk,d1_dif_et_dx_0_ip0jk,d1_dif_et_dx_0_ip1jk &
            ,d1_dif_et_dx_0_ijk &
            ,d2_dif_et_dydy_0_0_ijm1k_ijm1m1k,d2_dif_et_dydy_0_0_ijm1k_ijm1p0k,d2_dif_et_dydy_0_0_ijm1k_ijm1p1k &
            ,d2_dif_et_dydy_0_0_ijm1k &
            ,d2_dif_et_dydy_0_0_ijp1k_ijp1m1k,d2_dif_et_dydy_0_0_ijp1k_ijp1p0k,d2_dif_et_dydy_0_0_ijp1k_ijp1p1k &
            ,d2_dif_et_dydy_0_0_ijp1k &
            ,d1_dif_et_dy_0_ijm1k,d1_dif_et_dy_0_ijp0k,d1_dif_et_dy_0_ijp1k &
            ,d1_dif_et_dy_0_ijk &
            ,d1_dif_nut_dx_0_im1jk,d1_dif_nut_dx_0_ip0jk,d1_dif_nut_dx_0_ip1jk &
            ,d1_dif_nut_dx_0_ijk &
            ,d1_dif_nut_dx_1_im1jk,d1_dif_nut_dx_1_ip0jk,d1_dif_nut_dx_1_ip1jk &
            ,d1_dif_nut_dx_1_ijk &
            ,d1_dif_nut_dy_0_ijm1k,d1_dif_nut_dy_0_ijp0k,d1_dif_nut_dy_0_ijp1k &
            ,d1_dif_nut_dy_0_ijk &
            ,d1_dif_nut_dy_1_ijm1k,d1_dif_nut_dy_1_ijp0k,d1_dif_nut_dy_1_ijp1k &
            ,d1_dif_nut_dy_1_ijk 

integer :: i,j,k,nq
integer :: bi,bj,bk
integer :: biend,bjend,bkend
integer :: size_bi,size_bj,size_bk

integer :: indvars(neq),indvarsst(neqst)


  indvars   = ind(1:neq)
  indvarsst = ind(1+neq:neq+neqst) 



!***********************************************************
!                                                           
! Start building RHS with source terms (2D) ****************
!                                                           
!***********************************************************


     do j=idloop(3),idloop(4) 
 
      do i=idloop(1),idloop(2) 


!***********************************************************
!                                                           
! building source terms in RHS for d(rho)/dt ***************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([rho*u]_1x)+deltayI*([rho*v]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_conv_rho_dx_0_im2jk = q(i-2,j,indvars(1))*q(i-2,j,indvars(2))

d1_conv_rho_dx_0_im1jk = q(i-1,j,indvars(1))*q(i-1,j,indvars(2))

d1_conv_rho_dx_0_ip1jk = q(i+1,j,indvars(1))*q(i+1,j,indvars(2))

d1_conv_rho_dx_0_ip2jk = q(i+2,j,indvars(1))*q(i+2,j,indvars(2))

d1_conv_rho_dx_0_ijk = 0.08333333333333333_wp*d1_conv_rho_dx_0_im2jk-&
          0.666666666666667_wp*d1_conv_rho_dx_0_im1jk+&
          0.666666666666667_wp*d1_conv_rho_dx_0_ip1jk-&
          0.08333333333333333_wp*d1_conv_rho_dx_0_ip2jk

d1_conv_rho_dx_0_ijk = d1_conv_rho_dx_0_ijk*param_float(1)

d1_conv_rho_dy_0_ijm2k = q(i,j-2,indvars(1))*q(i,j-2,indvars(3))

d1_conv_rho_dy_0_ijm1k = q(i,j-1,indvars(1))*q(i,j-1,indvars(3))

d1_conv_rho_dy_0_ijp1k = q(i,j+1,indvars(1))*q(i,j+1,indvars(3))

d1_conv_rho_dy_0_ijp2k = q(i,j+2,indvars(1))*q(i,j+2,indvars(3))

d1_conv_rho_dy_0_ijk = 0.08333333333333333_wp*d1_conv_rho_dy_0_ijm2k-&
          0.666666666666667_wp*d1_conv_rho_dy_0_ijm1k+&
          0.666666666666667_wp*d1_conv_rho_dy_0_ijp1k-&
          0.08333333333333333_wp*d1_conv_rho_dy_0_ijp2k

d1_conv_rho_dy_0_ijk = d1_conv_rho_dy_0_ijk*param_float(2)



!***********************************************************
!                                                           
! Update RHS terms for d(rho)/dt ***************************
!                                                           
!***********************************************************


rhs(i,j,indvars(1)) =   -  ( qst(i,j,indvarsst(10))*(d1_conv_rho_dx_0_ijk)+&
                    qst(i,j,indvarsst(11))*(d1_conv_rho_dy_0_ijk) ) 



!***********************************************************
!                                                           
! building source terms in RHS for d(rho u)/dt *************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([rho*u*u+p]_1x)+deltayI*([rho*v*u]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_conv_rhou_dx_0_im2jk = q(i-2,j,indvars(1))*q(i-2,j,indvars(2))*q(i-2,j,indvars(2))+(param_float(3 + 5))*q(i-2,j,indvars(1))*((q(i-2,j,indvars(4))-&
                    0.5_wp*(q(i-2,j,indvars(2))*q(i-2,j,indvars(2))+&
                    q(i-2,j,indvars(3))*q(i-2,j,indvars(3)))))

d1_conv_rhou_dx_0_im1jk = q(i-1,j,indvars(1))*q(i-1,j,indvars(2))*q(i-1,j,indvars(2))+(param_float(3 + 5))*q(i-1,j,indvars(1))*((q(i-1,j,indvars(4))-&
                    0.5_wp*(q(i-1,j,indvars(2))*q(i-1,j,indvars(2))+&
                    q(i-1,j,indvars(3))*q(i-1,j,indvars(3)))))

d1_conv_rhou_dx_0_ip1jk = q(i+1,j,indvars(1))*q(i+1,j,indvars(2))*q(i+1,j,indvars(2))+(param_float(3 + 5))*q(i+1,j,indvars(1))*((q(i+1,j,indvars(4))-&
                    0.5_wp*(q(i+1,j,indvars(2))*q(i+1,j,indvars(2))+&
                    q(i+1,j,indvars(3))*q(i+1,j,indvars(3)))))

d1_conv_rhou_dx_0_ip2jk = q(i+2,j,indvars(1))*q(i+2,j,indvars(2))*q(i+2,j,indvars(2))+(param_float(3 + 5))*q(i+2,j,indvars(1))*((q(i+2,j,indvars(4))-&
                    0.5_wp*(q(i+2,j,indvars(2))*q(i+2,j,indvars(2))+&
                    q(i+2,j,indvars(3))*q(i+2,j,indvars(3)))))

d1_conv_rhou_dx_0_ijk = 0.08333333333333333_wp*d1_conv_rhou_dx_0_im2jk-&
          0.666666666666667_wp*d1_conv_rhou_dx_0_im1jk+&
          0.666666666666667_wp*d1_conv_rhou_dx_0_ip1jk-&
          0.08333333333333333_wp*d1_conv_rhou_dx_0_ip2jk

d1_conv_rhou_dx_0_ijk = d1_conv_rhou_dx_0_ijk*param_float(1)

d1_conv_rhou_dy_0_ijm2k = q(i,j-2,indvars(1))*q(i,j-2,indvars(3))*q(i,j-2,indvars(2))

d1_conv_rhou_dy_0_ijm1k = q(i,j-1,indvars(1))*q(i,j-1,indvars(3))*q(i,j-1,indvars(2))

d1_conv_rhou_dy_0_ijp1k = q(i,j+1,indvars(1))*q(i,j+1,indvars(3))*q(i,j+1,indvars(2))

d1_conv_rhou_dy_0_ijp2k = q(i,j+2,indvars(1))*q(i,j+2,indvars(3))*q(i,j+2,indvars(2))

d1_conv_rhou_dy_0_ijk = 0.08333333333333333_wp*d1_conv_rhou_dy_0_ijm2k-&
          0.666666666666667_wp*d1_conv_rhou_dy_0_ijm1k+&
          0.666666666666667_wp*d1_conv_rhou_dy_0_ijp1k-&
          0.08333333333333333_wp*d1_conv_rhou_dy_0_ijp2k

d1_conv_rhou_dy_0_ijk = d1_conv_rhou_dy_0_ijk*param_float(2)



!***********************************************************
!                                                           
! Update RHS terms for d(rho u)/dt *************************
!                                                           
!***********************************************************


rhs(i,j,indvars(2)) =   -  ( qst(i,j,indvarsst(10))*(d1_conv_rhou_dx_0_ijk)+&
                    qst(i,j,indvarsst(11))*(d1_conv_rhou_dy_0_ijk) ) 



!***********************************************************
!                                                           
! building source terms in RHS for d(rho v)/dt *************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([rho*u*v]_1x)+deltayI*([rho*v*v+p]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_conv_rhov_dx_0_im2jk = q(i-2,j,indvars(1))*q(i-2,j,indvars(2))*q(i-2,j,indvars(3))

d1_conv_rhov_dx_0_im1jk = q(i-1,j,indvars(1))*q(i-1,j,indvars(2))*q(i-1,j,indvars(3))

d1_conv_rhov_dx_0_ip1jk = q(i+1,j,indvars(1))*q(i+1,j,indvars(2))*q(i+1,j,indvars(3))

d1_conv_rhov_dx_0_ip2jk = q(i+2,j,indvars(1))*q(i+2,j,indvars(2))*q(i+2,j,indvars(3))

d1_conv_rhov_dx_0_ijk = 0.08333333333333333_wp*d1_conv_rhov_dx_0_im2jk-&
          0.666666666666667_wp*d1_conv_rhov_dx_0_im1jk+&
          0.666666666666667_wp*d1_conv_rhov_dx_0_ip1jk-&
          0.08333333333333333_wp*d1_conv_rhov_dx_0_ip2jk

d1_conv_rhov_dx_0_ijk = d1_conv_rhov_dx_0_ijk*param_float(1)

d1_conv_rhov_dy_0_ijm2k = q(i,j-2,indvars(1))*q(i,j-2,indvars(3))*q(i,j-2,indvars(3))+(param_float(3 + 5))*q(i,j-2,indvars(1))*((q(i,j-2,indvars(4))-&
                    0.5_wp*(q(i,j-2,indvars(2))*q(i,j-2,indvars(2))+&
                    q(i,j-2,indvars(3))*q(i,j-2,indvars(3)))))

d1_conv_rhov_dy_0_ijm1k = q(i,j-1,indvars(1))*q(i,j-1,indvars(3))*q(i,j-1,indvars(3))+(param_float(3 + 5))*q(i,j-1,indvars(1))*((q(i,j-1,indvars(4))-&
                    0.5_wp*(q(i,j-1,indvars(2))*q(i,j-1,indvars(2))+&
                    q(i,j-1,indvars(3))*q(i,j-1,indvars(3)))))

d1_conv_rhov_dy_0_ijp1k = q(i,j+1,indvars(1))*q(i,j+1,indvars(3))*q(i,j+1,indvars(3))+(param_float(3 + 5))*q(i,j+1,indvars(1))*((q(i,j+1,indvars(4))-&
                    0.5_wp*(q(i,j+1,indvars(2))*q(i,j+1,indvars(2))+&
                    q(i,j+1,indvars(3))*q(i,j+1,indvars(3)))))

d1_conv_rhov_dy_0_ijp2k = q(i,j+2,indvars(1))*q(i,j+2,indvars(3))*q(i,j+2,indvars(3))+(param_float(3 + 5))*q(i,j+2,indvars(1))*((q(i,j+2,indvars(4))-&
                    0.5_wp*(q(i,j+2,indvars(2))*q(i,j+2,indvars(2))+&
                    q(i,j+2,indvars(3))*q(i,j+2,indvars(3)))))

d1_conv_rhov_dy_0_ijk = 0.08333333333333333_wp*d1_conv_rhov_dy_0_ijm2k-&
          0.666666666666667_wp*d1_conv_rhov_dy_0_ijm1k+&
          0.666666666666667_wp*d1_conv_rhov_dy_0_ijp1k-&
          0.08333333333333333_wp*d1_conv_rhov_dy_0_ijp2k

d1_conv_rhov_dy_0_ijk = d1_conv_rhov_dy_0_ijk*param_float(2)



!***********************************************************
!                                                           
! Update RHS terms for d(rho v)/dt *************************
!                                                           
!***********************************************************


rhs(i,j,indvars(3)) =   -  ( qst(i,j,indvarsst(10))*(d1_conv_rhov_dx_0_ijk)+&
                    qst(i,j,indvarsst(11))*(d1_conv_rhov_dy_0_ijk) ) 



!***********************************************************
!                                                           
! building source terms in RHS for d(rho et)/dt ************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([(rho*et+p)*u]_1x)+deltayI*([(rho*et+p)*v]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_conv_et_dx_0_im2jk = (q(i-2,j,indvars(1))*q(i-2,j,indvars(4))+&
                    (param_float(3 + 5))*q(i-2,j,indvars(1))*((q(i-2,j,indvars(4))-&
                    0.5_wp*(q(i-2,j,indvars(2))*q(i-2,j,indvars(2))+&
                    q(i-2,j,indvars(3))*q(i-2,j,indvars(3))))))*q(i-2,j,indvars(2))

d1_conv_et_dx_0_im1jk = (q(i-1,j,indvars(1))*q(i-1,j,indvars(4))+&
                    (param_float(3 + 5))*q(i-1,j,indvars(1))*((q(i-1,j,indvars(4))-&
                    0.5_wp*(q(i-1,j,indvars(2))*q(i-1,j,indvars(2))+&
                    q(i-1,j,indvars(3))*q(i-1,j,indvars(3))))))*q(i-1,j,indvars(2))

d1_conv_et_dx_0_ip1jk = (q(i+1,j,indvars(1))*q(i+1,j,indvars(4))+&
                    (param_float(3 + 5))*q(i+1,j,indvars(1))*((q(i+1,j,indvars(4))-&
                    0.5_wp*(q(i+1,j,indvars(2))*q(i+1,j,indvars(2))+&
                    q(i+1,j,indvars(3))*q(i+1,j,indvars(3))))))*q(i+1,j,indvars(2))

d1_conv_et_dx_0_ip2jk = (q(i+2,j,indvars(1))*q(i+2,j,indvars(4))+&
                    (param_float(3 + 5))*q(i+2,j,indvars(1))*((q(i+2,j,indvars(4))-&
                    0.5_wp*(q(i+2,j,indvars(2))*q(i+2,j,indvars(2))+&
                    q(i+2,j,indvars(3))*q(i+2,j,indvars(3))))))*q(i+2,j,indvars(2))

d1_conv_et_dx_0_ijk = 0.08333333333333333_wp*d1_conv_et_dx_0_im2jk-&
          0.666666666666667_wp*d1_conv_et_dx_0_im1jk+&
          0.666666666666667_wp*d1_conv_et_dx_0_ip1jk-&
          0.08333333333333333_wp*d1_conv_et_dx_0_ip2jk

d1_conv_et_dx_0_ijk = d1_conv_et_dx_0_ijk*param_float(1)

d1_conv_et_dy_0_ijm2k = (q(i,j-2,indvars(1))*q(i,j-2,indvars(4))+&
                    (param_float(3 + 5))*q(i,j-2,indvars(1))*((q(i,j-2,indvars(4))-&
                    0.5_wp*(q(i,j-2,indvars(2))*q(i,j-2,indvars(2))+&
                    q(i,j-2,indvars(3))*q(i,j-2,indvars(3))))))*q(i,j-2,indvars(3))

d1_conv_et_dy_0_ijm1k = (q(i,j-1,indvars(1))*q(i,j-1,indvars(4))+&
                    (param_float(3 + 5))*q(i,j-1,indvars(1))*((q(i,j-1,indvars(4))-&
                    0.5_wp*(q(i,j-1,indvars(2))*q(i,j-1,indvars(2))+&
                    q(i,j-1,indvars(3))*q(i,j-1,indvars(3))))))*q(i,j-1,indvars(3))

d1_conv_et_dy_0_ijp1k = (q(i,j+1,indvars(1))*q(i,j+1,indvars(4))+&
                    (param_float(3 + 5))*q(i,j+1,indvars(1))*((q(i,j+1,indvars(4))-&
                    0.5_wp*(q(i,j+1,indvars(2))*q(i,j+1,indvars(2))+&
                    q(i,j+1,indvars(3))*q(i,j+1,indvars(3))))))*q(i,j+1,indvars(3))

d1_conv_et_dy_0_ijp2k = (q(i,j+2,indvars(1))*q(i,j+2,indvars(4))+&
                    (param_float(3 + 5))*q(i,j+2,indvars(1))*((q(i,j+2,indvars(4))-&
                    0.5_wp*(q(i,j+2,indvars(2))*q(i,j+2,indvars(2))+&
                    q(i,j+2,indvars(3))*q(i,j+2,indvars(3))))))*q(i,j+2,indvars(3))

d1_conv_et_dy_0_ijk = 0.08333333333333333_wp*d1_conv_et_dy_0_ijm2k-&
          0.666666666666667_wp*d1_conv_et_dy_0_ijm1k+&
          0.666666666666667_wp*d1_conv_et_dy_0_ijp1k-&
          0.08333333333333333_wp*d1_conv_et_dy_0_ijp2k

d1_conv_et_dy_0_ijk = d1_conv_et_dy_0_ijk*param_float(2)



!***********************************************************
!                                                           
! Update RHS terms for d(rho et)/dt ************************
!                                                           
!***********************************************************


rhs(i,j,indvars(4)) =   -  ( qst(i,j,indvarsst(10))*(d1_conv_et_dx_0_ijk)+&
                    qst(i,j,indvarsst(11))*(d1_conv_et_dy_0_ijk) ) 



!***********************************************************
!                                                           
! building source terms in RHS for d(rho nut)/dt ***********
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([rho*u*nut-visc_t*sigmaI*deltaxI*({nut}_1x)]_1x)+deltayI*([rho*v*nut-visc_t*sigmaI*deltayI*({nut}_1y)]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_conv_nut_dxdx_0_0_im2jk_im2m2jk = q(i-2-2,j,indvars(5))

d2_conv_nut_dxdx_0_0_im2jk_im2m1jk = q(i-2-1,j,indvars(5))

d2_conv_nut_dxdx_0_0_im2jk_im2p1jk = q(i-2+1,j,indvars(5))

d2_conv_nut_dxdx_0_0_im2jk_im2p2jk = q(i-2+2,j,indvars(5))

d2_conv_nut_dxdx_0_0_im2jk = 0.08333333333333333_wp*d2_conv_nut_dxdx_0_0_im2jk_im2m2jk-&
          0.666666666666667_wp*d2_conv_nut_dxdx_0_0_im2jk_im2m1jk+&
          0.666666666666667_wp*d2_conv_nut_dxdx_0_0_im2jk_im2p1jk-&
          0.08333333333333333_wp*d2_conv_nut_dxdx_0_0_im2jk_im2p2jk

d2_conv_nut_dxdx_0_0_im2jk = d2_conv_nut_dxdx_0_0_im2jk*param_float(1)

d2_conv_nut_dxdx_0_0_im1jk_im1m2jk = q(i-1-2,j,indvars(5))

d2_conv_nut_dxdx_0_0_im1jk_im1m1jk = q(i-1-1,j,indvars(5))

d2_conv_nut_dxdx_0_0_im1jk_im1p1jk = q(i-1+1,j,indvars(5))

d2_conv_nut_dxdx_0_0_im1jk_im1p2jk = q(i-1+2,j,indvars(5))

d2_conv_nut_dxdx_0_0_im1jk = 0.08333333333333333_wp*d2_conv_nut_dxdx_0_0_im1jk_im1m2jk-&
          0.666666666666667_wp*d2_conv_nut_dxdx_0_0_im1jk_im1m1jk+&
          0.666666666666667_wp*d2_conv_nut_dxdx_0_0_im1jk_im1p1jk-&
          0.08333333333333333_wp*d2_conv_nut_dxdx_0_0_im1jk_im1p2jk

d2_conv_nut_dxdx_0_0_im1jk = d2_conv_nut_dxdx_0_0_im1jk*param_float(1)

d2_conv_nut_dxdx_0_0_ip1jk_ip1m2jk = q(i+1-2,j,indvars(5))

d2_conv_nut_dxdx_0_0_ip1jk_ip1m1jk = q(i+1-1,j,indvars(5))

d2_conv_nut_dxdx_0_0_ip1jk_ip1p1jk = q(i+1+1,j,indvars(5))

d2_conv_nut_dxdx_0_0_ip1jk_ip1p2jk = q(i+1+2,j,indvars(5))

d2_conv_nut_dxdx_0_0_ip1jk = 0.08333333333333333_wp*d2_conv_nut_dxdx_0_0_ip1jk_ip1m2jk-&
          0.666666666666667_wp*d2_conv_nut_dxdx_0_0_ip1jk_ip1m1jk+&
          0.666666666666667_wp*d2_conv_nut_dxdx_0_0_ip1jk_ip1p1jk-&
          0.08333333333333333_wp*d2_conv_nut_dxdx_0_0_ip1jk_ip1p2jk

d2_conv_nut_dxdx_0_0_ip1jk = d2_conv_nut_dxdx_0_0_ip1jk*param_float(1)

d2_conv_nut_dxdx_0_0_ip2jk_ip2m2jk = q(i+2-2,j,indvars(5))

d2_conv_nut_dxdx_0_0_ip2jk_ip2m1jk = q(i+2-1,j,indvars(5))

d2_conv_nut_dxdx_0_0_ip2jk_ip2p1jk = q(i+2+1,j,indvars(5))

d2_conv_nut_dxdx_0_0_ip2jk_ip2p2jk = q(i+2+2,j,indvars(5))

d2_conv_nut_dxdx_0_0_ip2jk = 0.08333333333333333_wp*d2_conv_nut_dxdx_0_0_ip2jk_ip2m2jk-&
          0.666666666666667_wp*d2_conv_nut_dxdx_0_0_ip2jk_ip2m1jk+&
          0.666666666666667_wp*d2_conv_nut_dxdx_0_0_ip2jk_ip2p1jk-&
          0.08333333333333333_wp*d2_conv_nut_dxdx_0_0_ip2jk_ip2p2jk

d2_conv_nut_dxdx_0_0_ip2jk = d2_conv_nut_dxdx_0_0_ip2jk*param_float(1)

d1_conv_nut_dx_0_im2jk = q(i-2,j,indvars(1))*q(i-2,j,indvars(2))*q(i-2,j,indvars(5))-&
                    (1.0_wp)*(1.0_wp+&
                    ((q(i-2,j,indvars(5))/1.0_wp*q(i-2,j,indvars(1)))**3.0_wp/((q(i-2,j,indvars(5))/1.0_wp*q(i-2,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i-2,j,indvars(5))/1.0_wp*q(i-2,j,indvars(1))))*param_float(1 + 5)*param_float(18 + 5)*qst(i-2,j,indvarsst(10))*(d2_conv_nut_dxdx_0_0_im2jk)

d1_conv_nut_dx_0_im1jk = q(i-1,j,indvars(1))*q(i-1,j,indvars(2))*q(i-1,j,indvars(5))-&
                    (1.0_wp)*(1.0_wp+&
                    ((q(i-1,j,indvars(5))/1.0_wp*q(i-1,j,indvars(1)))**3.0_wp/((q(i-1,j,indvars(5))/1.0_wp*q(i-1,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i-1,j,indvars(5))/1.0_wp*q(i-1,j,indvars(1))))*param_float(1 + 5)*param_float(18 + 5)*qst(i-1,j,indvarsst(10))*(d2_conv_nut_dxdx_0_0_im1jk)

d1_conv_nut_dx_0_ip1jk = q(i+1,j,indvars(1))*q(i+1,j,indvars(2))*q(i+1,j,indvars(5))-&
                    (1.0_wp)*(1.0_wp+&
                    ((q(i+1,j,indvars(5))/1.0_wp*q(i+1,j,indvars(1)))**3.0_wp/((q(i+1,j,indvars(5))/1.0_wp*q(i+1,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i+1,j,indvars(5))/1.0_wp*q(i+1,j,indvars(1))))*param_float(1 + 5)*param_float(18 + 5)*qst(i+1,j,indvarsst(10))*(d2_conv_nut_dxdx_0_0_ip1jk)

d1_conv_nut_dx_0_ip2jk = q(i+2,j,indvars(1))*q(i+2,j,indvars(2))*q(i+2,j,indvars(5))-&
                    (1.0_wp)*(1.0_wp+&
                    ((q(i+2,j,indvars(5))/1.0_wp*q(i+2,j,indvars(1)))**3.0_wp/((q(i+2,j,indvars(5))/1.0_wp*q(i+2,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i+2,j,indvars(5))/1.0_wp*q(i+2,j,indvars(1))))*param_float(1 + 5)*param_float(18 + 5)*qst(i+2,j,indvarsst(10))*(d2_conv_nut_dxdx_0_0_ip2jk)

d1_conv_nut_dx_0_ijk = 0.08333333333333333_wp*d1_conv_nut_dx_0_im2jk-&
          0.666666666666667_wp*d1_conv_nut_dx_0_im1jk+&
          0.666666666666667_wp*d1_conv_nut_dx_0_ip1jk-&
          0.08333333333333333_wp*d1_conv_nut_dx_0_ip2jk

d1_conv_nut_dx_0_ijk = d1_conv_nut_dx_0_ijk*param_float(1)

d2_conv_nut_dydy_0_0_ijm2k_ijm2m2k = q(i,j-2-2,indvars(5))

d2_conv_nut_dydy_0_0_ijm2k_ijm2m1k = q(i,j-2-1,indvars(5))

d2_conv_nut_dydy_0_0_ijm2k_ijm2p1k = q(i,j-2+1,indvars(5))

d2_conv_nut_dydy_0_0_ijm2k_ijm2p2k = q(i,j-2+2,indvars(5))

d2_conv_nut_dydy_0_0_ijm2k = 0.08333333333333333_wp*d2_conv_nut_dydy_0_0_ijm2k_ijm2m2k-&
          0.666666666666667_wp*d2_conv_nut_dydy_0_0_ijm2k_ijm2m1k+&
          0.666666666666667_wp*d2_conv_nut_dydy_0_0_ijm2k_ijm2p1k-&
          0.08333333333333333_wp*d2_conv_nut_dydy_0_0_ijm2k_ijm2p2k

d2_conv_nut_dydy_0_0_ijm2k = d2_conv_nut_dydy_0_0_ijm2k*param_float(2)

d2_conv_nut_dydy_0_0_ijm1k_ijm1m2k = q(i,j-1-2,indvars(5))

d2_conv_nut_dydy_0_0_ijm1k_ijm1m1k = q(i,j-1-1,indvars(5))

d2_conv_nut_dydy_0_0_ijm1k_ijm1p1k = q(i,j-1+1,indvars(5))

d2_conv_nut_dydy_0_0_ijm1k_ijm1p2k = q(i,j-1+2,indvars(5))

d2_conv_nut_dydy_0_0_ijm1k = 0.08333333333333333_wp*d2_conv_nut_dydy_0_0_ijm1k_ijm1m2k-&
          0.666666666666667_wp*d2_conv_nut_dydy_0_0_ijm1k_ijm1m1k+&
          0.666666666666667_wp*d2_conv_nut_dydy_0_0_ijm1k_ijm1p1k-&
          0.08333333333333333_wp*d2_conv_nut_dydy_0_0_ijm1k_ijm1p2k

d2_conv_nut_dydy_0_0_ijm1k = d2_conv_nut_dydy_0_0_ijm1k*param_float(2)

d2_conv_nut_dydy_0_0_ijp1k_ijp1m2k = q(i,j+1-2,indvars(5))

d2_conv_nut_dydy_0_0_ijp1k_ijp1m1k = q(i,j+1-1,indvars(5))

d2_conv_nut_dydy_0_0_ijp1k_ijp1p1k = q(i,j+1+1,indvars(5))

d2_conv_nut_dydy_0_0_ijp1k_ijp1p2k = q(i,j+1+2,indvars(5))

d2_conv_nut_dydy_0_0_ijp1k = 0.08333333333333333_wp*d2_conv_nut_dydy_0_0_ijp1k_ijp1m2k-&
          0.666666666666667_wp*d2_conv_nut_dydy_0_0_ijp1k_ijp1m1k+&
          0.666666666666667_wp*d2_conv_nut_dydy_0_0_ijp1k_ijp1p1k-&
          0.08333333333333333_wp*d2_conv_nut_dydy_0_0_ijp1k_ijp1p2k

d2_conv_nut_dydy_0_0_ijp1k = d2_conv_nut_dydy_0_0_ijp1k*param_float(2)

d2_conv_nut_dydy_0_0_ijp2k_ijp2m2k = q(i,j+2-2,indvars(5))

d2_conv_nut_dydy_0_0_ijp2k_ijp2m1k = q(i,j+2-1,indvars(5))

d2_conv_nut_dydy_0_0_ijp2k_ijp2p1k = q(i,j+2+1,indvars(5))

d2_conv_nut_dydy_0_0_ijp2k_ijp2p2k = q(i,j+2+2,indvars(5))

d2_conv_nut_dydy_0_0_ijp2k = 0.08333333333333333_wp*d2_conv_nut_dydy_0_0_ijp2k_ijp2m2k-&
          0.666666666666667_wp*d2_conv_nut_dydy_0_0_ijp2k_ijp2m1k+&
          0.666666666666667_wp*d2_conv_nut_dydy_0_0_ijp2k_ijp2p1k-&
          0.08333333333333333_wp*d2_conv_nut_dydy_0_0_ijp2k_ijp2p2k

d2_conv_nut_dydy_0_0_ijp2k = d2_conv_nut_dydy_0_0_ijp2k*param_float(2)

d1_conv_nut_dy_0_ijm2k = q(i,j-2,indvars(1))*q(i,j-2,indvars(3))*q(i,j-2,indvars(5))-&
                    (1.0_wp)*(1.0_wp+&
                    ((q(i,j-2,indvars(5))/1.0_wp*q(i,j-2,indvars(1)))**3.0_wp/((q(i,j-2,indvars(5))/1.0_wp*q(i,j-2,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i,j-2,indvars(5))/1.0_wp*q(i,j-2,indvars(1))))*param_float(1 + 5)*param_float(18 + 5)*qst(i,j-2,indvarsst(11))*(d2_conv_nut_dydy_0_0_ijm2k)

d1_conv_nut_dy_0_ijm1k = q(i,j-1,indvars(1))*q(i,j-1,indvars(3))*q(i,j-1,indvars(5))-&
                    (1.0_wp)*(1.0_wp+&
                    ((q(i,j-1,indvars(5))/1.0_wp*q(i,j-1,indvars(1)))**3.0_wp/((q(i,j-1,indvars(5))/1.0_wp*q(i,j-1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i,j-1,indvars(5))/1.0_wp*q(i,j-1,indvars(1))))*param_float(1 + 5)*param_float(18 + 5)*qst(i,j-1,indvarsst(11))*(d2_conv_nut_dydy_0_0_ijm1k)

d1_conv_nut_dy_0_ijp1k = q(i,j+1,indvars(1))*q(i,j+1,indvars(3))*q(i,j+1,indvars(5))-&
                    (1.0_wp)*(1.0_wp+&
                    ((q(i,j+1,indvars(5))/1.0_wp*q(i,j+1,indvars(1)))**3.0_wp/((q(i,j+1,indvars(5))/1.0_wp*q(i,j+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i,j+1,indvars(5))/1.0_wp*q(i,j+1,indvars(1))))*param_float(1 + 5)*param_float(18 + 5)*qst(i,j+1,indvarsst(11))*(d2_conv_nut_dydy_0_0_ijp1k)

d1_conv_nut_dy_0_ijp2k = q(i,j+2,indvars(1))*q(i,j+2,indvars(3))*q(i,j+2,indvars(5))-&
                    (1.0_wp)*(1.0_wp+&
                    ((q(i,j+2,indvars(5))/1.0_wp*q(i,j+2,indvars(1)))**3.0_wp/((q(i,j+2,indvars(5))/1.0_wp*q(i,j+2,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i,j+2,indvars(5))/1.0_wp*q(i,j+2,indvars(1))))*param_float(1 + 5)*param_float(18 + 5)*qst(i,j+2,indvarsst(11))*(d2_conv_nut_dydy_0_0_ijp2k)

d1_conv_nut_dy_0_ijk = 0.08333333333333333_wp*d1_conv_nut_dy_0_ijm2k-&
          0.666666666666667_wp*d1_conv_nut_dy_0_ijm1k+&
          0.666666666666667_wp*d1_conv_nut_dy_0_ijp1k-&
          0.08333333333333333_wp*d1_conv_nut_dy_0_ijp2k

d1_conv_nut_dy_0_ijk = d1_conv_nut_dy_0_ijk*param_float(2)



!***********************************************************
!                                                           
! Update RHS terms for d(rho nut)/dt ***********************
!                                                           
!***********************************************************


rhs(i,j,indvars(5)) =   -  ( qst(i,j,indvarsst(10))*(d1_conv_nut_dx_0_ijk)+&
                    qst(i,j,indvarsst(11))*(d1_conv_nut_dy_0_ijk) ) 

     enddo
   enddo


!***********************************************************
!                                                           
! Start building RHS with source terms (2D) ****************
!                                                           
!***********************************************************


     do j=idloop(3),idloop(4) 
 
      do i=idloop(1),idloop(2) 


!***********************************************************
!                                                           
! building source terms in RHS for d(rho u)/dt *************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([-visc_t*(2.0_wp*deltaxI*({u}_1x)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y)))]_1x)+deltayI*([-visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_dif_rhou_dxdx_0_0_im1jk_im1m1jk = q(i-1-1,j,indvars(2))

d2_dif_rhou_dxdx_0_0_im1jk_im1p1jk = q(i-1+1,j,indvars(2))

d2_dif_rhou_dxdx_0_0_im1jk = -&
          0.5_wp*d2_dif_rhou_dxdx_0_0_im1jk_im1m1jk+&
          0.5_wp*d2_dif_rhou_dxdx_0_0_im1jk_im1p1jk

d2_dif_rhou_dxdx_0_0_im1jk = d2_dif_rhou_dxdx_0_0_im1jk*param_float(1)

d2_dif_rhou_dxdx_0_0_ip1jk_ip1m1jk = q(i+1-1,j,indvars(2))

d2_dif_rhou_dxdx_0_0_ip1jk_ip1p1jk = q(i+1+1,j,indvars(2))

d2_dif_rhou_dxdx_0_0_ip1jk = -&
          0.5_wp*d2_dif_rhou_dxdx_0_0_ip1jk_ip1m1jk+&
          0.5_wp*d2_dif_rhou_dxdx_0_0_ip1jk_ip1p1jk

d2_dif_rhou_dxdx_0_0_ip1jk = d2_dif_rhou_dxdx_0_0_ip1jk*param_float(1)

d2_dif_rhou_dxdy_0_0_im1jk_im1jm1k = q(i-1,j-1,indvars(3))

d2_dif_rhou_dxdy_0_0_im1jk_im1jp1k = q(i-1,j+1,indvars(3))

d2_dif_rhou_dxdy_0_0_im1jk = -&
          0.5_wp*d2_dif_rhou_dxdy_0_0_im1jk_im1jm1k+&
          0.5_wp*d2_dif_rhou_dxdy_0_0_im1jk_im1jp1k

d2_dif_rhou_dxdy_0_0_im1jk = d2_dif_rhou_dxdy_0_0_im1jk*param_float(2)

d2_dif_rhou_dxdy_0_0_ip1jk_ip1jm1k = q(i+1,j-1,indvars(3))

d2_dif_rhou_dxdy_0_0_ip1jk_ip1jp1k = q(i+1,j+1,indvars(3))

d2_dif_rhou_dxdy_0_0_ip1jk = -&
          0.5_wp*d2_dif_rhou_dxdy_0_0_ip1jk_ip1jm1k+&
          0.5_wp*d2_dif_rhou_dxdy_0_0_ip1jk_ip1jp1k

d2_dif_rhou_dxdy_0_0_ip1jk = d2_dif_rhou_dxdy_0_0_ip1jk*param_float(2)

d1_dif_rhou_dx_0_im1jk = -(1.0_wp)*(1.0_wp+&
                    ((q(i-1,j,indvars(5))/1.0_wp*q(i-1,j,indvars(1)))**3.0_wp/((q(i-1,j,indvars(5))/1.0_wp*q(i-1,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i-1,j,indvars(5))/1.0_wp*q(i-1,j,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(i-1,j,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_im1jk)-&
                    2.0_wp/3.0_wp*(qst(i-1,j,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_im1jk)+&
                    qst(i-1,j,indvarsst(11))*(d2_dif_rhou_dxdy_0_0_im1jk)))

d1_dif_rhou_dx_0_ip1jk = -(1.0_wp)*(1.0_wp+&
                    ((q(i+1,j,indvars(5))/1.0_wp*q(i+1,j,indvars(1)))**3.0_wp/((q(i+1,j,indvars(5))/1.0_wp*q(i+1,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i+1,j,indvars(5))/1.0_wp*q(i+1,j,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(i+1,j,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_ip1jk)-&
                    2.0_wp/3.0_wp*(qst(i+1,j,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_ip1jk)+&
                    qst(i+1,j,indvarsst(11))*(d2_dif_rhou_dxdy_0_0_ip1jk)))

d1_dif_rhou_dx_0_ijk = -&
          0.5_wp*d1_dif_rhou_dx_0_im1jk+&
          0.5_wp*d1_dif_rhou_dx_0_ip1jk

d1_dif_rhou_dx_0_ijk = d1_dif_rhou_dx_0_ijk*param_float(1)

d2_dif_rhou_dydx_0_0_ijm1k_im1jm1k = q(i-1,j-1,indvars(3))

d2_dif_rhou_dydx_0_0_ijm1k_ip1jm1k = q(i+1,j-1,indvars(3))

d2_dif_rhou_dydx_0_0_ijm1k = -&
          0.5_wp*d2_dif_rhou_dydx_0_0_ijm1k_im1jm1k+&
          0.5_wp*d2_dif_rhou_dydx_0_0_ijm1k_ip1jm1k

d2_dif_rhou_dydx_0_0_ijm1k = d2_dif_rhou_dydx_0_0_ijm1k*param_float(1)

d2_dif_rhou_dydx_0_0_ijp1k_im1jp1k = q(i-1,j+1,indvars(3))

d2_dif_rhou_dydx_0_0_ijp1k_ip1jp1k = q(i+1,j+1,indvars(3))

d2_dif_rhou_dydx_0_0_ijp1k = -&
          0.5_wp*d2_dif_rhou_dydx_0_0_ijp1k_im1jp1k+&
          0.5_wp*d2_dif_rhou_dydx_0_0_ijp1k_ip1jp1k

d2_dif_rhou_dydx_0_0_ijp1k = d2_dif_rhou_dydx_0_0_ijp1k*param_float(1)

d2_dif_rhou_dydy_0_0_ijm1k_ijm1m1k = q(i,j-1-1,indvars(2))

d2_dif_rhou_dydy_0_0_ijm1k_ijm1p1k = q(i,j-1+1,indvars(2))

d2_dif_rhou_dydy_0_0_ijm1k = -&
          0.5_wp*d2_dif_rhou_dydy_0_0_ijm1k_ijm1m1k+&
          0.5_wp*d2_dif_rhou_dydy_0_0_ijm1k_ijm1p1k

d2_dif_rhou_dydy_0_0_ijm1k = d2_dif_rhou_dydy_0_0_ijm1k*param_float(2)

d2_dif_rhou_dydy_0_0_ijp1k_ijp1m1k = q(i,j+1-1,indvars(2))

d2_dif_rhou_dydy_0_0_ijp1k_ijp1p1k = q(i,j+1+1,indvars(2))

d2_dif_rhou_dydy_0_0_ijp1k = -&
          0.5_wp*d2_dif_rhou_dydy_0_0_ijp1k_ijp1m1k+&
          0.5_wp*d2_dif_rhou_dydy_0_0_ijp1k_ijp1p1k

d2_dif_rhou_dydy_0_0_ijp1k = d2_dif_rhou_dydy_0_0_ijp1k*param_float(2)

d1_dif_rhou_dy_0_ijm1k = -(1.0_wp)*(1.0_wp+&
                    ((q(i,j-1,indvars(5))/1.0_wp*q(i,j-1,indvars(1)))**3.0_wp/((q(i,j-1,indvars(5))/1.0_wp*q(i,j-1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i,j-1,indvars(5))/1.0_wp*q(i,j-1,indvars(1))))*param_float(1 + 5)*(qst(i,j-1,indvarsst(11))*(d2_dif_rhou_dydy_0_0_ijm1k)+&
                    qst(i,j-1,indvarsst(10))*(d2_dif_rhou_dydx_0_0_ijm1k))

d1_dif_rhou_dy_0_ijp1k = -(1.0_wp)*(1.0_wp+&
                    ((q(i,j+1,indvars(5))/1.0_wp*q(i,j+1,indvars(1)))**3.0_wp/((q(i,j+1,indvars(5))/1.0_wp*q(i,j+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i,j+1,indvars(5))/1.0_wp*q(i,j+1,indvars(1))))*param_float(1 + 5)*(qst(i,j+1,indvarsst(11))*(d2_dif_rhou_dydy_0_0_ijp1k)+&
                    qst(i,j+1,indvarsst(10))*(d2_dif_rhou_dydx_0_0_ijp1k))

d1_dif_rhou_dy_0_ijk = -&
          0.5_wp*d1_dif_rhou_dy_0_ijm1k+&
          0.5_wp*d1_dif_rhou_dy_0_ijp1k

d1_dif_rhou_dy_0_ijk = d1_dif_rhou_dy_0_ijk*param_float(2)



!***********************************************************
!                                                           
! Update RHS terms for d(rho u)/dt *************************
!                                                           
!***********************************************************


rhs(i,j,indvars(2)) = rhs(i,j,indvars(2))  -  ( qst(i,j,indvarsst(10))*(d1_dif_rhou_dx_0_ijk)+&
                    qst(i,j,indvarsst(11))*(d1_dif_rhou_dy_0_ijk) ) 



!***********************************************************
!                                                           
! building source terms in RHS for d(rho v)/dt *************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([-visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x))]_1x)+deltayI*([-visc_t*(2.0_wp*deltayI*({v}_1y)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y)))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_dif_rhov_dxdx_0_0_im1jk_im1m1jk = q(i-1-1,j,indvars(3))

d2_dif_rhov_dxdx_0_0_im1jk_im1p1jk = q(i-1+1,j,indvars(3))

d2_dif_rhov_dxdx_0_0_im1jk = -&
          0.5_wp*d2_dif_rhov_dxdx_0_0_im1jk_im1m1jk+&
          0.5_wp*d2_dif_rhov_dxdx_0_0_im1jk_im1p1jk

d2_dif_rhov_dxdx_0_0_im1jk = d2_dif_rhov_dxdx_0_0_im1jk*param_float(1)

d2_dif_rhov_dxdx_0_0_ip1jk_ip1m1jk = q(i+1-1,j,indvars(3))

d2_dif_rhov_dxdx_0_0_ip1jk_ip1p1jk = q(i+1+1,j,indvars(3))

d2_dif_rhov_dxdx_0_0_ip1jk = -&
          0.5_wp*d2_dif_rhov_dxdx_0_0_ip1jk_ip1m1jk+&
          0.5_wp*d2_dif_rhov_dxdx_0_0_ip1jk_ip1p1jk

d2_dif_rhov_dxdx_0_0_ip1jk = d2_dif_rhov_dxdx_0_0_ip1jk*param_float(1)

d2_dif_rhov_dxdy_0_0_im1jk_im1jm1k = q(i-1,j-1,indvars(2))

d2_dif_rhov_dxdy_0_0_im1jk_im1jp1k = q(i-1,j+1,indvars(2))

d2_dif_rhov_dxdy_0_0_im1jk = -&
          0.5_wp*d2_dif_rhov_dxdy_0_0_im1jk_im1jm1k+&
          0.5_wp*d2_dif_rhov_dxdy_0_0_im1jk_im1jp1k

d2_dif_rhov_dxdy_0_0_im1jk = d2_dif_rhov_dxdy_0_0_im1jk*param_float(2)

d2_dif_rhov_dxdy_0_0_ip1jk_ip1jm1k = q(i+1,j-1,indvars(2))

d2_dif_rhov_dxdy_0_0_ip1jk_ip1jp1k = q(i+1,j+1,indvars(2))

d2_dif_rhov_dxdy_0_0_ip1jk = -&
          0.5_wp*d2_dif_rhov_dxdy_0_0_ip1jk_ip1jm1k+&
          0.5_wp*d2_dif_rhov_dxdy_0_0_ip1jk_ip1jp1k

d2_dif_rhov_dxdy_0_0_ip1jk = d2_dif_rhov_dxdy_0_0_ip1jk*param_float(2)

d1_dif_rhov_dx_0_im1jk = -(1.0_wp)*(1.0_wp+&
                    ((q(i-1,j,indvars(5))/1.0_wp*q(i-1,j,indvars(1)))**3.0_wp/((q(i-1,j,indvars(5))/1.0_wp*q(i-1,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i-1,j,indvars(5))/1.0_wp*q(i-1,j,indvars(1))))*param_float(1 + 5)*(qst(i-1,j,indvarsst(11))*(d2_dif_rhov_dxdy_0_0_im1jk)+&
                    qst(i-1,j,indvarsst(10))*(d2_dif_rhov_dxdx_0_0_im1jk))

d1_dif_rhov_dx_0_ip1jk = -(1.0_wp)*(1.0_wp+&
                    ((q(i+1,j,indvars(5))/1.0_wp*q(i+1,j,indvars(1)))**3.0_wp/((q(i+1,j,indvars(5))/1.0_wp*q(i+1,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i+1,j,indvars(5))/1.0_wp*q(i+1,j,indvars(1))))*param_float(1 + 5)*(qst(i+1,j,indvarsst(11))*(d2_dif_rhov_dxdy_0_0_ip1jk)+&
                    qst(i+1,j,indvarsst(10))*(d2_dif_rhov_dxdx_0_0_ip1jk))

d1_dif_rhov_dx_0_ijk = -&
          0.5_wp*d1_dif_rhov_dx_0_im1jk+&
          0.5_wp*d1_dif_rhov_dx_0_ip1jk

d1_dif_rhov_dx_0_ijk = d1_dif_rhov_dx_0_ijk*param_float(1)

d2_dif_rhov_dydx_0_0_ijm1k_im1jm1k = q(i-1,j-1,indvars(2))

d2_dif_rhov_dydx_0_0_ijm1k_ip1jm1k = q(i+1,j-1,indvars(2))

d2_dif_rhov_dydx_0_0_ijm1k = -&
          0.5_wp*d2_dif_rhov_dydx_0_0_ijm1k_im1jm1k+&
          0.5_wp*d2_dif_rhov_dydx_0_0_ijm1k_ip1jm1k

d2_dif_rhov_dydx_0_0_ijm1k = d2_dif_rhov_dydx_0_0_ijm1k*param_float(1)

d2_dif_rhov_dydx_0_0_ijp1k_im1jp1k = q(i-1,j+1,indvars(2))

d2_dif_rhov_dydx_0_0_ijp1k_ip1jp1k = q(i+1,j+1,indvars(2))

d2_dif_rhov_dydx_0_0_ijp1k = -&
          0.5_wp*d2_dif_rhov_dydx_0_0_ijp1k_im1jp1k+&
          0.5_wp*d2_dif_rhov_dydx_0_0_ijp1k_ip1jp1k

d2_dif_rhov_dydx_0_0_ijp1k = d2_dif_rhov_dydx_0_0_ijp1k*param_float(1)

d2_dif_rhov_dydy_0_0_ijm1k_ijm1m1k = q(i,j-1-1,indvars(3))

d2_dif_rhov_dydy_0_0_ijm1k_ijm1p1k = q(i,j-1+1,indvars(3))

d2_dif_rhov_dydy_0_0_ijm1k = -&
          0.5_wp*d2_dif_rhov_dydy_0_0_ijm1k_ijm1m1k+&
          0.5_wp*d2_dif_rhov_dydy_0_0_ijm1k_ijm1p1k

d2_dif_rhov_dydy_0_0_ijm1k = d2_dif_rhov_dydy_0_0_ijm1k*param_float(2)

d2_dif_rhov_dydy_0_0_ijp1k_ijp1m1k = q(i,j+1-1,indvars(3))

d2_dif_rhov_dydy_0_0_ijp1k_ijp1p1k = q(i,j+1+1,indvars(3))

d2_dif_rhov_dydy_0_0_ijp1k = -&
          0.5_wp*d2_dif_rhov_dydy_0_0_ijp1k_ijp1m1k+&
          0.5_wp*d2_dif_rhov_dydy_0_0_ijp1k_ijp1p1k

d2_dif_rhov_dydy_0_0_ijp1k = d2_dif_rhov_dydy_0_0_ijp1k*param_float(2)

d1_dif_rhov_dy_0_ijm1k = -(1.0_wp)*(1.0_wp+&
                    ((q(i,j-1,indvars(5))/1.0_wp*q(i,j-1,indvars(1)))**3.0_wp/((q(i,j-1,indvars(5))/1.0_wp*q(i,j-1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i,j-1,indvars(5))/1.0_wp*q(i,j-1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(i,j-1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_ijm1k)-&
                    2.0_wp/3.0_wp*(qst(i,j-1,indvarsst(10))*(d2_dif_rhov_dydx_0_0_ijm1k)+&
                    qst(i,j-1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_ijm1k)))

d1_dif_rhov_dy_0_ijp1k = -(1.0_wp)*(1.0_wp+&
                    ((q(i,j+1,indvars(5))/1.0_wp*q(i,j+1,indvars(1)))**3.0_wp/((q(i,j+1,indvars(5))/1.0_wp*q(i,j+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i,j+1,indvars(5))/1.0_wp*q(i,j+1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(i,j+1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_ijp1k)-&
                    2.0_wp/3.0_wp*(qst(i,j+1,indvarsst(10))*(d2_dif_rhov_dydx_0_0_ijp1k)+&
                    qst(i,j+1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_ijp1k)))

d1_dif_rhov_dy_0_ijk = -&
          0.5_wp*d1_dif_rhov_dy_0_ijm1k+&
          0.5_wp*d1_dif_rhov_dy_0_ijp1k

d1_dif_rhov_dy_0_ijk = d1_dif_rhov_dy_0_ijk*param_float(2)



!***********************************************************
!                                                           
! Update RHS terms for d(rho v)/dt *************************
!                                                           
!***********************************************************


rhs(i,j,indvars(3)) = rhs(i,j,indvars(3))  -  ( qst(i,j,indvarsst(10))*(d1_dif_rhov_dx_0_ijk)+&
                    qst(i,j,indvarsst(11))*(d1_dif_rhov_dy_0_ijk) ) 



!***********************************************************
!                                                           
! building source terms in RHS for d(rho et)/dt ************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([-kappa*deltaxI*({T}_1x)-u*(visc_t*(2.0_wp*deltaxI*({u}_1x)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y))))-v*(visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x)))]_1x)+deltayI*([-kappa*deltayI*({T}_1y)-u*(visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x)))-v*(visc_t*(2.0_wp*deltayI*({v}_1y)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y))))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_dif_et_dxdx_0_0_im1jk_im1m1jk = ((q(i-1-1,j,indvars(4))-&
                    0.5_wp*(q(i-1-1,j,indvars(2))*q(i-1-1,j,indvars(2))+&
                    q(i-1-1,j,indvars(3))*q(i-1-1,j,indvars(3)))))/(param_float(4 + 5))

d2_dif_et_dxdx_0_0_im1jk_im1p1jk = ((q(i-1+1,j,indvars(4))-&
                    0.5_wp*(q(i-1+1,j,indvars(2))*q(i-1+1,j,indvars(2))+&
                    q(i-1+1,j,indvars(3))*q(i-1+1,j,indvars(3)))))/(param_float(4 + 5))

d2_dif_et_dxdx_0_0_im1jk = -&
          0.5_wp*d2_dif_et_dxdx_0_0_im1jk_im1m1jk+&
          0.5_wp*d2_dif_et_dxdx_0_0_im1jk_im1p1jk

d2_dif_et_dxdx_0_0_im1jk = d2_dif_et_dxdx_0_0_im1jk*param_float(1)

d2_dif_et_dxdx_0_0_ip1jk_ip1m1jk = ((q(i+1-1,j,indvars(4))-&
                    0.5_wp*(q(i+1-1,j,indvars(2))*q(i+1-1,j,indvars(2))+&
                    q(i+1-1,j,indvars(3))*q(i+1-1,j,indvars(3)))))/(param_float(4 + 5))

d2_dif_et_dxdx_0_0_ip1jk_ip1p1jk = ((q(i+1+1,j,indvars(4))-&
                    0.5_wp*(q(i+1+1,j,indvars(2))*q(i+1+1,j,indvars(2))+&
                    q(i+1+1,j,indvars(3))*q(i+1+1,j,indvars(3)))))/(param_float(4 + 5))

d2_dif_et_dxdx_0_0_ip1jk = -&
          0.5_wp*d2_dif_et_dxdx_0_0_ip1jk_ip1m1jk+&
          0.5_wp*d2_dif_et_dxdx_0_0_ip1jk_ip1p1jk

d2_dif_et_dxdx_0_0_ip1jk = d2_dif_et_dxdx_0_0_ip1jk*param_float(1)

d1_dif_et_dx_0_im1jk = -param_float(2 + 5)*qst(i-1,j,indvarsst(10))*(d2_dif_et_dxdx_0_0_im1jk)-&
                    q(i-1,j,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(i-1,j,indvars(5))/1.0_wp*q(i-1,j,indvars(1)))**3.0_wp/((q(i-1,j,indvars(5))/1.0_wp*q(i-1,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i-1,j,indvars(5))/1.0_wp*q(i-1,j,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(i-1,j,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_im1jk)-&
                    2.0_wp/3.0_wp*(qst(i-1,j,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_im1jk)+&
                    qst(i-1,j,indvarsst(11))*(d2_dif_rhou_dxdy_0_0_im1jk))))-&
                    q(i-1,j,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(i-1,j,indvars(5))/1.0_wp*q(i-1,j,indvars(1)))**3.0_wp/((q(i-1,j,indvars(5))/1.0_wp*q(i-1,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i-1,j,indvars(5))/1.0_wp*q(i-1,j,indvars(1))))*param_float(1 + 5)*(qst(i-1,j,indvarsst(11))*(d2_dif_rhov_dxdy_0_0_im1jk)+&
                    qst(i-1,j,indvarsst(10))*(d2_dif_rhov_dxdx_0_0_im1jk)))

d1_dif_et_dx_0_ip1jk = -param_float(2 + 5)*qst(i+1,j,indvarsst(10))*(d2_dif_et_dxdx_0_0_ip1jk)-&
                    q(i+1,j,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(i+1,j,indvars(5))/1.0_wp*q(i+1,j,indvars(1)))**3.0_wp/((q(i+1,j,indvars(5))/1.0_wp*q(i+1,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i+1,j,indvars(5))/1.0_wp*q(i+1,j,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(i+1,j,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_ip1jk)-&
                    2.0_wp/3.0_wp*(qst(i+1,j,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_ip1jk)+&
                    qst(i+1,j,indvarsst(11))*(d2_dif_rhou_dxdy_0_0_ip1jk))))-&
                    q(i+1,j,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(i+1,j,indvars(5))/1.0_wp*q(i+1,j,indvars(1)))**3.0_wp/((q(i+1,j,indvars(5))/1.0_wp*q(i+1,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i+1,j,indvars(5))/1.0_wp*q(i+1,j,indvars(1))))*param_float(1 + 5)*(qst(i+1,j,indvarsst(11))*(d2_dif_rhov_dxdy_0_0_ip1jk)+&
                    qst(i+1,j,indvarsst(10))*(d2_dif_rhov_dxdx_0_0_ip1jk)))

d1_dif_et_dx_0_ijk = -&
          0.5_wp*d1_dif_et_dx_0_im1jk+&
          0.5_wp*d1_dif_et_dx_0_ip1jk

d1_dif_et_dx_0_ijk = d1_dif_et_dx_0_ijk*param_float(1)

d2_dif_et_dydy_0_0_ijm1k_ijm1m1k = ((q(i,j-1-1,indvars(4))-&
                    0.5_wp*(q(i,j-1-1,indvars(2))*q(i,j-1-1,indvars(2))+&
                    q(i,j-1-1,indvars(3))*q(i,j-1-1,indvars(3)))))/(param_float(4 + 5))

d2_dif_et_dydy_0_0_ijm1k_ijm1p1k = ((q(i,j-1+1,indvars(4))-&
                    0.5_wp*(q(i,j-1+1,indvars(2))*q(i,j-1+1,indvars(2))+&
                    q(i,j-1+1,indvars(3))*q(i,j-1+1,indvars(3)))))/(param_float(4 + 5))

d2_dif_et_dydy_0_0_ijm1k = -&
          0.5_wp*d2_dif_et_dydy_0_0_ijm1k_ijm1m1k+&
          0.5_wp*d2_dif_et_dydy_0_0_ijm1k_ijm1p1k

d2_dif_et_dydy_0_0_ijm1k = d2_dif_et_dydy_0_0_ijm1k*param_float(2)

d2_dif_et_dydy_0_0_ijp1k_ijp1m1k = ((q(i,j+1-1,indvars(4))-&
                    0.5_wp*(q(i,j+1-1,indvars(2))*q(i,j+1-1,indvars(2))+&
                    q(i,j+1-1,indvars(3))*q(i,j+1-1,indvars(3)))))/(param_float(4 + 5))

d2_dif_et_dydy_0_0_ijp1k_ijp1p1k = ((q(i,j+1+1,indvars(4))-&
                    0.5_wp*(q(i,j+1+1,indvars(2))*q(i,j+1+1,indvars(2))+&
                    q(i,j+1+1,indvars(3))*q(i,j+1+1,indvars(3)))))/(param_float(4 + 5))

d2_dif_et_dydy_0_0_ijp1k = -&
          0.5_wp*d2_dif_et_dydy_0_0_ijp1k_ijp1m1k+&
          0.5_wp*d2_dif_et_dydy_0_0_ijp1k_ijp1p1k

d2_dif_et_dydy_0_0_ijp1k = d2_dif_et_dydy_0_0_ijp1k*param_float(2)

d1_dif_et_dy_0_ijm1k = -param_float(2 + 5)*qst(i,j-1,indvarsst(11))*(d2_dif_et_dydy_0_0_ijm1k)-&
                    q(i,j-1,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(i,j-1,indvars(5))/1.0_wp*q(i,j-1,indvars(1)))**3.0_wp/((q(i,j-1,indvars(5))/1.0_wp*q(i,j-1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i,j-1,indvars(5))/1.0_wp*q(i,j-1,indvars(1))))*param_float(1 + 5)*(qst(i,j-1,indvarsst(11))*(d2_dif_rhou_dydy_0_0_ijm1k)+&
                    qst(i,j-1,indvarsst(10))*(d2_dif_rhou_dydx_0_0_ijm1k)))-&
                    q(i,j-1,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(i,j-1,indvars(5))/1.0_wp*q(i,j-1,indvars(1)))**3.0_wp/((q(i,j-1,indvars(5))/1.0_wp*q(i,j-1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i,j-1,indvars(5))/1.0_wp*q(i,j-1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(i,j-1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_ijm1k)-&
                    2.0_wp/3.0_wp*(qst(i,j-1,indvarsst(10))*(d2_dif_rhov_dydx_0_0_ijm1k)+&
                    qst(i,j-1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_ijm1k))))

d1_dif_et_dy_0_ijp1k = -param_float(2 + 5)*qst(i,j+1,indvarsst(11))*(d2_dif_et_dydy_0_0_ijp1k)-&
                    q(i,j+1,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(i,j+1,indvars(5))/1.0_wp*q(i,j+1,indvars(1)))**3.0_wp/((q(i,j+1,indvars(5))/1.0_wp*q(i,j+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i,j+1,indvars(5))/1.0_wp*q(i,j+1,indvars(1))))*param_float(1 + 5)*(qst(i,j+1,indvarsst(11))*(d2_dif_rhou_dydy_0_0_ijp1k)+&
                    qst(i,j+1,indvarsst(10))*(d2_dif_rhou_dydx_0_0_ijp1k)))-&
                    q(i,j+1,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(i,j+1,indvars(5))/1.0_wp*q(i,j+1,indvars(1)))**3.0_wp/((q(i,j+1,indvars(5))/1.0_wp*q(i,j+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i,j+1,indvars(5))/1.0_wp*q(i,j+1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(i,j+1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_ijp1k)-&
                    2.0_wp/3.0_wp*(qst(i,j+1,indvarsst(10))*(d2_dif_rhov_dydx_0_0_ijp1k)+&
                    qst(i,j+1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_ijp1k))))

d1_dif_et_dy_0_ijk = -&
          0.5_wp*d1_dif_et_dy_0_ijm1k+&
          0.5_wp*d1_dif_et_dy_0_ijp1k

d1_dif_et_dy_0_ijk = d1_dif_et_dy_0_ijk*param_float(2)



!***********************************************************
!                                                           
! Update RHS terms for d(rho et)/dt ************************
!                                                           
!***********************************************************


rhs(i,j,indvars(4)) = rhs(i,j,indvars(4))  -  ( qst(i,j,indvarsst(10))*(d1_dif_et_dx_0_ijk)+&
                    qst(i,j,indvarsst(11))*(d1_dif_et_dy_0_ijk) ) 



!***********************************************************
!                                                           
! building source terms in RHS for d(rho nut)/dt ***********
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! -ReI*Cb2*sigmaI*((deltaxI)**2*([rho*nut]_1x)*([nut]_1x)+(deltayI)**2*([rho*nut]_1y)*([nut]_1y))-Cb1*(1.0_wp-ft2)*SS*rho*nut+ReI*(Cw1*fw-Cb1/k**2.0_wp*ft2)*rho*nut**2.0_wp/eta**2.0_wp
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_dif_nut_dx_0_im1jk = q(i-1,j,indvars(1))*q(i-1,j,indvars(5))

d1_dif_nut_dx_0_ip1jk = q(i+1,j,indvars(1))*q(i+1,j,indvars(5))

d1_dif_nut_dx_0_ijk = -&
          0.5_wp*d1_dif_nut_dx_0_im1jk+&
          0.5_wp*d1_dif_nut_dx_0_ip1jk

d1_dif_nut_dx_0_ijk = d1_dif_nut_dx_0_ijk*param_float(1)

d1_dif_nut_dx_1_im1jk = q(i-1,j,indvars(5))

d1_dif_nut_dx_1_ip1jk = q(i+1,j,indvars(5))

d1_dif_nut_dx_1_ijk = -&
          0.5_wp*d1_dif_nut_dx_1_im1jk+&
          0.5_wp*d1_dif_nut_dx_1_ip1jk

d1_dif_nut_dx_1_ijk = d1_dif_nut_dx_1_ijk*param_float(1)

d1_dif_nut_dy_0_ijm1k = q(i,j-1,indvars(1))*q(i,j-1,indvars(5))

d1_dif_nut_dy_0_ijp1k = q(i,j+1,indvars(1))*q(i,j+1,indvars(5))

d1_dif_nut_dy_0_ijk = -&
          0.5_wp*d1_dif_nut_dy_0_ijm1k+&
          0.5_wp*d1_dif_nut_dy_0_ijp1k

d1_dif_nut_dy_0_ijk = d1_dif_nut_dy_0_ijk*param_float(2)

d1_dif_nut_dy_1_ijm1k = q(i,j-1,indvars(5))

d1_dif_nut_dy_1_ijp1k = q(i,j+1,indvars(5))

d1_dif_nut_dy_1_ijk = -&
          0.5_wp*d1_dif_nut_dy_1_ijm1k+&
          0.5_wp*d1_dif_nut_dy_1_ijp1k

d1_dif_nut_dy_1_ijk = d1_dif_nut_dy_1_ijk*param_float(2)



!***********************************************************
!                                                           
! Update RHS terms for d(rho nut)/dt ***********************
!                                                           
!***********************************************************


rhs(i,j,indvars(5)) = rhs(i,j,indvars(5))  -  ( -param_float(1 + 5)*param_float(7 + 5)*param_float(18 + 5)*((qst(i,j,indvarsst(10)))**2*(d1_dif_nut_dx_0_ijk)*(d1_dif_nut_dx_1_ijk)+&
                    (qst(i,j,indvarsst(11)))**2*(d1_dif_nut_dy_0_ijk)*(d1_dif_nut_dy_1_ijk))-&
                    param_float(6 + 5)*(1.0_wp-&
                    (param_float(16 + 5)*exp(-&
                    param_float(17 + 5)*(q(i,j,indvars(5))/1.0_wp*q(i,j,indvars(1)))**2.0_wp)))*qst(i,j,indvarsst(12))*q(i,j,indvars(1))*q(i,j,indvars(5))+&
                    param_float(1 + 5)*(param_float(10 + 5)*(((min(param_float(1 + 5)*(q(i,j,indvars(5))/(qst(i,j,indvarsst(12))*param_float(9 + 5)**2.0_wp*qst(i,j,indvarsst(2))**2.0_wp)),10.0_wp))+&
                    param_float(11 + 5)*((min(param_float(1 + 5)*(q(i,j,indvars(5))/(qst(i,j,indvarsst(12))*param_float(9 + 5)**2.0_wp*qst(i,j,indvarsst(2))**2.0_wp)),10.0_wp))**6.0_wp-&
                    (min(param_float(1 + 5)*(q(i,j,indvars(5))/(qst(i,j,indvarsst(12))*param_float(9 + 5)**2.0_wp*qst(i,j,indvarsst(2))**2.0_wp)),10.0_wp))))*((1.0_wp+&
                    param_float(12 + 5)**6.0_wp)/(((min(param_float(1 + 5)*(q(i,j,indvars(5))/(qst(i,j,indvarsst(12))*param_float(9 + 5)**2.0_wp*qst(i,j,indvarsst(2))**2.0_wp)),10.0_wp))+&
                    param_float(11 + 5)*((min(param_float(1 + 5)*(q(i,j,indvars(5))/(qst(i,j,indvarsst(12))*param_float(9 + 5)**2.0_wp*qst(i,j,indvarsst(2))**2.0_wp)),10.0_wp))**6.0_wp-&
                    (min(param_float(1 + 5)*(q(i,j,indvars(5))/(qst(i,j,indvarsst(12))*param_float(9 + 5)**2.0_wp*qst(i,j,indvarsst(2))**2.0_wp)),10.0_wp))))**6.0_wp+&
                    param_float(12 + 5)**6.0_wp))**(1.0_wp/6.0_wp))-&
                    param_float(6 + 5)/param_float(9 + 5)**2.0_wp*(param_float(16 + 5)*exp(-&
                    param_float(17 + 5)*(q(i,j,indvars(5))/1.0_wp*q(i,j,indvars(1)))**2.0_wp)))*q(i,j,indvars(1))*q(i,j,indvars(5))**2.0_wp/qst(i,j,indvarsst(2))**2.0_wp ) 

     enddo
   enddo

end subroutine cmprhs

subroutine cmpstored(param_float,ind,idloop,idarray,neq,neqst,sizeblck,q,qst,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

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


integer, intent(in) :: neq,neqst
integer, intent(in), dimension(6)   :: idloop, idarray
integer, intent(in), dimension(neq+neqst) :: ind
integer, intent(in), dimension(3) :: sizeblck
integer, intent(in) :: nvar_f(3),nvar_e(3)

real(wp), intent(in)    :: param_float(*)

real(wp),intent(inout) :: q(idarray(1):idarray(2),idarray(3):idarray(4),neq),&
                        qst(idarray(1):idarray(2),idarray(3):idarray(4),neqst)

real(wp),intent(inout) :: qface_i(idarray(3):idarray(4),nvar_f(1)),&
                       qface_j(idarray(1):idarray(2),nvar_f(2)),&
                       qface_k(1),&
                       qedge_ij(1),&
                       qedge_jk(1),&
                       qedge_ik(1)

! LOCAL VARIABLES



 real(wp) ::  d1_stemp_dx_0_im2jk,d1_stemp_dx_0_im1jk,d1_stemp_dx_0_ip0jk,d1_stemp_dx_0_ip1jk,d1_stemp_dx_0_ip2jk &
            ,d1_stemp_dx_0_ijk &
            ,d1_stemp_dy_0_ijm2k,d1_stemp_dy_0_ijm1k,d1_stemp_dy_0_ijp0k,d1_stemp_dy_0_ijp1k,d1_stemp_dy_0_ijp2k &
            ,d1_stemp_dy_0_ijk 

integer :: i,j,k,nq
integer :: bi,bj,bk
integer :: biend,bjend,bkend
integer :: size_bi,size_bj,size_bk

integer :: indvars(neq),indvarsst(neqst)


  indvars   = ind(1:neq)
  indvarsst = ind(1+neq:neq+neqst) 



!***********************************************************
!                                                           
! Start computing stored variables *************************
!                                                           
!***********************************************************


     do j=idloop(3),idloop(4) 
 
      do i=idloop(1),idloop(2) 


!***********************************************************
!                                                           
! building stored variable stemp ***************************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (2.0_wp*(dabs(0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_im2jk = q(i-2,j,indvars(3))

d1_stemp_dx_0_im1jk = q(i-1,j,indvars(3))

d1_stemp_dx_0_ip1jk = q(i+1,j,indvars(3))

d1_stemp_dx_0_ip2jk = q(i+2,j,indvars(3))

d1_stemp_dx_0_ijk = 0.08333333333333333_wp*d1_stemp_dx_0_im2jk-&
          0.666666666666667_wp*d1_stemp_dx_0_im1jk+&
          0.666666666666667_wp*d1_stemp_dx_0_ip1jk-&
          0.08333333333333333_wp*d1_stemp_dx_0_ip2jk

d1_stemp_dx_0_ijk = d1_stemp_dx_0_ijk*param_float(1)

d1_stemp_dy_0_ijm2k = q(i,j-2,indvars(2))

d1_stemp_dy_0_ijm1k = q(i,j-1,indvars(2))

d1_stemp_dy_0_ijp1k = q(i,j+1,indvars(2))

d1_stemp_dy_0_ijp2k = q(i,j+2,indvars(2))

d1_stemp_dy_0_ijk = 0.08333333333333333_wp*d1_stemp_dy_0_ijm2k-&
          0.666666666666667_wp*d1_stemp_dy_0_ijm1k+&
          0.666666666666667_wp*d1_stemp_dy_0_ijp1k-&
          0.08333333333333333_wp*d1_stemp_dy_0_ijp2k

d1_stemp_dy_0_ijk = d1_stemp_dy_0_ijk*param_float(2)



!***********************************************************
!                                                           
! Update stored variables  ( 2.0_wp*( dabs(  0.5_wp*( deltayI*( [u]_1y ) - deltaxI*( [v]_1x ) ) ) ) )  
!                                                           
!***********************************************************


qst(i,j,indvarsst(4)) =  (2.0_wp*(dabs(0.5_wp*(qst(i,j,indvarsst(11))*(d1_stemp_dy_0_ijk)-&
                    qst(i,j,indvarsst(10))*(d1_stemp_dx_0_ijk)))))



!***********************************************************
!                                                           
! building stored variable SS ******************************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (stemp+fv2*ReI*nut/(k**2.0_wp*eta**2.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update stored variables  ( stemp + fv2*ReI*nut/(k**2.0_wp*eta**2.0_wp) )  
!                                                           
!***********************************************************


qst(i,j,indvarsst(12)) =  (qst(i,j,indvarsst(4))+&
                    (1.0_wp-&
                    (q(i,j,indvars(5))/1.0_wp*q(i,j,indvars(1)))/(1.0_wp+&
                    (q(i,j,indvars(5))/1.0_wp*q(i,j,indvars(1)))*((q(i,j,indvars(5))/1.0_wp*q(i,j,indvars(1)))**3.0_wp/((q(i,j,indvars(5))/1.0_wp*q(i,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))))*param_float(1 + 5)*q(i,j,indvars(5))/(param_float(9 + 5)**2.0_wp*qst(i,j,indvarsst(2))**2.0_wp))



!***********************************************************
!                                                           
! building stored variable tau_wall ************************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (visc_t*([u]_1y)*deltayI)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update stored variables  ( visc_t*( [u]_1y )*deltayI ) ***
!                                                           
!***********************************************************


qst(i,j,indvarsst(13)) =  ((1.0_wp)*(1.0_wp+&
                    ((q(i,j,indvars(5))/1.0_wp*q(i,j,indvars(1)))**3.0_wp/((q(i,j,indvars(5))/1.0_wp*q(i,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i,j,indvars(5))/1.0_wp*q(i,j,indvars(1))))*param_float(1 + 5)*(d1_stemp_dy_0_ijk)*qst(i,j,indvarsst(11)))



!***********************************************************
!                                                           
! building stored variable visc_SA *************************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! nut*rho
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update stored variables nut*rho **************************
!                                                           
!***********************************************************


qst(i,j,indvarsst(14)) =  q(i,j,indvars(5))*q(i,j,indvars(1))



!***********************************************************
!                                                           
! building stored variable visc_turb ***********************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! nut*rho*fv1
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update stored variables nut*rho*fv1 **********************
!                                                           
!***********************************************************


qst(i,j,indvarsst(15)) =  q(i,j,indvars(5))*q(i,j,indvars(1))*((q(i,j,indvars(5))/1.0_wp*q(i,j,indvars(1)))**3.0_wp/((q(i,j,indvars(5))/1.0_wp*q(i,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))



!***********************************************************
!                                                           
! building stored variable Pressure ************************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (gamma_m1)*rho*(e)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update stored variables (gamma_m1)*rho*(e) ***************
!                                                           
!***********************************************************


qst(i,j,indvarsst(16)) =  (param_float(3 + 5))*q(i,j,indvars(1))*((q(i,j,indvars(4))-&
                    0.5_wp*(q(i,j,indvars(2))*q(i,j,indvars(2))+&
                    q(i,j,indvars(3))*q(i,j,indvars(3)))))

     enddo
   enddo

end subroutine cmpstored

subroutine cmpstoredstatic(param_float,ind,idloop,idarray,neq,neqst,sizeblck,q,qst,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

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


integer, intent(in) :: neq,neqst
integer, intent(in), dimension(6)   :: idloop, idarray
integer, intent(in), dimension(neq+neqst) :: ind
integer, intent(in), dimension(3) :: sizeblck
integer, intent(in) :: nvar_f(3),nvar_e(3)

real(wp), intent(in)    :: param_float(*)

real(wp),intent(inout) :: q(idarray(1):idarray(2),idarray(3):idarray(4),neq),&
                        qst(idarray(1):idarray(2),idarray(3):idarray(4),neqst)

real(wp),intent(inout) :: qface_i(idarray(3):idarray(4),nvar_f(1)),&
                       qface_j(idarray(1):idarray(2),nvar_f(2)),&
                       qface_k(1),&
                       qedge_ij(1),&
                       qedge_jk(1),&
                       qedge_ik(1)

! LOCAL VARIABLES



 real(wp) ::  d1_detady_dy_0_ijm2k,d1_detady_dy_0_ijm1k,d1_detady_dy_0_ijp0k,d1_detady_dy_0_ijp1k,d1_detady_dy_0_ijp2k &
            ,d1_detady_dy_0_ijk &
            ,d1_dksidy_dy_0_ijm2k,d1_dksidy_dy_0_ijm1k,d1_dksidy_dy_0_ijp0k,d1_dksidy_dy_0_ijp1k,d1_dksidy_dy_0_ijp2k &
            ,d1_dksidy_dy_0_ijk &
            ,d1_detadx_dx_0_im2jk,d1_detadx_dx_0_im1jk,d1_detadx_dx_0_ip0jk,d1_detadx_dx_0_ip1jk,d1_detadx_dx_0_ip2jk &
            ,d1_detadx_dx_0_ijk &
            ,d1_dksidx_dx_0_im2jk,d1_dksidx_dx_0_im1jk,d1_dksidx_dx_0_ip0jk,d1_dksidx_dx_0_ip1jk,d1_dksidx_dx_0_ip2jk &
            ,d1_dksidx_dx_0_ijk 

integer :: i,j,k,nq
integer :: bi,bj,bk
integer :: biend,bjend,bkend
integer :: size_bi,size_bj,size_bk

integer :: indvars(neq),indvarsst(neqst)


  indvars   = ind(1:neq)
  indvarsst = ind(1+neq:neq+neqst) 



!***********************************************************
!                                                           
! Start computing stored variables *************************
!                                                           
!***********************************************************


     do j=idloop(3),idloop(4) 
 
      do i=idloop(1),idloop(2) 


!***********************************************************
!                                                           
! building stored variable d *******************************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! d
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update stored variables  d  ******************************
!                                                           
!***********************************************************


qst(i,j,indvarsst(1)) =  qst(i,j,indvarsst(1))



!***********************************************************
!                                                           
! building stored variable eta *****************************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! eta
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update stored variables  eta  ****************************
!                                                           
!***********************************************************


qst(i,j,indvarsst(2)) =  qst(i,j,indvarsst(2))



!***********************************************************
!                                                           
! building stored variable ksi *****************************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ksi
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update stored variables  ksi  ****************************
!                                                           
!***********************************************************


qst(i,j,indvarsst(3)) =  qst(i,j,indvarsst(3))



!***********************************************************
!                                                           
! building stored variable symm ****************************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ((sign(1.0_wp,ksi)-1.0_wp)/(-2.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update stored variables  ( ( sign(1.0_wp, ksi) - 1.0_wp ) /(-2.0_wp) ) 
!                                                           
!***********************************************************


qst(i,j,indvarsst(5)) =  ((sign(1.0_wp,qst(i,j,indvarsst(3)))-&
                    1.0_wp)/(-&
                    2.0_wp))



!***********************************************************
!                                                           
! building stored variable detady **************************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [eta]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_detady_dy_0_ijm2k = qst(i,j-2,indvarsst(2))

d1_detady_dy_0_ijm1k = qst(i,j-1,indvarsst(2))

d1_detady_dy_0_ijp1k = qst(i,j+1,indvarsst(2))

d1_detady_dy_0_ijp2k = qst(i,j+2,indvarsst(2))

d1_detady_dy_0_ijk = 0.08333333333333333_wp*d1_detady_dy_0_ijm2k-&
          0.666666666666667_wp*d1_detady_dy_0_ijm1k+&
          0.666666666666667_wp*d1_detady_dy_0_ijp1k-&
          0.08333333333333333_wp*d1_detady_dy_0_ijp2k

d1_detady_dy_0_ijk = d1_detady_dy_0_ijk*param_float(2)



!***********************************************************
!                                                           
! Update stored variables  [ eta ]_1y  *********************
!                                                           
!***********************************************************


qst(i,j,indvarsst(6)) =  d1_detady_dy_0_ijk



!***********************************************************
!                                                           
! building stored variable dksidy **************************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [ksi]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_dksidy_dy_0_ijm2k = qst(i,j-2,indvarsst(3))

d1_dksidy_dy_0_ijm1k = qst(i,j-1,indvarsst(3))

d1_dksidy_dy_0_ijp1k = qst(i,j+1,indvarsst(3))

d1_dksidy_dy_0_ijp2k = qst(i,j+2,indvarsst(3))

d1_dksidy_dy_0_ijk = 0.08333333333333333_wp*d1_dksidy_dy_0_ijm2k-&
          0.666666666666667_wp*d1_dksidy_dy_0_ijm1k+&
          0.666666666666667_wp*d1_dksidy_dy_0_ijp1k-&
          0.08333333333333333_wp*d1_dksidy_dy_0_ijp2k

d1_dksidy_dy_0_ijk = d1_dksidy_dy_0_ijk*param_float(2)



!***********************************************************
!                                                           
! Update stored variables  [ ksi ]_1y  *********************
!                                                           
!***********************************************************


qst(i,j,indvarsst(7)) =  d1_dksidy_dy_0_ijk



!***********************************************************
!                                                           
! building stored variable detadx **************************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [eta]_1x
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_detadx_dx_0_im2jk = qst(i-2,j,indvarsst(2))

d1_detadx_dx_0_im1jk = qst(i-1,j,indvarsst(2))

d1_detadx_dx_0_ip1jk = qst(i+1,j,indvarsst(2))

d1_detadx_dx_0_ip2jk = qst(i+2,j,indvarsst(2))

d1_detadx_dx_0_ijk = 0.08333333333333333_wp*d1_detadx_dx_0_im2jk-&
          0.666666666666667_wp*d1_detadx_dx_0_im1jk+&
          0.666666666666667_wp*d1_detadx_dx_0_ip1jk-&
          0.08333333333333333_wp*d1_detadx_dx_0_ip2jk

d1_detadx_dx_0_ijk = d1_detadx_dx_0_ijk*param_float(1)



!***********************************************************
!                                                           
! Update stored variables  [ eta ]_1x  *********************
!                                                           
!***********************************************************


qst(i,j,indvarsst(8)) =  d1_detadx_dx_0_ijk



!***********************************************************
!                                                           
! building stored variable dksidx **************************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ([ksi]_1x)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_dksidx_dx_0_im2jk = qst(i-2,j,indvarsst(3))

d1_dksidx_dx_0_im1jk = qst(i-1,j,indvarsst(3))

d1_dksidx_dx_0_ip1jk = qst(i+1,j,indvarsst(3))

d1_dksidx_dx_0_ip2jk = qst(i+2,j,indvarsst(3))

d1_dksidx_dx_0_ijk = 0.08333333333333333_wp*d1_dksidx_dx_0_im2jk-&
          0.666666666666667_wp*d1_dksidx_dx_0_im1jk+&
          0.666666666666667_wp*d1_dksidx_dx_0_ip1jk-&
          0.08333333333333333_wp*d1_dksidx_dx_0_ip2jk

d1_dksidx_dx_0_ijk = d1_dksidx_dx_0_ijk*param_float(1)



!***********************************************************
!                                                           
! Update stored variables  ( [ ksi ]_1x )  *****************
!                                                           
!***********************************************************


qst(i,j,indvarsst(9)) =  (d1_dksidx_dx_0_ijk)



!***********************************************************
!                                                           
! building stored variable deltaxI *************************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! 1.0_wp/(dksidx)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update stored variables  1.0_wp / ( dksidx )  ************
!                                                           
!***********************************************************


qst(i,j,indvarsst(10)) =  1.0_wp/(qst(i,j,indvarsst(9)))



!***********************************************************
!                                                           
! building stored variable deltayI *************************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! 1.0_wp/(detady)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update stored variables  1.0_wp / ( detady )  ************
!                                                           
!***********************************************************


qst(i,j,indvarsst(11)) =  1.0_wp/(qst(i,j,indvarsst(6)))

     enddo
   enddo

end subroutine cmpstoredstatic
