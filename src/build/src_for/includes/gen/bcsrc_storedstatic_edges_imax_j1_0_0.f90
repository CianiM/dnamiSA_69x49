














!====================================================================================================
!
! General Boundary Conditions: compute Eqns with Boundary Scheme
!     
!====================================================================================================

subroutine boundarySchemestoredstatic_edges_imax_j1_0_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

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



 real(wp) ::  d1_detady_dy_0_nxp2p01m2p0p0k,d1_detady_dy_0_nxp2p01m2p0p1k,d1_detady_dy_0_nxp2p01m2p0p2k &
            ,d1_detady_dy_0_nxp2p01m2p0k &
            ,d1_dksidy_dy_0_nxp2p01m2p0p0k,d1_dksidy_dy_0_nxp2p01m2p0p1k,d1_dksidy_dy_0_nxp2p01m2p0p2k &
            ,d1_dksidy_dy_0_nxp2p01m2p0k &
            ,d1_detadx_dx_0_nxp2p0p01m2p0k,d1_detadx_dx_0_nxp2p0m11m2p0k,d1_detadx_dx_0_nxp2p0m21m2p0k &
            ,d1_detadx_dx_0_nxp2p01m2p0k &
            ,d1_dksidx_dx_0_nxp2p0p01m2p0k,d1_dksidx_dx_0_nxp2p0m11m2p0k,d1_dksidx_dx_0_nxp2p0m21m2p0k &
            ,d1_dksidx_dx_0_nxp2p01m2p0k 

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
! BC layer: 0 0 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None d ********
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! d
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None d *********************
!                                                           
!***********************************************************


qst(nx+2+0,1-2+0,indvarsst(1)) =  qst(nx+2+0,1-2+0,indvarsst(1))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None eta ******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! eta
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None eta *******************
!                                                           
!***********************************************************


qst(nx+2+0,1-2+0,indvarsst(2)) =  qst(nx+2+0,1-2+0,indvarsst(2))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None ksi ******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ksi
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None ksi *******************
!                                                           
!***********************************************************


qst(nx+2+0,1-2+0,indvarsst(3)) =  qst(nx+2+0,1-2+0,indvarsst(3))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None symm *****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ((sign(1.0_wp,ksi)-1.0_wp)/(-2.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None symm ******************
!                                                           
!***********************************************************


qst(nx+2+0,1-2+0,indvarsst(5)) =  ((sign(1.0_wp,qst(nx+2+0,1-2+0,indvarsst(3)))-&
                    1.0_wp)/(-&
                    2.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None detady ***
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [eta]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_detady_dy_0_nxp2p01m2p0p0k = qst(nx+2+0,1-2+0+0,indvarsst(2))

d1_detady_dy_0_nxp2p01m2p0p1k = qst(nx+2+0,1-2+0+1,indvarsst(2))

d1_detady_dy_0_nxp2p01m2p0p2k = qst(nx+2+0,1-2+0+2,indvarsst(2))

d1_detady_dy_0_nxp2p01m2p0k = -&
          1.5_wp*d1_detady_dy_0_nxp2p01m2p0p0k+&
          2.0_wp*d1_detady_dy_0_nxp2p01m2p0p1k-&
          0.5_wp*d1_detady_dy_0_nxp2p01m2p0p2k

d1_detady_dy_0_nxp2p01m2p0k = d1_detady_dy_0_nxp2p01m2p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None detady ****************
!                                                           
!***********************************************************


qst(nx+2+0,1-2+0,indvarsst(6)) =  d1_detady_dy_0_nxp2p01m2p0k



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None dksidy ***
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [ksi]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_dksidy_dy_0_nxp2p01m2p0p0k = qst(nx+2+0,1-2+0+0,indvarsst(3))

d1_dksidy_dy_0_nxp2p01m2p0p1k = qst(nx+2+0,1-2+0+1,indvarsst(3))

d1_dksidy_dy_0_nxp2p01m2p0p2k = qst(nx+2+0,1-2+0+2,indvarsst(3))

d1_dksidy_dy_0_nxp2p01m2p0k = -&
          1.5_wp*d1_dksidy_dy_0_nxp2p01m2p0p0k+&
          2.0_wp*d1_dksidy_dy_0_nxp2p01m2p0p1k-&
          0.5_wp*d1_dksidy_dy_0_nxp2p01m2p0p2k

d1_dksidy_dy_0_nxp2p01m2p0k = d1_dksidy_dy_0_nxp2p01m2p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None dksidy ****************
!                                                           
!***********************************************************


qst(nx+2+0,1-2+0,indvarsst(7)) =  d1_dksidy_dy_0_nxp2p01m2p0k



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None detadx ***
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [eta]_1x
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_detadx_dx_0_nxp2p0p01m2p0k = qst(nx+2+0+0,1-2+0,indvarsst(2))

d1_detadx_dx_0_nxp2p0m11m2p0k = qst(nx+2+0-1,1-2+0,indvarsst(2))

d1_detadx_dx_0_nxp2p0m21m2p0k = qst(nx+2+0-2,1-2+0,indvarsst(2))

d1_detadx_dx_0_nxp2p01m2p0k = 1.5_wp*d1_detadx_dx_0_nxp2p0p01m2p0k-&
          2.0_wp*d1_detadx_dx_0_nxp2p0m11m2p0k+&
          0.5_wp*d1_detadx_dx_0_nxp2p0m21m2p0k

d1_detadx_dx_0_nxp2p01m2p0k = d1_detadx_dx_0_nxp2p01m2p0k*param_float(1)



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None detadx ****************
!                                                           
!***********************************************************


qst(nx+2+0,1-2+0,indvarsst(8)) =  d1_detadx_dx_0_nxp2p01m2p0k



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None dksidx ***
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ([ksi]_1x)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_dksidx_dx_0_nxp2p0p01m2p0k = qst(nx+2+0+0,1-2+0,indvarsst(3))

d1_dksidx_dx_0_nxp2p0m11m2p0k = qst(nx+2+0-1,1-2+0,indvarsst(3))

d1_dksidx_dx_0_nxp2p0m21m2p0k = qst(nx+2+0-2,1-2+0,indvarsst(3))

d1_dksidx_dx_0_nxp2p01m2p0k = 1.5_wp*d1_dksidx_dx_0_nxp2p0p01m2p0k-&
          2.0_wp*d1_dksidx_dx_0_nxp2p0m11m2p0k+&
          0.5_wp*d1_dksidx_dx_0_nxp2p0m21m2p0k

d1_dksidx_dx_0_nxp2p01m2p0k = d1_dksidx_dx_0_nxp2p01m2p0k*param_float(1)



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None dksidx ****************
!                                                           
!***********************************************************


qst(nx+2+0,1-2+0,indvarsst(9)) =  (d1_dksidx_dx_0_nxp2p01m2p0k)



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None deltaxI **
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (1.0_wp/(dksidx))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None deltaxI ***************
!                                                           
!***********************************************************


qst(nx+2+0,1-2+0,indvarsst(10)) =  (1.0_wp/(qst(nx+2+0,1-2+0,indvarsst(9))))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None deltayI **
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (1.0_wp/(detady))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None deltayI ***************
!                                                           
!***********************************************************


qst(nx+2+0,1-2+0,indvarsst(11)) =  (1.0_wp/(qst(nx+2+0,1-2+0,indvarsst(6))))


    enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO

 end subroutine boundarySchemestoredstatic_edges_imax_j1_0_0


