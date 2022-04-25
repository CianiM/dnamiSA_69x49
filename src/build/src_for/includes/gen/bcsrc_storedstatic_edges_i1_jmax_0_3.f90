














!====================================================================================================
!
! General Boundary Conditions: compute Eqns with Boundary Scheme
!     
!====================================================================================================

subroutine boundarySchemestoredstatic_edges_i1_jmax_0_3(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

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



 real(wp) ::  d1_detady_dy_0_1m4p0nyp4m3m1k,d1_detady_dy_0_1m4p0nyp4m3p0k,d1_detady_dy_0_1m4p0nyp4m3p1k &
            ,d1_detady_dy_0_1m4p0nyp4m3k &
            ,d1_dksidy_dy_0_1m4p0nyp4m3m1k,d1_dksidy_dy_0_1m4p0nyp4m3p0k,d1_dksidy_dy_0_1m4p0nyp4m3p1k &
            ,d1_dksidy_dy_0_1m4p0nyp4m3k &
            ,d1_detadx_dx_0_1m4p0p0nyp4m3k,d1_detadx_dx_0_1m4p0p1nyp4m3k,d1_detadx_dx_0_1m4p0p2nyp4m3k &
            ,d1_detadx_dx_0_1m4p0nyp4m3k &
            ,d1_dksidx_dx_0_1m4p0p0nyp4m3k,d1_dksidx_dx_0_1m4p0p1nyp4m3k,d1_dksidx_dx_0_1m4p0p2nyp4m3k &
            ,d1_dksidx_dx_0_1m4p0nyp4m3k 

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
! Start building layers for BC : i1 jmax None **************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 0 3 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 0 3 None d ********
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! d
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 3 None d *********************
!                                                           
!***********************************************************


qst(1-4+0,ny+4-3,indvarsst(1)) =  qst(1-4+0,ny+4-3,indvarsst(1))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 3 None eta ******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! eta
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 3 None eta *******************
!                                                           
!***********************************************************


qst(1-4+0,ny+4-3,indvarsst(2)) =  qst(1-4+0,ny+4-3,indvarsst(2))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 3 None ksi ******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ksi
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 3 None ksi *******************
!                                                           
!***********************************************************


qst(1-4+0,ny+4-3,indvarsst(3)) =  qst(1-4+0,ny+4-3,indvarsst(3))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 3 None symm *****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ((sign(1.0_wp,ksi)-1.0_wp)/(-2.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 3 None symm ******************
!                                                           
!***********************************************************


qst(1-4+0,ny+4-3,indvarsst(5)) =  ((sign(1.0_wp,qst(1-4+0,ny+4-3,indvarsst(3)))-&
                    1.0_wp)/(-&
                    2.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 3 None detady ***
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [eta]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_detady_dy_0_1m4p0nyp4m3m1k = qst(1-4+0,ny+4-3-1,indvarsst(2))

d1_detady_dy_0_1m4p0nyp4m3p1k = qst(1-4+0,ny+4-3+1,indvarsst(2))

d1_detady_dy_0_1m4p0nyp4m3k = -&
          0.5_wp*d1_detady_dy_0_1m4p0nyp4m3m1k+&
          0.5_wp*d1_detady_dy_0_1m4p0nyp4m3p1k

d1_detady_dy_0_1m4p0nyp4m3k = d1_detady_dy_0_1m4p0nyp4m3k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 3 None detady ****************
!                                                           
!***********************************************************


qst(1-4+0,ny+4-3,indvarsst(6)) =  d1_detady_dy_0_1m4p0nyp4m3k



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 3 None dksidy ***
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [ksi]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_dksidy_dy_0_1m4p0nyp4m3m1k = qst(1-4+0,ny+4-3-1,indvarsst(3))

d1_dksidy_dy_0_1m4p0nyp4m3p1k = qst(1-4+0,ny+4-3+1,indvarsst(3))

d1_dksidy_dy_0_1m4p0nyp4m3k = -&
          0.5_wp*d1_dksidy_dy_0_1m4p0nyp4m3m1k+&
          0.5_wp*d1_dksidy_dy_0_1m4p0nyp4m3p1k

d1_dksidy_dy_0_1m4p0nyp4m3k = d1_dksidy_dy_0_1m4p0nyp4m3k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 3 None dksidy ****************
!                                                           
!***********************************************************


qst(1-4+0,ny+4-3,indvarsst(7)) =  d1_dksidy_dy_0_1m4p0nyp4m3k



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 3 None detadx ***
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! [eta]_1x
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_detadx_dx_0_1m4p0p0nyp4m3k = qst(1-4+0+0,ny+4-3,indvarsst(2))

d1_detadx_dx_0_1m4p0p1nyp4m3k = qst(1-4+0+1,ny+4-3,indvarsst(2))

d1_detadx_dx_0_1m4p0p2nyp4m3k = qst(1-4+0+2,ny+4-3,indvarsst(2))

d1_detadx_dx_0_1m4p0nyp4m3k = -&
          1.5_wp*d1_detadx_dx_0_1m4p0p0nyp4m3k+&
          2.0_wp*d1_detadx_dx_0_1m4p0p1nyp4m3k-&
          0.5_wp*d1_detadx_dx_0_1m4p0p2nyp4m3k

d1_detadx_dx_0_1m4p0nyp4m3k = d1_detadx_dx_0_1m4p0nyp4m3k*param_float(1)



!***********************************************************
!                                                           
! Update BC terms for layer 0 3 None detadx ****************
!                                                           
!***********************************************************


qst(1-4+0,ny+4-3,indvarsst(8)) =  d1_detadx_dx_0_1m4p0nyp4m3k



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 3 None dksidx ***
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ([ksi]_1x)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_dksidx_dx_0_1m4p0p0nyp4m3k = qst(1-4+0+0,ny+4-3,indvarsst(3))

d1_dksidx_dx_0_1m4p0p1nyp4m3k = qst(1-4+0+1,ny+4-3,indvarsst(3))

d1_dksidx_dx_0_1m4p0p2nyp4m3k = qst(1-4+0+2,ny+4-3,indvarsst(3))

d1_dksidx_dx_0_1m4p0nyp4m3k = -&
          1.5_wp*d1_dksidx_dx_0_1m4p0p0nyp4m3k+&
          2.0_wp*d1_dksidx_dx_0_1m4p0p1nyp4m3k-&
          0.5_wp*d1_dksidx_dx_0_1m4p0p2nyp4m3k

d1_dksidx_dx_0_1m4p0nyp4m3k = d1_dksidx_dx_0_1m4p0nyp4m3k*param_float(1)



!***********************************************************
!                                                           
! Update BC terms for layer 0 3 None dksidx ****************
!                                                           
!***********************************************************


qst(1-4+0,ny+4-3,indvarsst(9)) =  (d1_dksidx_dx_0_1m4p0nyp4m3k)



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 3 None deltaxI **
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! 1.0_wp/(dksidx)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 3 None deltaxI ***************
!                                                           
!***********************************************************


qst(1-4+0,ny+4-3,indvarsst(10)) =  1.0_wp/(qst(1-4+0,ny+4-3,indvarsst(9)))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 3 None deltayI **
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! 1.0_wp/(detady)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 3 None deltayI ***************
!                                                           
!***********************************************************


qst(1-4+0,ny+4-3,indvarsst(11)) =  1.0_wp/(qst(1-4+0,ny+4-3,indvarsst(6)))


    enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO

 end subroutine boundarySchemestoredstatic_edges_i1_jmax_0_3


