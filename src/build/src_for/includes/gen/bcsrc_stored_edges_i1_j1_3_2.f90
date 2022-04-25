














!====================================================================================================
!
! General Boundary Conditions: compute Eqns with Boundary Scheme
!     
!====================================================================================================

subroutine boundarySchemestored_edges_i1_j1_3_2(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

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



 real(wp) ::  d1_stemp_dx_0_1m4p3m11m4p2k,d1_stemp_dx_0_1m4p3p01m4p2k,d1_stemp_dx_0_1m4p3p11m4p2k &
            ,d1_stemp_dx_0_1m4p31m4p2k &
            ,d1_stemp_dy_0_1m4p31m4p2m1k,d1_stemp_dy_0_1m4p31m4p2p0k,d1_stemp_dy_0_1m4p31m4p2p1k &
            ,d1_stemp_dy_0_1m4p31m4p2k 

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
! BC layer: 3 2 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 3 2 None stemp ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (((0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))**2+(0.5_wp*(deltaxI*([v]_1x)-deltayI*([u]_1y)))**2)*2)**0.5
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_1m4p3m11m4p2k = q(1-4+3-1,1-4+2,indvars(3))

d1_stemp_dx_0_1m4p3p11m4p2k = q(1-4+3+1,1-4+2,indvars(3))

d1_stemp_dx_0_1m4p31m4p2k = -&
          0.5_wp*d1_stemp_dx_0_1m4p3m11m4p2k+&
          0.5_wp*d1_stemp_dx_0_1m4p3p11m4p2k

d1_stemp_dx_0_1m4p31m4p2k = d1_stemp_dx_0_1m4p31m4p2k*param_float(1)

d1_stemp_dy_0_1m4p31m4p2m1k = q(1-4+3,1-4+2-1,indvars(2))

d1_stemp_dy_0_1m4p31m4p2p1k = q(1-4+3,1-4+2+1,indvars(2))

d1_stemp_dy_0_1m4p31m4p2k = -&
          0.5_wp*d1_stemp_dy_0_1m4p31m4p2m1k+&
          0.5_wp*d1_stemp_dy_0_1m4p31m4p2p1k

d1_stemp_dy_0_1m4p31m4p2k = d1_stemp_dy_0_1m4p31m4p2k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 3 2 None stemp *****************
!                                                           
!***********************************************************


qst(1-4+3,1-4+2,indvarsst(4)) =  (((0.5_wp*(qst(1-4+3,1-4+2,indvarsst(11))*(d1_stemp_dy_0_1m4p31m4p2k)-&
                    qst(1-4+3,1-4+2,indvarsst(10))*(d1_stemp_dx_0_1m4p31m4p2k)))**2+&
                    (0.5_wp*(qst(1-4+3,1-4+2,indvarsst(10))*(d1_stemp_dx_0_1m4p31m4p2k)-&
                    qst(1-4+3,1-4+2,indvarsst(11))*(d1_stemp_dy_0_1m4p31m4p2k)))**2)*2)**0.5



!***********************************************************
!                                                           
! building source terms in RHS for layer 3 2 None SS *******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (stemp+ReI*nut/(k**2*eta**2))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 3 2 None SS ********************
!                                                           
!***********************************************************


qst(1-4+3,1-4+2,indvarsst(12)) =  (qst(1-4+3,1-4+2,indvarsst(4))+&
                    param_float(1 + 5)*q(1-4+3,1-4+2,indvars(5))/(param_float(9 + 5)**2*qst(1-4+3,1-4+2,indvarsst(2))**2))



!***********************************************************
!                                                           
! building source terms in RHS for layer 3 2 None tau_wall *
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (visc_t*([u]_1y)*deltayI)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 3 2 None tau_wall **************
!                                                           
!***********************************************************


qst(1-4+3,1-4+2,indvarsst(13)) =  (((1+&
                    param_float(21 + 5))/(((q(1-4+3,1-4+2,indvars(4))-&
                    0.5_wp*(q(1-4+3,1-4+2,indvars(2))*q(1-4+3,1-4+2,indvars(2))+&
                    q(1-4+3,1-4+2,indvars(3))*q(1-4+3,1-4+2,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-4+3,1-4+2,indvars(4))-&
                    0.5_wp*(q(1-4+3,1-4+2,indvars(2))*q(1-4+3,1-4+2,indvars(2))+&
                    q(1-4+3,1-4+2,indvars(3))*q(1-4+3,1-4+2,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(1-4+3,1-4+2,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-4+3,1-4+2,indvars(4))-&
                    0.5_wp*(q(1-4+3,1-4+2,indvars(2))*q(1-4+3,1-4+2,indvars(2))+&
                    q(1-4+3,1-4+2,indvars(3))*q(1-4+3,1-4+2,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-4+3,1-4+2,indvars(4))-&
                    0.5_wp*(q(1-4+3,1-4+2,indvars(2))*q(1-4+3,1-4+2,indvars(2))+&
                    q(1-4+3,1-4+2,indvars(3))*q(1-4+3,1-4+2,indvars(3)))))/param_float(4 + 5)**1.5*q(1-4+3,1-4+2,indvars(1)))**3/((q(1-4+3,1-4+2,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-4+3,1-4+2,indvars(4))-&
                    0.5_wp*(q(1-4+3,1-4+2,indvars(2))*q(1-4+3,1-4+2,indvars(2))+&
                    q(1-4+3,1-4+2,indvars(3))*q(1-4+3,1-4+2,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-4+3,1-4+2,indvars(4))-&
                    0.5_wp*(q(1-4+3,1-4+2,indvars(2))*q(1-4+3,1-4+2,indvars(2))+&
                    q(1-4+3,1-4+2,indvars(3))*q(1-4+3,1-4+2,indvars(3)))))/param_float(4 + 5)**1.5*q(1-4+3,1-4+2,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(1-4+3,1-4+2,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(1-4+3,1-4+2,indvars(4))-&
                    0.5_wp*(q(1-4+3,1-4+2,indvars(2))*q(1-4+3,1-4+2,indvars(2))+&
                    q(1-4+3,1-4+2,indvars(3))*q(1-4+3,1-4+2,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(1-4+3,1-4+2,indvars(4))-&
                    0.5_wp*(q(1-4+3,1-4+2,indvars(2))*q(1-4+3,1-4+2,indvars(2))+&
                    q(1-4+3,1-4+2,indvars(3))*q(1-4+3,1-4+2,indvars(3)))))/param_float(4 + 5)**1.5*q(1-4+3,1-4+2,indvars(1))))*param_float(1 + 5)*(d1_stemp_dy_0_1m4p31m4p2k)*qst(1-4+3,1-4+2,indvarsst(11)))


    enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO

 end subroutine boundarySchemestored_edges_i1_j1_3_2


