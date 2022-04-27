














!====================================================================================================
!
! General Boundary Conditions: compute Eqns with Boundary Scheme
!     
!====================================================================================================

subroutine boundarySchemestored_edges_i1_j1_2_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

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



 real(wp) ::  d1_stemp_dx_0_1m4p2m11m4p0k,d1_stemp_dx_0_1m4p2p01m4p0k,d1_stemp_dx_0_1m4p2p11m4p0k &
            ,d1_stemp_dx_0_1m4p21m4p0k &
            ,d1_stemp_dy_0_1m4p21m4p0p0k,d1_stemp_dy_0_1m4p21m4p0p1k,d1_stemp_dy_0_1m4p21m4p0p2k &
            ,d1_stemp_dy_0_1m4p21m4p0k 

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
! BC layer: 2 0 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 2 0 None stemp ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (2.0_wp*(dabs(0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_1m4p2m11m4p0k = q(1-4+2-1,1-4+0,indvars(3))

d1_stemp_dx_0_1m4p2p11m4p0k = q(1-4+2+1,1-4+0,indvars(3))

d1_stemp_dx_0_1m4p21m4p0k = -&
          0.5_wp*d1_stemp_dx_0_1m4p2m11m4p0k+&
          0.5_wp*d1_stemp_dx_0_1m4p2p11m4p0k

d1_stemp_dx_0_1m4p21m4p0k = d1_stemp_dx_0_1m4p21m4p0k*param_float(1)

d1_stemp_dy_0_1m4p21m4p0p0k = q(1-4+2,1-4+0+0,indvars(2))

d1_stemp_dy_0_1m4p21m4p0p1k = q(1-4+2,1-4+0+1,indvars(2))

d1_stemp_dy_0_1m4p21m4p0p2k = q(1-4+2,1-4+0+2,indvars(2))

d1_stemp_dy_0_1m4p21m4p0k = -&
          1.5_wp*d1_stemp_dy_0_1m4p21m4p0p0k+&
          2.0_wp*d1_stemp_dy_0_1m4p21m4p0p1k-&
          0.5_wp*d1_stemp_dy_0_1m4p21m4p0p2k

d1_stemp_dy_0_1m4p21m4p0k = d1_stemp_dy_0_1m4p21m4p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 2 0 None stemp *****************
!                                                           
!***********************************************************


qst(1-4+2,1-4+0,indvarsst(4)) =  (2.0_wp*(dabs(0.5_wp*(qst(1-4+2,1-4+0,indvarsst(11))*(d1_stemp_dy_0_1m4p21m4p0k)-&
                    qst(1-4+2,1-4+0,indvarsst(10))*(d1_stemp_dx_0_1m4p21m4p0k)))))



!***********************************************************
!                                                           
! building source terms in RHS for layer 2 0 None SS *******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (stemp+fv2*ReI*nut/(k**2.0_wp*eta**2.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 2 0 None SS ********************
!                                                           
!***********************************************************


qst(1-4+2,1-4+0,indvarsst(12)) =  (qst(1-4+2,1-4+0,indvarsst(4))+&
                    (1.0_wp-&
                    (q(1-4+2,1-4+0,indvars(5))/1.0_wp*q(1-4+2,1-4+0,indvars(1)))/(1.0_wp+&
                    (q(1-4+2,1-4+0,indvars(5))/1.0_wp*q(1-4+2,1-4+0,indvars(1)))*((q(1-4+2,1-4+0,indvars(5))/1.0_wp*q(1-4+2,1-4+0,indvars(1)))**3.0_wp/((q(1-4+2,1-4+0,indvars(5))/1.0_wp*q(1-4+2,1-4+0,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))))*param_float(1 + 5)*q(1-4+2,1-4+0,indvars(5))/(param_float(9 + 5)**2.0_wp*qst(1-4+2,1-4+0,indvarsst(2))**2.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 2 0 None tau_wall *
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (visc_t*([u]_1y)*deltayI)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 2 0 None tau_wall **************
!                                                           
!***********************************************************


qst(1-4+2,1-4+0,indvarsst(13)) =  ((1.0_wp)*(1.0_wp+&
                    ((q(1-4+2,1-4+0,indvars(5))/1.0_wp*q(1-4+2,1-4+0,indvars(1)))**3.0_wp/((q(1-4+2,1-4+0,indvars(5))/1.0_wp*q(1-4+2,1-4+0,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-4+2,1-4+0,indvars(5))/1.0_wp*q(1-4+2,1-4+0,indvars(1))))*param_float(1 + 5)*(d1_stemp_dy_0_1m4p21m4p0k)*qst(1-4+2,1-4+0,indvarsst(11)))



!***********************************************************
!                                                           
! building source terms in RHS for layer 2 0 None visc_SA **
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! nut*rho
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 2 0 None visc_SA ***************
!                                                           
!***********************************************************


qst(1-4+2,1-4+0,indvarsst(14)) =  q(1-4+2,1-4+0,indvars(5))*q(1-4+2,1-4+0,indvars(1))



!***********************************************************
!                                                           
! building source terms in RHS for layer 2 0 None visc_turb 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! nut*rho*fv1
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 2 0 None visc_turb *************
!                                                           
!***********************************************************


qst(1-4+2,1-4+0,indvarsst(15)) =  q(1-4+2,1-4+0,indvars(5))*q(1-4+2,1-4+0,indvars(1))*((q(1-4+2,1-4+0,indvars(5))/1.0_wp*q(1-4+2,1-4+0,indvars(1)))**3.0_wp/((q(1-4+2,1-4+0,indvars(5))/1.0_wp*q(1-4+2,1-4+0,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 2 0 None Pressure *
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (gamma_m1)*rho*(e)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 2 0 None Pressure **************
!                                                           
!***********************************************************


qst(1-4+2,1-4+0,indvarsst(16)) =  (param_float(3 + 5))*q(1-4+2,1-4+0,indvars(1))*((q(1-4+2,1-4+0,indvars(4))-&
                    0.5_wp*(q(1-4+2,1-4+0,indvars(2))*q(1-4+2,1-4+0,indvars(2))+&
                    q(1-4+2,1-4+0,indvars(3))*q(1-4+2,1-4+0,indvars(3)))))


    enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO

 end subroutine boundarySchemestored_edges_i1_j1_2_0


