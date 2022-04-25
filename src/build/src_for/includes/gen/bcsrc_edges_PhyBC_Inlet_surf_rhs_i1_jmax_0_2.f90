














!====================================================================================================
!
! General Boundary Conditions: compute Eqns with Boundary Scheme
!     
!====================================================================================================

subroutine boundaryScheme_edges_PhyBC_Inlet_surf_rhs_i1_jmax_0_2(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

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



 real(wp) ::  d1_rhs_rho_dx_0_1m4p0p0nyp4m2k,d1_rhs_rho_dx_0_1m4p0p1nyp4m2k,d1_rhs_rho_dx_0_1m4p0p2nyp4m2k &
            ,d1_rhs_rho_dx_0_1m4p0nyp4m2k &
            ,d1_rhs_rho_dx_1_1m4p0p0nyp4m2k,d1_rhs_rho_dx_1_1m4p0p1nyp4m2k,d1_rhs_rho_dx_1_1m4p0p2nyp4m2k &
            ,d1_rhs_rho_dx_1_1m4p0nyp4m2k &
            ,d1_rhs_rho_dy_0_1m4p0nyp4m2m1k,d1_rhs_rho_dy_0_1m4p0nyp4m2p0k,d1_rhs_rho_dy_0_1m4p0nyp4m2p1k &
            ,d1_rhs_rho_dy_0_1m4p0nyp4m2k 

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
! BC layer: 0 2 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 0 2 None d(rho)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (1.0_wp/c**2*((((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI)*(gamma_m1)+0.5_wp*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI+((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI)))+[rho*v]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_rhs_rho_dx_0_1m4p0p0nyp4m2k = q(1-4+0+0,ny+4-2,indvars(2))

d1_rhs_rho_dx_0_1m4p0p1nyp4m2k = q(1-4+0+1,ny+4-2,indvars(2))

d1_rhs_rho_dx_0_1m4p0p2nyp4m2k = q(1-4+0+2,ny+4-2,indvars(2))

d1_rhs_rho_dx_0_1m4p0nyp4m2k = -&
          1.5_wp*d1_rhs_rho_dx_0_1m4p0p0nyp4m2k+&
          2.0_wp*d1_rhs_rho_dx_0_1m4p0p1nyp4m2k-&
          0.5_wp*d1_rhs_rho_dx_0_1m4p0p2nyp4m2k

d1_rhs_rho_dx_0_1m4p0nyp4m2k = d1_rhs_rho_dx_0_1m4p0nyp4m2k*param_float(1)

d1_rhs_rho_dx_1_1m4p0p0nyp4m2k = (param_float(3 + 5))*q(1-4+0+0,ny+4-2,indvars(1))*((q(1-4+0+0,ny+4-2,indvars(4))-&
                    0.5_wp*(q(1-4+0+0,ny+4-2,indvars(2))*q(1-4+0+0,ny+4-2,indvars(2))+&
                    q(1-4+0+0,ny+4-2,indvars(3))*q(1-4+0+0,ny+4-2,indvars(3)))))

d1_rhs_rho_dx_1_1m4p0p1nyp4m2k = (param_float(3 + 5))*q(1-4+0+1,ny+4-2,indvars(1))*((q(1-4+0+1,ny+4-2,indvars(4))-&
                    0.5_wp*(q(1-4+0+1,ny+4-2,indvars(2))*q(1-4+0+1,ny+4-2,indvars(2))+&
                    q(1-4+0+1,ny+4-2,indvars(3))*q(1-4+0+1,ny+4-2,indvars(3)))))

d1_rhs_rho_dx_1_1m4p0p2nyp4m2k = (param_float(3 + 5))*q(1-4+0+2,ny+4-2,indvars(1))*((q(1-4+0+2,ny+4-2,indvars(4))-&
                    0.5_wp*(q(1-4+0+2,ny+4-2,indvars(2))*q(1-4+0+2,ny+4-2,indvars(2))+&
                    q(1-4+0+2,ny+4-2,indvars(3))*q(1-4+0+2,ny+4-2,indvars(3)))))

d1_rhs_rho_dx_1_1m4p0nyp4m2k = -&
          1.5_wp*d1_rhs_rho_dx_1_1m4p0p0nyp4m2k+&
          2.0_wp*d1_rhs_rho_dx_1_1m4p0p1nyp4m2k-&
          0.5_wp*d1_rhs_rho_dx_1_1m4p0p2nyp4m2k

d1_rhs_rho_dx_1_1m4p0nyp4m2k = d1_rhs_rho_dx_1_1m4p0nyp4m2k*param_float(1)

d1_rhs_rho_dy_0_1m4p0nyp4m2m1k = q(1-4+0,ny+4-2-1,indvars(1))*q(1-4+0,ny+4-2-1,indvars(3))

d1_rhs_rho_dy_0_1m4p0nyp4m2p1k = q(1-4+0,ny+4-2+1,indvars(1))*q(1-4+0,ny+4-2+1,indvars(3))

d1_rhs_rho_dy_0_1m4p0nyp4m2k = -&
          0.5_wp*d1_rhs_rho_dy_0_1m4p0nyp4m2m1k+&
          0.5_wp*d1_rhs_rho_dy_0_1m4p0nyp4m2p1k

d1_rhs_rho_dy_0_1m4p0nyp4m2k = d1_rhs_rho_dy_0_1m4p0nyp4m2k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 2 None d(rho)/dt *************
!                                                           
!***********************************************************


rhs(1-4+0,ny+4-2,indvars(1)) =   -  ( (1.0_wp/(param_float(23 + 5)*(param_float(3 + 5))*q(1-4+0,ny+4-2,indvars(1))*((q(1-4+0,ny+4-2,indvars(4))-&
                    0.5_wp*(q(1-4+0,ny+4-2,indvars(2))*q(1-4+0,ny+4-2,indvars(2))+&
                    q(1-4+0,ny+4-2,indvars(3))*q(1-4+0,ny+4-2,indvars(3)))))/q(1-4+0,ny+4-2,indvars(1)))**2*(((((param_float(23 + 5)*(param_float(3 + 5))*q(1-4+0,ny+4-2,indvars(1))*((q(1-4+0,ny+4-2,indvars(4))-&
                    0.5_wp*(q(1-4+0,ny+4-2,indvars(2))*q(1-4+0,ny+4-2,indvars(2))+&
                    q(1-4+0,ny+4-2,indvars(3))*q(1-4+0,ny+4-2,indvars(3)))))/q(1-4+0,ny+4-2,indvars(1)))+&
                    q(1-4+0,ny+4-2,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(1-4+0,ny+4-2,indvars(1))*((q(1-4+0,ny+4-2,indvars(4))-&
                    0.5_wp*(q(1-4+0,ny+4-2,indvars(2))*q(1-4+0,ny+4-2,indvars(2))+&
                    q(1-4+0,ny+4-2,indvars(3))*q(1-4+0,ny+4-2,indvars(3)))))/q(1-4+0,ny+4-2,indvars(1)))*q(1-4+0,ny+4-2,indvars(1))*d1_rhs_rho_dx_0_1m4p0nyp4m2k+&
                    1.0_wp*d1_rhs_rho_dx_1_1m4p0nyp4m2k))*qst(1-4+0,ny+4-2,indvarsst(10)))*(param_float(3 + 5))+&
                    0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(1-4+0,ny+4-2,indvars(1))*((q(1-4+0,ny+4-2,indvars(4))-&
                    0.5_wp*(q(1-4+0,ny+4-2,indvars(2))*q(1-4+0,ny+4-2,indvars(2))+&
                    q(1-4+0,ny+4-2,indvars(3))*q(1-4+0,ny+4-2,indvars(3)))))/q(1-4+0,ny+4-2,indvars(1)))+&
                    q(1-4+0,ny+4-2,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(1-4+0,ny+4-2,indvars(1))*((q(1-4+0,ny+4-2,indvars(4))-&
                    0.5_wp*(q(1-4+0,ny+4-2,indvars(2))*q(1-4+0,ny+4-2,indvars(2))+&
                    q(1-4+0,ny+4-2,indvars(3))*q(1-4+0,ny+4-2,indvars(3)))))/q(1-4+0,ny+4-2,indvars(1)))*q(1-4+0,ny+4-2,indvars(1))*d1_rhs_rho_dx_0_1m4p0nyp4m2k+&
                    1.0_wp*d1_rhs_rho_dx_1_1m4p0nyp4m2k))*qst(1-4+0,ny+4-2,indvarsst(10))+&
                    (((param_float(23 + 5)*(param_float(3 + 5))*q(1-4+0,ny+4-2,indvars(1))*((q(1-4+0,ny+4-2,indvars(4))-&
                    0.5_wp*(q(1-4+0,ny+4-2,indvars(2))*q(1-4+0,ny+4-2,indvars(2))+&
                    q(1-4+0,ny+4-2,indvars(3))*q(1-4+0,ny+4-2,indvars(3)))))/q(1-4+0,ny+4-2,indvars(1)))+&
                    q(1-4+0,ny+4-2,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(1-4+0,ny+4-2,indvars(1))*((q(1-4+0,ny+4-2,indvars(4))-&
                    0.5_wp*(q(1-4+0,ny+4-2,indvars(2))*q(1-4+0,ny+4-2,indvars(2))+&
                    q(1-4+0,ny+4-2,indvars(3))*q(1-4+0,ny+4-2,indvars(3)))))/q(1-4+0,ny+4-2,indvars(1)))*q(1-4+0,ny+4-2,indvars(1))*d1_rhs_rho_dx_0_1m4p0nyp4m2k+&
                    1.0_wp*d1_rhs_rho_dx_1_1m4p0nyp4m2k))*qst(1-4+0,ny+4-2,indvarsst(10)))))+d1_rhs_rho_dy_0_1m4p0nyp4m2k ) 


    enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO

 end subroutine boundaryScheme_edges_PhyBC_Inlet_surf_rhs_i1_jmax_0_2


