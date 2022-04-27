














!====================================================================================================
!
! General Boundary Conditions: enforce physical BC and compute Eqns with Boundary Scheme
!                              Both steady and unsteady BC can be generated (i.e. on vector q or rhs)
!     
!====================================================================================================

subroutine phyInlet_surfbc_faces_i1_8(param_float,hlo,neq,neqst,ind,idloop,idarray,nx,ny,nz,sizeblck,q,qst,q1,q2,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

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



 real(wp) ::  d1_rhs_rho_dx_0_1m4p0p0jk,d1_rhs_rho_dx_0_1m4p0p1jk,d1_rhs_rho_dx_0_1m4p0p2jk &
            ,d1_rhs_rho_dx_0_1m4p0jk &
            ,d1_rhs_rho_dx_1_1m4p0p0jk,d1_rhs_rho_dx_1_1m4p0p1jk,d1_rhs_rho_dx_1_1m4p0p2jk &
            ,d1_rhs_rho_dx_1_1m4p0jk &
            ,d1_rhs_rho_dy_0_1m4p0jm1k,d1_rhs_rho_dy_0_1m4p0jp0k,d1_rhs_rho_dy_0_1m4p0jp1k &
            ,d1_rhs_rho_dy_0_1m4p0jk 

  integer :: indvars(neq),indvarsst(neqst)
  
  indvars   = ind(1:neq)
  indvarsst = ind(1+neq:neq+neqst) 

!f2py intent(in)    :: q1,q2,qst,nx,ny,nz,rhs,nrk
!f2py intent(inout) :: q

      
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
! Start building layers for BC : i1 None None **************
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
! building source terms in RHS for layer 0 None None d(rho)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (1.0_wp/c**2*((((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI)*(gamma_m1)+0.5_wp*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI+((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI)))+[rho*v]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_rhs_rho_dx_0_1m4p0p0jk = q(1-4+0+0,j,indvars(2))

d1_rhs_rho_dx_0_1m4p0p1jk = q(1-4+0+1,j,indvars(2))

d1_rhs_rho_dx_0_1m4p0p2jk = q(1-4+0+2,j,indvars(2))

d1_rhs_rho_dx_0_1m4p0jk = -&
          1.5_wp*d1_rhs_rho_dx_0_1m4p0p0jk+&
          2.0_wp*d1_rhs_rho_dx_0_1m4p0p1jk-&
          0.5_wp*d1_rhs_rho_dx_0_1m4p0p2jk

d1_rhs_rho_dx_0_1m4p0jk = d1_rhs_rho_dx_0_1m4p0jk*param_float(1)

d1_rhs_rho_dx_1_1m4p0p0jk = (param_float(3 + 5))*q(1-4+0+0,j,indvars(1))*((q(1-4+0+0,j,indvars(4))-&
                    0.5_wp*(q(1-4+0+0,j,indvars(2))*q(1-4+0+0,j,indvars(2))+&
                    q(1-4+0+0,j,indvars(3))*q(1-4+0+0,j,indvars(3)))))

d1_rhs_rho_dx_1_1m4p0p1jk = (param_float(3 + 5))*q(1-4+0+1,j,indvars(1))*((q(1-4+0+1,j,indvars(4))-&
                    0.5_wp*(q(1-4+0+1,j,indvars(2))*q(1-4+0+1,j,indvars(2))+&
                    q(1-4+0+1,j,indvars(3))*q(1-4+0+1,j,indvars(3)))))

d1_rhs_rho_dx_1_1m4p0p2jk = (param_float(3 + 5))*q(1-4+0+2,j,indvars(1))*((q(1-4+0+2,j,indvars(4))-&
                    0.5_wp*(q(1-4+0+2,j,indvars(2))*q(1-4+0+2,j,indvars(2))+&
                    q(1-4+0+2,j,indvars(3))*q(1-4+0+2,j,indvars(3)))))

d1_rhs_rho_dx_1_1m4p0jk = -&
          1.5_wp*d1_rhs_rho_dx_1_1m4p0p0jk+&
          2.0_wp*d1_rhs_rho_dx_1_1m4p0p1jk-&
          0.5_wp*d1_rhs_rho_dx_1_1m4p0p2jk

d1_rhs_rho_dx_1_1m4p0jk = d1_rhs_rho_dx_1_1m4p0jk*param_float(1)

d1_rhs_rho_dy_0_1m4p0jm1k = q(1-4+0,j-1,indvars(1))*q(1-4+0,j-1,indvars(3))

d1_rhs_rho_dy_0_1m4p0jp1k = q(1-4+0,j+1,indvars(1))*q(1-4+0,j+1,indvars(3))

d1_rhs_rho_dy_0_1m4p0jk = -&
          0.5_wp*d1_rhs_rho_dy_0_1m4p0jm1k+&
          0.5_wp*d1_rhs_rho_dy_0_1m4p0jp1k

d1_rhs_rho_dy_0_1m4p0jk = d1_rhs_rho_dy_0_1m4p0jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None d(rho)/dt **********
!                                                           
!***********************************************************


rhs(1-4+0,j,indvars(1)) =   -  ( (1.0_wp/(param_float(23 + 5)*(param_float(3 + 5))*q(1-4+0,j,indvars(1))*((q(1-4+0,j,indvars(4))-&
                    0.5_wp*(q(1-4+0,j,indvars(2))*q(1-4+0,j,indvars(2))+&
                    q(1-4+0,j,indvars(3))*q(1-4+0,j,indvars(3)))))/q(1-4+0,j,indvars(1)))**0.5**2*(((((param_float(23 + 5)*(param_float(3 + 5))*q(1-4+0,j,indvars(1))*((q(1-4+0,j,indvars(4))-&
                    0.5_wp*(q(1-4+0,j,indvars(2))*q(1-4+0,j,indvars(2))+&
                    q(1-4+0,j,indvars(3))*q(1-4+0,j,indvars(3)))))/q(1-4+0,j,indvars(1)))**0.5+&
                    q(1-4+0,j,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(1-4+0,j,indvars(1))*((q(1-4+0,j,indvars(4))-&
                    0.5_wp*(q(1-4+0,j,indvars(2))*q(1-4+0,j,indvars(2))+&
                    q(1-4+0,j,indvars(3))*q(1-4+0,j,indvars(3)))))/q(1-4+0,j,indvars(1)))**0.5*q(1-4+0,j,indvars(1))*d1_rhs_rho_dx_0_1m4p0jk+&
                    1.0_wp*d1_rhs_rho_dx_1_1m4p0jk))*qst(1-4+0,j,indvarsst(10)))*(param_float(3 + 5))+&
                    0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(1-4+0,j,indvars(1))*((q(1-4+0,j,indvars(4))-&
                    0.5_wp*(q(1-4+0,j,indvars(2))*q(1-4+0,j,indvars(2))+&
                    q(1-4+0,j,indvars(3))*q(1-4+0,j,indvars(3)))))/q(1-4+0,j,indvars(1)))**0.5+&
                    q(1-4+0,j,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(1-4+0,j,indvars(1))*((q(1-4+0,j,indvars(4))-&
                    0.5_wp*(q(1-4+0,j,indvars(2))*q(1-4+0,j,indvars(2))+&
                    q(1-4+0,j,indvars(3))*q(1-4+0,j,indvars(3)))))/q(1-4+0,j,indvars(1)))**0.5*q(1-4+0,j,indvars(1))*d1_rhs_rho_dx_0_1m4p0jk+&
                    1.0_wp*d1_rhs_rho_dx_1_1m4p0jk))*qst(1-4+0,j,indvarsst(10))+&
                    (((param_float(23 + 5)*(param_float(3 + 5))*q(1-4+0,j,indvars(1))*((q(1-4+0,j,indvars(4))-&
                    0.5_wp*(q(1-4+0,j,indvars(2))*q(1-4+0,j,indvars(2))+&
                    q(1-4+0,j,indvars(3))*q(1-4+0,j,indvars(3)))))/q(1-4+0,j,indvars(1)))**0.5+&
                    q(1-4+0,j,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(1-4+0,j,indvars(1))*((q(1-4+0,j,indvars(4))-&
                    0.5_wp*(q(1-4+0,j,indvars(2))*q(1-4+0,j,indvars(2))+&
                    q(1-4+0,j,indvars(3))*q(1-4+0,j,indvars(3)))))/q(1-4+0,j,indvars(1)))**0.5*q(1-4+0,j,indvars(1))*d1_rhs_rho_dx_0_1m4p0jk+&
                    1.0_wp*d1_rhs_rho_dx_1_1m4p0jk))*qst(1-4+0,j,indvarsst(10)))))+d1_rhs_rho_dy_0_1m4p0jk ) 

     enddo

    enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO

 end subroutine phyInlet_surfbc_faces_i1_8


