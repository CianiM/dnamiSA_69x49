













  subroutine time_march(param_int,param_float,data_float)

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
  real(wp), intent(inout) ::  data_float(*)
  integer,  intent(in)    ::   param_int(*)

  ! LOCAL VARIABLES

  integer :: shiftstore
  integer :: shiftfacei,shiftfacej,shiftfacek
  integer :: shiftedgeij,shiftedgejk,shiftedgeik
  integer :: nbc
  integer :: nfacei,nfacej,nfacek
  integer :: nedgeij,nedgejk
  integer :: nvar_f(3),nvar_e(3)

  nbc       = param_int(13 + param_int(5) + param_int(6)  )
  
  nvar_f(1) = param_int(13 + param_int(5) + param_int(6) + 1 + nbc      )
  nvar_f(2) = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 1  )
  nvar_f(3) = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 2  ) 

  nvar_e(1) = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 3      )
  nvar_e(2) = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 3 + 1  )
  nvar_e(3) = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 3 + 2  ) 
   
  nfacei    = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 6  )
  nfacej    = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 7  )
  nfacek    = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 8  )
  nedgeij   = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 9  )
  nedgejk   = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 10 )

  if ( param_int(6) .gt. 0 ) then
   shiftstore = param_int(6)*param_int(9)
  else
   shiftstore = 1
  endif

  shiftfacei  = shiftstore  + nfacei
  shiftfacej  = shiftfacei  + nfacej
  shiftfacek  = shiftfacej  + nfacek
  shiftedgeij = shiftfacek  + nedgeij
  shiftedgejk = shiftedgeij + nedgejk

  call rk3_stepk(param_float                                                          ,&
                 param_int(1), param_int(8)                               ,&
                 param_int(5)                                                  ,&
                 param_int(6)                                                ,&
                 param_int(13)                                                  ,&
                 param_int(2),param_int(3) , param_int(4)              ,&
                 param_int(10)                                             ,&
                 param_int(13 + param_int(5) + param_int(6))    ,&
                 param_int(13 + param_int(5) + param_int(6) + 1),&
                 nvar_f(1),nvar_f(2),nvar_f(3),&
                 nvar_e(1),nvar_e(2),nvar_e(3),&
                 data_float(1)                                                        ,& ! q
                 data_float(1+param_int(7)  )                               ,& ! q1
                 data_float(1+param_int(7)*2)                               ,& ! q2 
                 data_float(1+param_int(7)*3)                               ,& ! rhs
                 data_float(1+param_int(7)*4)                               ,& ! stored     (if any)
                 data_float(1+param_int(7)*4 + shiftstore)                  ,& ! qbcface_i  (if any)
                 data_float(1+param_int(7)*4 + shiftfacei)                  ,& ! qbcface_j  (if any)
                 data_float(1+param_int(7)*4 + shiftfacej)                  ,& ! qbcface_k  (if any)
                 data_float(1+param_int(7)*4 + shiftfacek)                  ,& ! qbcedge_ij (if any)
                 data_float(1+param_int(7)*4 + shiftedgeij)                 ,& ! qbcedge_jk (if any)
                 data_float(1+param_int(7)*4 + shiftedgejk)                 )  ! qbcedge_ik (if any)
      

  end subroutine time_march


  subroutine stored(param_int,param_float,data_float,type_st)

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
  real(wp), intent(inout) ::  data_float(*)
  integer,  intent(in)    ::   param_int(*)
  integer                 ::   type_st
 
 !f2py integer optional, intent(in)    :: type_st = 0

  integer :: shiftstore
  integer :: shiftfacei,shiftfacej,shiftfacek
  integer :: shiftedgeij,shiftedgejk,shiftedgeik
  integer :: nbc
  integer :: nfacei,nfacej,nfacek
  integer :: nedgeij,nedgejk
  integer :: nvar_f(3),nvar_e(3)

  nbc       = param_int(13 + param_int(5) + param_int(6)  )
  
  nvar_f(1) = param_int(13 + param_int(5) + param_int(6) + 1 + nbc      )
  nvar_f(2) = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 1  )
  nvar_f(3) = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 2  ) 

  nvar_e(1) = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 3      )
  nvar_e(2) = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 3 + 1  )
  nvar_e(3) = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 3 + 2  ) 
   
  nfacei    = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 6  )
  nfacej    = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 7  )
  nfacek    = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 8  )
  nedgeij   = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 9  )
  nedgejk   = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 10 )

  if ( param_int(6) .gt. 0 ) then
   shiftstore = param_int(6)*param_int(9)
  else
   shiftstore = 1
  endif

  shiftfacei  = shiftstore  + nfacei
  shiftfacej  = shiftfacei  + nfacej
  shiftfacek  = shiftfacej  + nfacek
  shiftedgeij = shiftfacek  + nedgeij
  shiftedgejk = shiftedgeij + nedgejk

  call cmp_stored(param_float                                                          ,&
                  type_st                                                              ,&
                  param_int(1), param_int(8)                               ,&
                  param_int(5)                                                  ,&
                  param_int(6)                                                ,&
                  param_int(13)                                                  ,&
                  param_int(2),param_int(3) , param_int(4)              ,&
                  param_int(10)                                             ,&
                  param_int(13 + param_int(5) + param_int(6))    ,&
                  param_int(13 + param_int(5) + param_int(6) + 1),&                  
                  nvar_f(1),nvar_f(2),nvar_f(3),&
                  nvar_e(1),nvar_e(2),nvar_e(3),&
                  data_float(1)                                                        ,& ! q
                  data_float(1+param_int(7)*3)                               ,& ! rhs
                  data_float(1+param_int(7)*4)                               ,& ! stored (if any)
                  data_float(1+param_int(7)*4 + shiftstore)                  ,& ! qbcface_i  (if any)
                  data_float(1+param_int(7)*4 + shiftfacei)                  ,& ! qbcface_j  (if any)
                  data_float(1+param_int(7)*4 + shiftfacej)                  ,& ! qbcface_k  (if any)
                  data_float(1+param_int(7)*4 + shiftfacek)                  ,& ! qbcedge_ij (if any)
                  data_float(1+param_int(7)*4 + shiftedgeij)                 ,& ! qbcedge_jk (if any)
                  data_float(1+param_int(7)*4 + shiftedgejk)                 )  ! qbcedge_ik (if any)         


  end subroutine stored  

  subroutine cmp_stored(param_float,type_st,hlo,nrk,neq,neqst,ind,nx,ny,nz,sizeblck,&
                        nbc,bc,&
                        nfacei,nfacej,nfacek,&
                        nedgeij,nedgejk,nedgeik,&
                        q,rhs,qst,& 
                        qface_i , qface_j ,qface_k,&
                        qedge_ij, qedge_jk,qedge_ik)

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


  
  integer :: i,j,k
  integer :: bi,bj,bk
  integer :: biend,bjend,bkend
  integer :: size_bi,size_bj,size_bk
  integer :: ibc

  real(wp), intent(in)    :: param_float(*)

  integer, intent(in) :: nx,ny,nz,hlo,nrk,neq,neqst,type_st
  integer, intent(in) :: nfacei,nfacej,nfacek,nedgeij,nedgejk,nedgeik
  integer, intent(in) :: ind(1:neq+neqst)
  integer, intent(in) :: sizeblck(3)
  integer, intent(in) :: bc(nbc),nbc

  
real(wp),intent(inout) :: q(1-hlo:nx+hlo,1-hlo:ny+hlo,neq),&
                        rhs(1-hlo:nx+hlo,1-hlo:ny+hlo,neq),&
                        qst(1-hlo:nx+hlo,1-hlo:ny+hlo,neqst)

real(wp),intent(inout) :: qface_i(1-hlo:ny+hlo,nfacei),&
                       qface_j(1-hlo:nx+hlo,nfacej),&
                       qface_k(1),&
                       qedge_ij(1),&
                       qedge_jk(1),&
                       qedge_ik(1)

  
  integer :: idloop(6),idarray(6),indvars(neq),indvarsst(neqst)
  integer :: nvar_f(3),nvar_e(3)

  indvars   = ind(1:neq)
  indvarsst = ind(1+neq:neq+neqst) 

!f2py intent(in)    :: nx,ny,nz,nrk,q,nvar_f,nvar_e,neq,nbc
!f2py intent(inout) :: qst

nvar_f(1) = nfacei
nvar_f(2) = nfacej
nvar_f(3) = nfacek

nvar_e(1) = nedgeij
nvar_e(2) = nedgejk
nvar_e(3) = nedgeik

!$OMP PARALLEL DEFAULT(SHARED) private(idloop)

      
  idarray(1) = 1 -hlo
  idarray(2) = nx+hlo

  idarray(3) = 1 -hlo
  idarray(4) = ny+hlo

  idarray(5) = 1 -hlo
  idarray(6) = nz+hlo

size_bk = sizeblck(3)
size_bj = sizeblck(2)
size_bi = sizeblck(1)


    SELECT CASE (type_st)

    CASE(0)

!$OMP DO SCHEDULE(GUIDED,4) COLLAPSE(2)  
   do bk=1,nz,size_bk  
    do bj=1,ny,size_bj 
      do bi=1,nx,size_bi 
    
idloop(6) = min( bk+size_bk, nz+1)-1
idloop(4) = min( bj+size_bj, ny+1)-1
idloop(2) = min( bi+size_bi, nx+1)-1


idloop(5) = bk
idloop(3) = bj
idloop(1) = bi 

call cmpstored(param_float,ind,idloop,idarray,neq,neqst,sizeblck,q,qst,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

    enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO NOWAIT



do ibc=1,nbc

    SELECT CASE (bc(ibc))

CASE (1)
      call  boundarySchemestored_edges_i1_j1_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemestored_edges_i1_j1_1(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemestored_edges_i1_j1_2(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemestored_edges_i1_j1_3(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (2)
      call  boundarySchemestored_edges_i1_jmax_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemestored_edges_i1_jmax_1(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemestored_edges_i1_jmax_2(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemestored_edges_i1_jmax_3(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (3)
      call  boundarySchemestored_edges_imax_j1_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemestored_edges_imax_j1_1(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemestored_edges_imax_j1_2(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemestored_edges_imax_j1_3(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (5)
      call  boundarySchemestored_faces_i1(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (6)
      call  boundarySchemestored_faces_imax(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (7)
      call  boundarySchemestored_faces_j1(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (8)
      call  boundarySchemestored_faces_jmax(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (4)
      call  boundarySchemestored_edges_imax_jmax_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemestored_edges_imax_jmax_1(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemestored_edges_imax_jmax_2(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemestored_edges_imax_jmax_3(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

    END SELECT

enddo



do ibc=1,nbc

    SELECT CASE (bc(ibc))

CASE (3)
      call  boundarySchemevarbc_edges_imax_j1_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemevarbc_edges_imax_j1_1(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemevarbc_edges_imax_j1_2(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemevarbc_edges_imax_j1_3(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (6)
      call  boundarySchemevarbc_faces_imax(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (8)
      call  boundarySchemevarbc_faces_jmax(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (4)
      call  boundarySchemevarbc_edges_imax_jmax_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemevarbc_edges_imax_jmax_1(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemevarbc_edges_imax_jmax_2(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemevarbc_edges_imax_jmax_3(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

    END SELECT

enddo

    
    CASE(1)

!$OMP DO SCHEDULE(GUIDED,4) COLLAPSE(2)  
   do bk=1,nz,size_bk  
    do bj=1,ny,size_bj 
      do bi=1,nx,size_bi 
    
idloop(6) = min( bk+size_bk, nz+1)-1
idloop(4) = min( bj+size_bj, ny+1)-1
idloop(2) = min( bi+size_bi, nx+1)-1


idloop(5) = bk
idloop(3) = bj
idloop(1) = bi 

call cmpstoredstatic(param_float,ind,idloop,idarray,neq,neqst,sizeblck,q,qst,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

    enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO NOWAIT



do ibc=1,nbc

    SELECT CASE (bc(ibc))

CASE (1)
      call  boundarySchemestoredstatic_edges_i1_j1_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemestoredstatic_edges_i1_j1_1(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemestoredstatic_edges_i1_j1_2(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemestoredstatic_edges_i1_j1_3(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (2)
      call  boundarySchemestoredstatic_edges_i1_jmax_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemestoredstatic_edges_i1_jmax_1(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemestoredstatic_edges_i1_jmax_2(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemestoredstatic_edges_i1_jmax_3(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (3)
      call  boundarySchemestoredstatic_edges_imax_j1_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemestoredstatic_edges_imax_j1_1(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemestoredstatic_edges_imax_j1_2(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemestoredstatic_edges_imax_j1_3(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (5)
      call  boundarySchemestoredstatic_faces_i1(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (6)
      call  boundarySchemestoredstatic_faces_imax(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (7)
      call  boundarySchemestoredstatic_faces_j1(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (8)
      call  boundarySchemestoredstatic_faces_jmax(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (4)
      call  boundarySchemestoredstatic_edges_imax_jmax_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemestoredstatic_edges_imax_jmax_1(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemestoredstatic_edges_imax_jmax_2(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundarySchemestoredstatic_edges_imax_jmax_3(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

    END SELECT

enddo



do ibc=1,nbc

    SELECT CASE (bc(ibc))


    END SELECT

enddo

  END SELECT

!$OMP END PARALLEL

 end subroutine cmp_stored


 subroutine rk3_stepk(param_float,hlo,nrk,neq,neqst,ind,nx,ny,nz,sizeblck,&
                      nbc,bc,&
                      nfacei,nfacej,nfacek,&
                      nedgeij,nedgejk,nedgeik,&
                      q,q1,q2,rhs,qst,& 
                      qface_i , qface_j ,qface_k,&
                      qedge_ij, qedge_jk,qedge_ik)

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


  
  integer :: i,j,k
  integer :: bi,bj,bk
  integer :: biend,bjend,bkend
  integer :: size_bi,size_bj,size_bk
  integer :: ibc

  real(wp), intent(in)    :: param_float(*)

  integer, intent(in) :: nx,ny,nz,hlo,nrk,neq,neqst
  integer, intent(in) :: nfacei,nfacej,nfacek,nedgeij,nedgejk,nedgeik
  integer, intent(in) :: ind(1:neq+neqst)
  integer, intent(in) :: sizeblck(3)
  

  integer, intent(in) :: bc(nbc),nbc

  real(wp),  dimension(1:3) :: rk1 = (/ 2.0_wp/ 3.0_wp, &
                                        5.0_wp/12.0_wp, &
                                        3.0_wp/ 5.0_wp /)
  real(wp),  dimension(1:3) :: rk2 = (/ 1.0_wp/ 4.0_wp, &
                                        3.0_wp/20.0_wp, &
                                        3.0_wp/ 5.0_wp /)

  real(wp) :: one_over_rho
  
real(wp),intent(inout) ::  q(1-hlo:nx+hlo,1-hlo:ny+hlo,neq),&
                          q1(1-hlo:nx+hlo,1-hlo:ny+hlo,neq),&
                          q2(1-hlo:nx+hlo,1-hlo:ny+hlo,neq),&
                         rhs(1-hlo:nx+hlo,1-hlo:ny+hlo,neq),&
                        qst(1-hlo:nx+hlo,1-hlo:ny+hlo,neqst)

real(wp),intent(inout) :: qface_i(1-hlo:ny+hlo,nfacei),&
                       qface_j(1-hlo:nx+hlo,nfacej),&
                       qface_k(1),&
                       qedge_ij(1),&
                       qedge_jk(1),&
                       qedge_ik(1)
  
  integer :: idloop(6),idarray(6),idrhs(6),indvars(neq),indvarsst(neqst)
  integer :: nvar_f(3),nvar_e(3)



  indvars   = ind(1:neq)
  indvarsst = ind(1+neq:neq+neqst) 

!f2py intent(in)    :: q1,q2,nx,ny,nz,rhs,nrk,qst,nvar_f,nvar_e,neq,nbc
!f2py intent(inout) :: q

nvar_f(1) = nfacei
nvar_f(2) = nfacej
nvar_f(3) = nfacek

nvar_e(1) = nedgeij
nvar_e(2) = nedgejk
nvar_e(3) = nedgeik



!$OMP PARALLEL DEFAULT(SHARED) private(idloop)


      
  idrhs(1) = 1 
  idrhs(2) = nx

  idrhs(3) = 1 
  idrhs(4) = ny

  idrhs(5) = 1 
  idrhs(6) = nz
      
  idarray(1) = 1 -hlo
  idarray(2) = nx+hlo

  idarray(3) = 1 -hlo
  idarray(4) = ny+hlo

  idarray(5) = 1 -hlo
  idarray(6) = nz+hlo

size_bk = sizeblck(3)
size_bj = sizeblck(2)
size_bi = sizeblck(1)


! #include "include_bcq.f90"


!$OMP DO SCHEDULE(GUIDED,4) COLLAPSE(2)  
  do bk=1,nz,size_bk  
    do bj=1,ny,size_bj 
      do bi=1,nx,size_bi 
    
idloop(6) = min( bk+size_bk, nz+1)-1
idloop(4) = min( bj+size_bj, ny+1)-1
idloop(2) = min( bi+size_bi, nx+1)-1

idloop(5) = bk
idloop(3) = bj
idloop(1) = bi 

call cmprhs(param_float,ind,idloop,idarray,neq,neqst,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

    enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO



do ibc=1,nbc

    SELECT CASE (bc(ibc))

CASE (1)
      call  boundaryScheme_edges_i1_j1_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_edges_i1_j1_1(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_edges_i1_j1_2(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_edges_i1_j1_3(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      idrhs(1)=idarray(1)
      idrhs(3)=idarray(3)
      call  boundaryScheme_edges_PhyBC_Inlet_surf_rhs_i1_j1_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_edges_PhyBC_Low_surf_rhs_j1_i1_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (2)
      call  boundaryScheme_edges_i1_jmax_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_edges_i1_jmax_1(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_edges_i1_jmax_2(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_edges_i1_jmax_3(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      idrhs(1)=idarray(1)
      idrhs(4)=idarray(4)
      call  boundaryScheme_edges_PhyBC_Inlet_surf_rhs_i1_jmax_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_edges_PhyBC_Top_surf_rhs_jmax_i1_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (3)
      call  boundaryScheme_edges_imax_j1_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_edges_imax_j1_1(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_edges_imax_j1_2(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_edges_imax_j1_3(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      idrhs(2)=idarray(2)
      idrhs(3)=idarray(3)
      call  boundaryScheme_edges_PhyBC_Out_surf_rhs_imax_j1_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_edges_PhyBC_Low_surf_rhs_j1_imax_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (4)
      call  boundaryScheme_edges_imax_jmax_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_edges_imax_jmax_1(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_edges_imax_jmax_2(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_edges_imax_jmax_3(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      idrhs(2)=idarray(2)
      idrhs(4)=idarray(4)
      call  boundaryScheme_edges_PhyBC_Out_surf_rhs_imax_jmax_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_edges_PhyBC_Top_surf_rhs_jmax_imax_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (5)
      call  boundaryScheme_faces_i1_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_faces_i1_1(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_faces_i1_2(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_faces_i1_3(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      idrhs(1)=idarray(1)
      call  phyInlet_surfbc_faces_i1_8(param_float,hlo,neq,neqst,ind,idloop,idarray,nx,ny,nz,sizeblck,q,qst,q1,q2,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (6)
      call  boundaryScheme_faces_imax_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_faces_imax_1(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_faces_imax_2(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_faces_imax_3(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      idrhs(2)=idarray(2)
      call  phyOut_surfbc_faces_imax_9(param_float,hlo,neq,neqst,ind,idloop,idarray,nx,ny,nz,sizeblck,q,qst,q1,q2,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (7)
      call  boundaryScheme_faces_j1_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_faces_j1_1(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_faces_j1_2(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_faces_j1_3(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      idrhs(3)=idarray(3)
      call  phyLow_surfbc_faces_j1_6(param_float,hlo,neq,neqst,ind,idloop,idarray,nx,ny,nz,sizeblck,q,qst,q1,q2,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (8)
      call  boundaryScheme_faces_jmax_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_faces_jmax_1(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_faces_jmax_2(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_faces_jmax_3(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      idrhs(4)=idarray(4)
      call  phyTop_surfbc_faces_jmax_10(param_float,hlo,neq,neqst,ind,idloop,idarray,nx,ny,nz,sizeblck,q,qst,q1,q2,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

	END SELECT

enddo
 
!$OMP DO SCHEDULE(GUIDED,4) COLLAPSE(2)       
  do bk=idrhs(5),idrhs(6),size_bk  
    do bj=idrhs(3),idrhs(4),size_bj 
      do bi=idrhs(1),idrhs(2),size_bi 
 
idloop(6) = min( bk+size_bk, idrhs(6)+1)-1
idloop(4) = min( bj+size_bj, idrhs(4)+1)-1
idloop(2) = min( bi+size_bi, idrhs(2)+1)-1
    
idloop(5) = bk
idloop(3) = bj
idloop(1) = bi 


if (nrk == 1) then

     do j=idloop(3),idloop(4) 
 
      do i=idloop(1),idloop(2) 

q1(i,j,1) = q(i,j,1)
q1(i,j,2) = q(i,j,2)*q(i,j,1)
q1(i,j,3) = q(i,j,3)*q(i,j,1)
q1(i,j,4) = q(i,j,4)*q(i,j,1)
q1(i,j,5) = q(i,j,5)*q(i,j,1)

     enddo
   enddo

endif  

     do j=idloop(3),idloop(4) 
 
      do i=idloop(1),idloop(2) 

q2(i,j,1) = q1(i,j,1)
q1(i,j,1) = param_float(4)*rk2(nrk)*rhs(i,j,1) + q1(i,j,1)
q2(i,j,2) = q1(i,j,2)
q1(i,j,2) = param_float(4)*rk2(nrk)*rhs(i,j,2) + q1(i,j,2)
q2(i,j,3) = q1(i,j,3)
q1(i,j,3) = param_float(4)*rk2(nrk)*rhs(i,j,3) + q1(i,j,3)
q2(i,j,4) = q1(i,j,4)
q1(i,j,4) = param_float(4)*rk2(nrk)*rhs(i,j,4) + q1(i,j,4)
q2(i,j,5) = q1(i,j,5)
q1(i,j,5) = param_float(4)*rk2(nrk)*rhs(i,j,5) + q1(i,j,5)
     
     enddo
   enddo

     do j=idloop(3),idloop(4) 
 
      do i=idloop(1),idloop(2) 

q (i,j,1) = param_float(4)*rk1(nrk)*rhs(i,j,1) + q2(i,j,1)

q (i,j,2) = param_float(4)*rk1(nrk)*rhs(i,j,2) + q2(i,j,2)

q (i,j,3) = param_float(4)*rk1(nrk)*rhs(i,j,3) + q2(i,j,3)

q (i,j,4) = param_float(4)*rk1(nrk)*rhs(i,j,4) + q2(i,j,4)

q (i,j,5) = param_float(4)*rk1(nrk)*rhs(i,j,5) + q2(i,j,5)

one_over_rho = 1.0_wp/q(i,j,1)
q(i,j,2) = q(i,j,2)*one_over_rho
q(i,j,3) = q(i,j,3)*one_over_rho
q(i,j,4) = q(i,j,4)*one_over_rho
q(i,j,5) = q(i,j,5)*one_over_rho
     
     enddo
   enddo

    enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO NOWAIT





do ibc=1,nbc

    SELECT CASE (bc(ibc))

CASE (1)
      call  boundaryScheme_edges_PhyBC_Inlet_surf_q_i1_j1_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_edges_PhyBC_Low_surf_q_j1_i1_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (2)
      call  boundaryScheme_edges_PhyBC_Inlet_surf_q_i1_jmax_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (3)
      call  boundaryScheme_edges_PhyBC_Low_surf_q_j1_imax_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (5)
      call  phyInlet_surfbc_faces_i1_7(param_float,hlo,neq,neqst,ind,idloop,idarray,nx,ny,nz,sizeblck,q,qst,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (7)
      call  phyLow_surfbc_faces_j1_5(param_float,hlo,neq,neqst,ind,idloop,idarray,nx,ny,nz,sizeblck,q,qst,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

	END SELECT

enddo



!$OMP END PARALLEL



 end subroutine rk3_stepk

 subroutine applybc(param_int,param_float,data_float)

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
  real(wp), intent(inout) ::  data_float(*)
  integer,  intent(in)    ::   param_int(*)
 
 !f2py integer optional, intent(in)    :: type_st = 0

  integer :: shiftstore
  integer :: shiftfacei,shiftfacej,shiftfacek
  integer :: shiftedgeij,shiftedgejk,shiftedgeik
  integer :: nbc
  integer :: nfacei,nfacej,nfacek
  integer :: nedgeij,nedgejk
  integer :: nvar_f(3),nvar_e(3)

  nbc       = param_int(13 + param_int(5) + param_int(6)  )
  
  nvar_f(1) = param_int(13 + param_int(5) + param_int(6) + 1 + nbc      )
  nvar_f(2) = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 1  )
  nvar_f(3) = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 2  ) 

  nvar_e(1) = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 3      )
  nvar_e(2) = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 3 + 1  )
  nvar_e(3) = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 3 + 2  ) 
   
  nfacei    = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 6  )
  nfacej    = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 7  )
  nfacek    = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 8  )
  nedgeij   = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 9  )
  nedgejk   = param_int(13 + param_int(5) + param_int(6) + 1 + nbc + 10 )

  if ( param_int(6) .gt. 0 ) then
   shiftstore = param_int(6)*param_int(9)
  else
   shiftstore = 1
  endif

  shiftfacei  = shiftstore  + nfacei
  shiftfacej  = shiftfacei  + nfacej
  shiftfacek  = shiftfacej  + nfacek
  shiftedgeij = shiftfacek  + nedgeij
  shiftedgejk = shiftedgeij + nedgejk

  call applybcker(param_float                                                          ,&
                  param_int(1)                                                   ,&
                  param_int(5)                                                  ,&
                  param_int(6)                                                ,&
                  param_int(13)                                                  ,&
                  param_int(2),param_int(3) , param_int(4)              ,&
                  param_int(10)                                             ,&
                  param_int(13 + param_int(5) + param_int(6))    ,&
                  param_int(13 + param_int(5) + param_int(6) + 1),&                  
                  nvar_f(1),nvar_f(2),nvar_f(3),&
                  nvar_e(1),nvar_e(2),nvar_e(3),&
                  data_float(1)                                                        ,& ! q
                  data_float(1+param_int(7)*3)                               ,& ! rhs
                  data_float(1+param_int(7)*4)                               ,& ! stored (if any)
                  data_float(1+param_int(7)*4 + shiftstore)                  ,& ! qbcface_i  (if any)
                  data_float(1+param_int(7)*4 + shiftfacei)                  ,& ! qbcface_j  (if any)
                  data_float(1+param_int(7)*4 + shiftfacej)                  ,& ! qbcface_k  (if any)
                  data_float(1+param_int(7)*4 + shiftfacek)                  ,& ! qbcedge_ij (if any)
                  data_float(1+param_int(7)*4 + shiftedgeij)                 ,& ! qbcedge_jk (if any)
                  data_float(1+param_int(7)*4 + shiftedgejk)                 )  ! qbcedge_ik (if any)         


  end subroutine applybc

 subroutine applybcker(param_float,hlo,neq,neqst,ind,nx,ny,nz,sizeblck,&
                     nbc,bc,&
                     nfacei,nfacej,nfacek,&
                     nedgeij,nedgejk,nedgeik,&
                     q,rhs,qst,& 
                     qface_i , qface_j ,qface_k,&
                     qedge_ij, qedge_jk,qedge_ik)

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


  
  integer :: i,j,k
  integer :: bi,bj,bk
  integer :: biend,bjend,bkend
  integer :: size_bi,size_bj,size_bk
  integer :: ibc

  real(wp), intent(in)    :: param_float(*)

  integer, intent(in) :: nx,ny,nz,hlo,neq,neqst
  integer, intent(in) :: nfacei,nfacej,nfacek,nedgeij,nedgejk,nedgeik
  integer, intent(in) :: ind(1:neq+neqst)
  integer, intent(in) :: sizeblck(3)
  integer, intent(in) :: bc(nbc),nbc

real(wp),intent(inout) :: q(1-hlo:nx+hlo,1-hlo:ny+hlo,neq),&
                        rhs(1-hlo:nx+hlo,1-hlo:ny+hlo,neq),&
                        qst(1-hlo:nx+hlo,1-hlo:ny+hlo,neqst)

real(wp),intent(inout) :: qface_i(1-hlo:ny+hlo,nfacei),&
                       qface_j(1-hlo:nx+hlo,nfacej),&
                       qface_k(1),&
                       qedge_ij(1),&
                       qedge_jk(1),&
                       qedge_ik(1)
  
  integer :: idloop(6),idarray(6),indvars(neq),indvarsst(neqst)
  integer :: nvar_f(3),nvar_e(3)

  indvars   = ind(1:neq)
  indvarsst = ind(1+neq:neq+neqst) 

!f2py intent(in)    :: nx,ny,nz,nrk,q,nvar_f,nvar_e,neq,nbc
!f2py intent(inout) :: qst

nvar_f(1) = nfacei
nvar_f(2) = nfacej
nvar_f(3) = nfacek

nvar_e(1) = nedgeij
nvar_e(2) = nedgejk
nvar_e(3) = nedgeik

      
  idarray(1) = 1 -hlo
  idarray(2) = nx+hlo

  idarray(3) = 1 -hlo
  idarray(4) = ny+hlo

  idarray(5) = 1 -hlo
  idarray(6) = nz+hlo

size_bk = sizeblck(3)
size_bj = sizeblck(2)
size_bi = sizeblck(1)



do ibc=1,nbc

    SELECT CASE (bc(ibc))

CASE (1)
      call  boundaryScheme_edges_PhyBC_Inlet_surf_q_i1_j1_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
      call  boundaryScheme_edges_PhyBC_Low_surf_q_j1_i1_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (2)
      call  boundaryScheme_edges_PhyBC_Inlet_surf_q_i1_jmax_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (3)
      call  boundaryScheme_edges_PhyBC_Low_surf_q_j1_imax_0(param_float,hlo,ind,idarray,neq,neqst,nx,ny,nz,sizeblck,q,qst,rhs,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (5)
      call  phyInlet_surfbc_faces_i1_7(param_float,hlo,neq,neqst,ind,idloop,idarray,nx,ny,nz,sizeblck,q,qst,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)
CASE (7)
      call  phyLow_surfbc_faces_j1_5(param_float,hlo,neq,neqst,ind,idloop,idarray,nx,ny,nz,sizeblck,q,qst,nvar_f,nvar_e,qface_i,qface_j,qface_k,qedge_ij,qedge_jk,qedge_ik)

	END SELECT

enddo


end subroutine applybcker



subroutine filter(dir,param_int,param_float,data_float)

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
real(wp), intent(inout) ::  data_float(*)
integer,  intent(in)    :: param_int(*)

integer, intent(in)   :: dir

!f2py  intent(in)   :: param_int,param_float,dir
!f2py intent(inout) :: data_float    


 if     (dir == 1) then

    call flt_x(param_float                                                          ,&
               param_int(1)                                                   ,&
               param_int(5)                                                  ,&
               param_int(2),param_int(3) , param_int(4)              ,&
               param_int(10)                                             ,& 
               param_int(13 + param_int(5) + param_int(6) )   ,&
               param_int(13 + param_int(5) + param_int(6) + 1),&                             
               data_float(1)                                                        ,&
               data_float(1+param_int(7)*2)                                ) ! q2 is used for filter working array
    
 elseif (dir == 2)then

    call flt_y(param_float                                                          ,&
               param_int(1)                                                   ,&
               param_int(5)                                                  ,&
               param_int(2),param_int(3) , param_int(4)              ,&
               param_int(10)                                             ,& 
               param_int(13 + param_int(5) + param_int(6) )   ,&
               param_int(13 + param_int(5) + param_int(6) + 1),&                                              
               data_float(1)                                                        ,&
               data_float(1+param_int(7)*2)                                ) ! q2 is used for filter working array

 elseif (dir == 3)then


    call flt_z(param_float                                                          ,&
               param_int(1)                                                   ,&
               param_int(5)                                                  ,&
               param_int(2),param_int(3) , param_int(4)              ,&
               param_int(10)                                             ,&
               param_int(13 + param_int(5) + param_int(6) )   ,&
               param_int(13 + param_int(5) + param_int(6) + 1),&                                               
               data_float(1)                                                        ,&
               data_float(1+param_int(7)*2)                                ) ! q2 is used for filter working array


 endif

 end subroutine filter


subroutine flt_x(param_float,hlo,neq,nx,ny,nz,sizeblck,nbc,bc,q,q2)

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
! LOCAL VARIABLES

  integer, intent(in) :: nx,ny,nz,hlo,neq  
  integer, intent(in) :: sizeblck(3)

  integer :: i,j,k  
  integer :: bi,bj,bk
  integer :: biend,bjend,bkend
  integer :: size_bi,size_bj,size_bk
  integer :: ibc
  integer, intent(in) :: bc(nbc),nbc

real(wp),intent(inout) ::  q(1-hlo:nx+hlo,1-hlo:ny+hlo,neq),&
                          q2(1-hlo:nx+hlo,1-hlo:ny+hlo,neq)
  
  integer :: idloop(6),idarray(6),idflt(6)


!f2py intent(in)    :: nx,ny,nz,hlo,neq,indvars,param_float
!f2py intent(inout) :: q

!$OMP PARALLEL DEFAULT(SHARED) private(idloop)

idloop(1) = 1 
idloop(2) = nx
idloop(3) = 1 
idloop(4) = ny
idloop(5) = 1 
idloop(6) = nz
    
idflt(1) = 1 
idflt(2) = nx

idflt(3) = 1 
idflt(4) = ny

idflt(5) = 1 
idflt(6) = nz

idarray(1) = 1 -hlo
idarray(2) = nx+hlo

idarray(3) = 1 -hlo
idarray(4) = ny+hlo

idarray(5) = 1 -hlo
idarray(6) = nz+hlo

size_bk = sizeblck(3)
size_bj = sizeblck(2)
size_bi = sizeblck(1)




do ibc=1,nbc

    SELECT CASE (bc(ibc))


 CASE (1)

   idloop(3) = idarray(3)

 CASE (2)

   idloop(4) = idarray(4)

 CASE (3)

   idloop(3) = idarray(3)

 CASE (4)

   idloop(4) = idarray(4)

 CASE (5)

        do j=idloop(3),idloop(4) 
   
   q2(1-4+0,j,1) =0.1604411764705_wp*q(1-4+0,j,1) &
               -0.2325_wp*q(1-4+1,j,1) &
               +0.0895588235295_wp*q(1-4+2,j,1) &
               -0.0175_wp*q(1-4+3,j,1) 
   
   q2(1-4+0,j,2) =0.1604411764705_wp*q(1-4+0,j,2) &
               -0.2325_wp*q(1-4+1,j,2) &
               +0.0895588235295_wp*q(1-4+2,j,2) &
               -0.0175_wp*q(1-4+3,j,2) 
   
   q2(1-4+0,j,3) =0.1604411764705_wp*q(1-4+0,j,3) &
               -0.2325_wp*q(1-4+1,j,3) &
               +0.0895588235295_wp*q(1-4+2,j,3) &
               -0.0175_wp*q(1-4+3,j,3) 
   
   q2(1-4+0,j,4) =0.1604411764705_wp*q(1-4+0,j,4) &
               -0.2325_wp*q(1-4+1,j,4) &
               +0.0895588235295_wp*q(1-4+2,j,4) &
               -0.0175_wp*q(1-4+3,j,4) 
   
   q2(1-4+0,j,5) =0.1604411764705_wp*q(1-4+0,j,5) &
               -0.2325_wp*q(1-4+1,j,5) &
               +0.0895588235295_wp*q(1-4+2,j,5) &
               -0.0175_wp*q(1-4+3,j,5) 
   
   q2(1-4+1,j,1) =-0.042888704485_wp*q(1-4+0,j,1) &
               +0.138814085762_wp*q(1-4+1,j,1) &
               -0.1784240360865_wp*q(1-4+2,j,1) &
               +0.111559546536_wp*q(1-4+3,j,1) &
               -0.0286735324325_wp*q(1-4+4,j,1) &
               -0.000373632298_wp*q(1-4+5,j,1) &
               -1.37269965e-05_wp*q(1-4+6,j,1) 
   
   q2(1-4+1,j,2) =-0.042888704485_wp*q(1-4+0,j,2) &
               +0.138814085762_wp*q(1-4+1,j,2) &
               -0.1784240360865_wp*q(1-4+2,j,2) &
               +0.111559546536_wp*q(1-4+3,j,2) &
               -0.0286735324325_wp*q(1-4+4,j,2) &
               -0.000373632298_wp*q(1-4+5,j,2) &
               -1.37269965e-05_wp*q(1-4+6,j,2) 
   
   q2(1-4+1,j,3) =-0.042888704485_wp*q(1-4+0,j,3) &
               +0.138814085762_wp*q(1-4+1,j,3) &
               -0.1784240360865_wp*q(1-4+2,j,3) &
               +0.111559546536_wp*q(1-4+3,j,3) &
               -0.0286735324325_wp*q(1-4+4,j,3) &
               -0.000373632298_wp*q(1-4+5,j,3) &
               -1.37269965e-05_wp*q(1-4+6,j,3) 
   
   q2(1-4+1,j,4) =-0.042888704485_wp*q(1-4+0,j,4) &
               +0.138814085762_wp*q(1-4+1,j,4) &
               -0.1784240360865_wp*q(1-4+2,j,4) &
               +0.111559546536_wp*q(1-4+3,j,4) &
               -0.0286735324325_wp*q(1-4+4,j,4) &
               -0.000373632298_wp*q(1-4+5,j,4) &
               -1.37269965e-05_wp*q(1-4+6,j,4) 
   
   q2(1-4+1,j,5) =-0.042888704485_wp*q(1-4+0,j,5) &
               +0.138814085762_wp*q(1-4+1,j,5) &
               -0.1784240360865_wp*q(1-4+2,j,5) &
               +0.111559546536_wp*q(1-4+3,j,5) &
               -0.0286735324325_wp*q(1-4+4,j,5) &
               -0.000373632298_wp*q(1-4+5,j,5) &
               -1.37269965e-05_wp*q(1-4+6,j,5) 
   
   q2(1-4+2,j,1) =0.052523901012_wp*q(1-4+0,j,1) &
               -0.206299133811_wp*q(1-4+1,j,1) &
               +0.35352799825_wp*q(1-4+2,j,1) &
               -0.348142394842_wp*q(1-4+3,j,1) &
               +0.181481803619_wp*q(1-4+4,j,1) &
               +0.00944080437_wp*q(1-4+5,j,1) &
               -0.077675100452_wp*q(1-4+6,j,1) &
               +0.044887364863_wp*q(1-4+7,j,1) &
               -0.009971961849_wp*q(1-4+8,j,1) &
               +0.00011335942_wp*q(1-4+9,j,1) &
               +0.00011335942_wp*q(1-4+10,j,1) 
   
   q2(1-4+2,j,2) =0.052523901012_wp*q(1-4+0,j,2) &
               -0.206299133811_wp*q(1-4+1,j,2) &
               +0.35352799825_wp*q(1-4+2,j,2) &
               -0.348142394842_wp*q(1-4+3,j,2) &
               +0.181481803619_wp*q(1-4+4,j,2) &
               +0.00944080437_wp*q(1-4+5,j,2) &
               -0.077675100452_wp*q(1-4+6,j,2) &
               +0.044887364863_wp*q(1-4+7,j,2) &
               -0.009971961849_wp*q(1-4+8,j,2) &
               +0.00011335942_wp*q(1-4+9,j,2) &
               +0.00011335942_wp*q(1-4+10,j,2) 
   
   q2(1-4+2,j,3) =0.052523901012_wp*q(1-4+0,j,3) &
               -0.206299133811_wp*q(1-4+1,j,3) &
               +0.35352799825_wp*q(1-4+2,j,3) &
               -0.348142394842_wp*q(1-4+3,j,3) &
               +0.181481803619_wp*q(1-4+4,j,3) &
               +0.00944080437_wp*q(1-4+5,j,3) &
               -0.077675100452_wp*q(1-4+6,j,3) &
               +0.044887364863_wp*q(1-4+7,j,3) &
               -0.009971961849_wp*q(1-4+8,j,3) &
               +0.00011335942_wp*q(1-4+9,j,3) &
               +0.00011335942_wp*q(1-4+10,j,3) 
   
   q2(1-4+2,j,4) =0.052523901012_wp*q(1-4+0,j,4) &
               -0.206299133811_wp*q(1-4+1,j,4) &
               +0.35352799825_wp*q(1-4+2,j,4) &
               -0.348142394842_wp*q(1-4+3,j,4) &
               +0.181481803619_wp*q(1-4+4,j,4) &
               +0.00944080437_wp*q(1-4+5,j,4) &
               -0.077675100452_wp*q(1-4+6,j,4) &
               +0.044887364863_wp*q(1-4+7,j,4) &
               -0.009971961849_wp*q(1-4+8,j,4) &
               +0.00011335942_wp*q(1-4+9,j,4) &
               +0.00011335942_wp*q(1-4+10,j,4) 
   
   q2(1-4+2,j,5) =0.052523901012_wp*q(1-4+0,j,5) &
               -0.206299133811_wp*q(1-4+1,j,5) &
               +0.35352799825_wp*q(1-4+2,j,5) &
               -0.348142394842_wp*q(1-4+3,j,5) &
               +0.181481803619_wp*q(1-4+4,j,5) &
               +0.00944080437_wp*q(1-4+5,j,5) &
               -0.077675100452_wp*q(1-4+6,j,5) &
               +0.044887364863_wp*q(1-4+7,j,5) &
               -0.009971961849_wp*q(1-4+8,j,5) &
               +0.00011335942_wp*q(1-4+9,j,5) &
               +0.00011335942_wp*q(1-4+10,j,5) 
   
   q2(1-4+3,j,1) =-5.459601e-05_wp*q(1-4+0,j,1) &
               +0.042124772446_wp*q(1-4+1,j,1) &
               -0.173103107841_wp*q(1-4+2,j,1) &
               +0.299615871352_wp*q(1-4+3,j,1) &
               -0.276543612935_wp*q(1-4+4,j,1) &
               +0.131223506571_wp*q(1-4+5,j,1) &
               -0.023424966418_wp*q(1-4+6,j,1) &
               +0.013937561779_wp*q(1-4+7,j,1) &
               -0.024565095706_wp*q(1-4+8,j,1) &
               +0.013098287852_wp*q(1-4+9,j,1) &
               -0.00230862109_wp*q(1-4+10,j,1) 
   
   q2(1-4+3,j,2) =-5.459601e-05_wp*q(1-4+0,j,2) &
               +0.042124772446_wp*q(1-4+1,j,2) &
               -0.173103107841_wp*q(1-4+2,j,2) &
               +0.299615871352_wp*q(1-4+3,j,2) &
               -0.276543612935_wp*q(1-4+4,j,2) &
               +0.131223506571_wp*q(1-4+5,j,2) &
               -0.023424966418_wp*q(1-4+6,j,2) &
               +0.013937561779_wp*q(1-4+7,j,2) &
               -0.024565095706_wp*q(1-4+8,j,2) &
               +0.013098287852_wp*q(1-4+9,j,2) &
               -0.00230862109_wp*q(1-4+10,j,2) 
   
   q2(1-4+3,j,3) =-5.459601e-05_wp*q(1-4+0,j,3) &
               +0.042124772446_wp*q(1-4+1,j,3) &
               -0.173103107841_wp*q(1-4+2,j,3) &
               +0.299615871352_wp*q(1-4+3,j,3) &
               -0.276543612935_wp*q(1-4+4,j,3) &
               +0.131223506571_wp*q(1-4+5,j,3) &
               -0.023424966418_wp*q(1-4+6,j,3) &
               +0.013937561779_wp*q(1-4+7,j,3) &
               -0.024565095706_wp*q(1-4+8,j,3) &
               +0.013098287852_wp*q(1-4+9,j,3) &
               -0.00230862109_wp*q(1-4+10,j,3) 
   
   q2(1-4+3,j,4) =-5.459601e-05_wp*q(1-4+0,j,4) &
               +0.042124772446_wp*q(1-4+1,j,4) &
               -0.173103107841_wp*q(1-4+2,j,4) &
               +0.299615871352_wp*q(1-4+3,j,4) &
               -0.276543612935_wp*q(1-4+4,j,4) &
               +0.131223506571_wp*q(1-4+5,j,4) &
               -0.023424966418_wp*q(1-4+6,j,4) &
               +0.013937561779_wp*q(1-4+7,j,4) &
               -0.024565095706_wp*q(1-4+8,j,4) &
               +0.013098287852_wp*q(1-4+9,j,4) &
               -0.00230862109_wp*q(1-4+10,j,4) 
   
   q2(1-4+3,j,5) =-5.459601e-05_wp*q(1-4+0,j,5) &
               +0.042124772446_wp*q(1-4+1,j,5) &
               -0.173103107841_wp*q(1-4+2,j,5) &
               +0.299615871352_wp*q(1-4+3,j,5) &
               -0.276543612935_wp*q(1-4+4,j,5) &
               +0.131223506571_wp*q(1-4+5,j,5) &
               -0.023424966418_wp*q(1-4+6,j,5) &
               +0.013937561779_wp*q(1-4+7,j,5) &
               -0.024565095706_wp*q(1-4+8,j,5) &
               +0.013098287852_wp*q(1-4+9,j,5) &
               -0.00230862109_wp*q(1-4+10,j,5) 
   
        enddo

 CASE (6)

        do j=idloop(3),idloop(4) 
   
   q2(nx+4+0,j,1) =0.1604411764705_wp*q(nx+4+0,j,1) &
               -0.2325_wp*q(nx+4-1,j,1) &
               +0.0895588235295_wp*q(nx+4-2,j,1) &
               -0.0175_wp*q(nx+4-3,j,1) 
   
   q2(nx+4+0,j,2) =0.1604411764705_wp*q(nx+4+0,j,2) &
               -0.2325_wp*q(nx+4-1,j,2) &
               +0.0895588235295_wp*q(nx+4-2,j,2) &
               -0.0175_wp*q(nx+4-3,j,2) 
   
   q2(nx+4+0,j,3) =0.1604411764705_wp*q(nx+4+0,j,3) &
               -0.2325_wp*q(nx+4-1,j,3) &
               +0.0895588235295_wp*q(nx+4-2,j,3) &
               -0.0175_wp*q(nx+4-3,j,3) 
   
   q2(nx+4+0,j,4) =0.1604411764705_wp*q(nx+4+0,j,4) &
               -0.2325_wp*q(nx+4-1,j,4) &
               +0.0895588235295_wp*q(nx+4-2,j,4) &
               -0.0175_wp*q(nx+4-3,j,4) 
   
   q2(nx+4+0,j,5) =0.1604411764705_wp*q(nx+4+0,j,5) &
               -0.2325_wp*q(nx+4-1,j,5) &
               +0.0895588235295_wp*q(nx+4-2,j,5) &
               -0.0175_wp*q(nx+4-3,j,5) 
   
   q2(nx+4-1,j,1) =-0.042888704485_wp*q(nx+4+0,j,1) &
               +0.138814085762_wp*q(nx+4-1,j,1) &
               -0.1784240360865_wp*q(nx+4-2,j,1) &
               +0.111559546536_wp*q(nx+4-3,j,1) &
               -0.0286735324325_wp*q(nx+4-4,j,1) &
               -0.000373632298_wp*q(nx+4-5,j,1) &
               -1.37269965e-05_wp*q(nx+4-6,j,1) 
   
   q2(nx+4-1,j,2) =-0.042888704485_wp*q(nx+4+0,j,2) &
               +0.138814085762_wp*q(nx+4-1,j,2) &
               -0.1784240360865_wp*q(nx+4-2,j,2) &
               +0.111559546536_wp*q(nx+4-3,j,2) &
               -0.0286735324325_wp*q(nx+4-4,j,2) &
               -0.000373632298_wp*q(nx+4-5,j,2) &
               -1.37269965e-05_wp*q(nx+4-6,j,2) 
   
   q2(nx+4-1,j,3) =-0.042888704485_wp*q(nx+4+0,j,3) &
               +0.138814085762_wp*q(nx+4-1,j,3) &
               -0.1784240360865_wp*q(nx+4-2,j,3) &
               +0.111559546536_wp*q(nx+4-3,j,3) &
               -0.0286735324325_wp*q(nx+4-4,j,3) &
               -0.000373632298_wp*q(nx+4-5,j,3) &
               -1.37269965e-05_wp*q(nx+4-6,j,3) 
   
   q2(nx+4-1,j,4) =-0.042888704485_wp*q(nx+4+0,j,4) &
               +0.138814085762_wp*q(nx+4-1,j,4) &
               -0.1784240360865_wp*q(nx+4-2,j,4) &
               +0.111559546536_wp*q(nx+4-3,j,4) &
               -0.0286735324325_wp*q(nx+4-4,j,4) &
               -0.000373632298_wp*q(nx+4-5,j,4) &
               -1.37269965e-05_wp*q(nx+4-6,j,4) 
   
   q2(nx+4-1,j,5) =-0.042888704485_wp*q(nx+4+0,j,5) &
               +0.138814085762_wp*q(nx+4-1,j,5) &
               -0.1784240360865_wp*q(nx+4-2,j,5) &
               +0.111559546536_wp*q(nx+4-3,j,5) &
               -0.0286735324325_wp*q(nx+4-4,j,5) &
               -0.000373632298_wp*q(nx+4-5,j,5) &
               -1.37269965e-05_wp*q(nx+4-6,j,5) 
   
   q2(nx+4-2,j,1) =0.052523901012_wp*q(nx+4+0,j,1) &
               -0.206299133811_wp*q(nx+4-1,j,1) &
               +0.35352799825_wp*q(nx+4-2,j,1) &
               -0.348142394842_wp*q(nx+4-3,j,1) &
               +0.181481803619_wp*q(nx+4-4,j,1) &
               +0.00944080437_wp*q(nx+4-5,j,1) &
               -0.077675100452_wp*q(nx+4-6,j,1) &
               +0.044887364863_wp*q(nx+4-7,j,1) &
               -0.009971961849_wp*q(nx+4-8,j,1) &
               +0.00011335942_wp*q(nx+4-9,j,1) &
               +0.00011335942_wp*q(nx+4-10,j,1) 
   
   q2(nx+4-2,j,2) =0.052523901012_wp*q(nx+4+0,j,2) &
               -0.206299133811_wp*q(nx+4-1,j,2) &
               +0.35352799825_wp*q(nx+4-2,j,2) &
               -0.348142394842_wp*q(nx+4-3,j,2) &
               +0.181481803619_wp*q(nx+4-4,j,2) &
               +0.00944080437_wp*q(nx+4-5,j,2) &
               -0.077675100452_wp*q(nx+4-6,j,2) &
               +0.044887364863_wp*q(nx+4-7,j,2) &
               -0.009971961849_wp*q(nx+4-8,j,2) &
               +0.00011335942_wp*q(nx+4-9,j,2) &
               +0.00011335942_wp*q(nx+4-10,j,2) 
   
   q2(nx+4-2,j,3) =0.052523901012_wp*q(nx+4+0,j,3) &
               -0.206299133811_wp*q(nx+4-1,j,3) &
               +0.35352799825_wp*q(nx+4-2,j,3) &
               -0.348142394842_wp*q(nx+4-3,j,3) &
               +0.181481803619_wp*q(nx+4-4,j,3) &
               +0.00944080437_wp*q(nx+4-5,j,3) &
               -0.077675100452_wp*q(nx+4-6,j,3) &
               +0.044887364863_wp*q(nx+4-7,j,3) &
               -0.009971961849_wp*q(nx+4-8,j,3) &
               +0.00011335942_wp*q(nx+4-9,j,3) &
               +0.00011335942_wp*q(nx+4-10,j,3) 
   
   q2(nx+4-2,j,4) =0.052523901012_wp*q(nx+4+0,j,4) &
               -0.206299133811_wp*q(nx+4-1,j,4) &
               +0.35352799825_wp*q(nx+4-2,j,4) &
               -0.348142394842_wp*q(nx+4-3,j,4) &
               +0.181481803619_wp*q(nx+4-4,j,4) &
               +0.00944080437_wp*q(nx+4-5,j,4) &
               -0.077675100452_wp*q(nx+4-6,j,4) &
               +0.044887364863_wp*q(nx+4-7,j,4) &
               -0.009971961849_wp*q(nx+4-8,j,4) &
               +0.00011335942_wp*q(nx+4-9,j,4) &
               +0.00011335942_wp*q(nx+4-10,j,4) 
   
   q2(nx+4-2,j,5) =0.052523901012_wp*q(nx+4+0,j,5) &
               -0.206299133811_wp*q(nx+4-1,j,5) &
               +0.35352799825_wp*q(nx+4-2,j,5) &
               -0.348142394842_wp*q(nx+4-3,j,5) &
               +0.181481803619_wp*q(nx+4-4,j,5) &
               +0.00944080437_wp*q(nx+4-5,j,5) &
               -0.077675100452_wp*q(nx+4-6,j,5) &
               +0.044887364863_wp*q(nx+4-7,j,5) &
               -0.009971961849_wp*q(nx+4-8,j,5) &
               +0.00011335942_wp*q(nx+4-9,j,5) &
               +0.00011335942_wp*q(nx+4-10,j,5) 
   
   q2(nx+4-3,j,1) =-5.459601e-05_wp*q(nx+4+0,j,1) &
               +0.042124772446_wp*q(nx+4-1,j,1) &
               -0.173103107841_wp*q(nx+4-2,j,1) &
               +0.299615871352_wp*q(nx+4-3,j,1) &
               -0.276543612935_wp*q(nx+4-4,j,1) &
               +0.131223506571_wp*q(nx+4-5,j,1) &
               -0.023424966418_wp*q(nx+4-6,j,1) &
               +0.013937561779_wp*q(nx+4-7,j,1) &
               -0.024565095706_wp*q(nx+4-8,j,1) &
               +0.013098287852_wp*q(nx+4-9,j,1) &
               -0.00230862109_wp*q(nx+4-10,j,1) 
   
   q2(nx+4-3,j,2) =-5.459601e-05_wp*q(nx+4+0,j,2) &
               +0.042124772446_wp*q(nx+4-1,j,2) &
               -0.173103107841_wp*q(nx+4-2,j,2) &
               +0.299615871352_wp*q(nx+4-3,j,2) &
               -0.276543612935_wp*q(nx+4-4,j,2) &
               +0.131223506571_wp*q(nx+4-5,j,2) &
               -0.023424966418_wp*q(nx+4-6,j,2) &
               +0.013937561779_wp*q(nx+4-7,j,2) &
               -0.024565095706_wp*q(nx+4-8,j,2) &
               +0.013098287852_wp*q(nx+4-9,j,2) &
               -0.00230862109_wp*q(nx+4-10,j,2) 
   
   q2(nx+4-3,j,3) =-5.459601e-05_wp*q(nx+4+0,j,3) &
               +0.042124772446_wp*q(nx+4-1,j,3) &
               -0.173103107841_wp*q(nx+4-2,j,3) &
               +0.299615871352_wp*q(nx+4-3,j,3) &
               -0.276543612935_wp*q(nx+4-4,j,3) &
               +0.131223506571_wp*q(nx+4-5,j,3) &
               -0.023424966418_wp*q(nx+4-6,j,3) &
               +0.013937561779_wp*q(nx+4-7,j,3) &
               -0.024565095706_wp*q(nx+4-8,j,3) &
               +0.013098287852_wp*q(nx+4-9,j,3) &
               -0.00230862109_wp*q(nx+4-10,j,3) 
   
   q2(nx+4-3,j,4) =-5.459601e-05_wp*q(nx+4+0,j,4) &
               +0.042124772446_wp*q(nx+4-1,j,4) &
               -0.173103107841_wp*q(nx+4-2,j,4) &
               +0.299615871352_wp*q(nx+4-3,j,4) &
               -0.276543612935_wp*q(nx+4-4,j,4) &
               +0.131223506571_wp*q(nx+4-5,j,4) &
               -0.023424966418_wp*q(nx+4-6,j,4) &
               +0.013937561779_wp*q(nx+4-7,j,4) &
               -0.024565095706_wp*q(nx+4-8,j,4) &
               +0.013098287852_wp*q(nx+4-9,j,4) &
               -0.00230862109_wp*q(nx+4-10,j,4) 
   
   q2(nx+4-3,j,5) =-5.459601e-05_wp*q(nx+4+0,j,5) &
               +0.042124772446_wp*q(nx+4-1,j,5) &
               -0.173103107841_wp*q(nx+4-2,j,5) &
               +0.299615871352_wp*q(nx+4-3,j,5) &
               -0.276543612935_wp*q(nx+4-4,j,5) &
               +0.131223506571_wp*q(nx+4-5,j,5) &
               -0.023424966418_wp*q(nx+4-6,j,5) &
               +0.013937561779_wp*q(nx+4-7,j,5) &
               -0.024565095706_wp*q(nx+4-8,j,5) &
               +0.013098287852_wp*q(nx+4-9,j,5) &
               -0.00230862109_wp*q(nx+4-10,j,5) 
   
        enddo

 CASE (7)

   idflt(3) = idarray(3)

 CASE (8)

   idflt(4) = idarray(4)

    END SELECT

enddo

!$OMP DO SCHEDULE(GUIDED,4) COLLAPSE(2) 
  do bk=idflt(5),idflt(6),size_bk  
    do bj=idflt(3),idflt(4),size_bj 
      do bi=idflt(1),idflt(2),size_bi 
    
idloop(6) = min( bk+size_bk, idflt(6)+1)-1
idloop(4) = min( bj+size_bj, idflt(4)+1)-1
idloop(2) = min( bi+size_bi, idflt(2)+1)-1

idloop(5) = bk
idloop(3) = bj
idloop(1) = bi 

     do j=idloop(3),idloop(4) 
 
      do i=idloop(1),idloop(2) 

q2(i,j,1) =0.0625_wp*q(i-2,j,1) &
            -0.25_wp*q(i-1,j,1) &
            +0.375_wp*q(i+0,j,1) &
            -0.25_wp*q(i+1,j,1) &
            +0.0625_wp*q(i+2,j,1) 

q2(i,j,2) =0.0625_wp*q(i-2,j,2) &
            -0.25_wp*q(i-1,j,2) &
            +0.375_wp*q(i+0,j,2) &
            -0.25_wp*q(i+1,j,2) &
            +0.0625_wp*q(i+2,j,2) 

q2(i,j,3) =0.0625_wp*q(i-2,j,3) &
            -0.25_wp*q(i-1,j,3) &
            +0.375_wp*q(i+0,j,3) &
            -0.25_wp*q(i+1,j,3) &
            +0.0625_wp*q(i+2,j,3) 

q2(i,j,4) =0.0625_wp*q(i-2,j,4) &
            -0.25_wp*q(i-1,j,4) &
            +0.375_wp*q(i+0,j,4) &
            -0.25_wp*q(i+1,j,4) &
            +0.0625_wp*q(i+2,j,4) 

q2(i,j,5) =0.0625_wp*q(i-2,j,5) &
            -0.25_wp*q(i-1,j,5) &
            +0.375_wp*q(i+0,j,5) &
            -0.25_wp*q(i+1,j,5) &
            +0.0625_wp*q(i+2,j,5) 

     
     enddo
   enddo

    enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO
 
!$OMP DO SCHEDULE(GUIDED,4) COLLAPSE(2) 
  do bk=idflt(5),idflt(6),size_bk  
    do bj=idflt(3),idflt(4),size_bj 
      do bi=idflt(1),idflt(2),size_bi 
    
idloop(6) = min( bk+size_bk, idflt(6)+1)-1
idloop(4) = min( bj+size_bj, idflt(4)+1)-1
idloop(2) = min( bi+size_bi, idflt(2)+1)-1

idloop(5) = bk
idloop(3) = bj
idloop(1) = bi 

     do j=idloop(3),idloop(4) 
 
      do i=idloop(1),idloop(2) 
q(i,j,1) = q(i,j,1) - param_float(5)*q2(i,j,1)

q(i,j,2) = q(i,j,2) - param_float(5)*q2(i,j,2)

q(i,j,3) = q(i,j,3) - param_float(5)*q2(i,j,3)

q(i,j,4) = q(i,j,4) - param_float(5)*q2(i,j,4)

q(i,j,5) = q(i,j,5) - param_float(5)*q2(i,j,5)

     enddo
   enddo
  
     enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO

idloop(1) = 1 
idloop(2) = nx
idloop(3) = 1 
idloop(4) = ny
idloop(5) = 1 
idloop(6) = nz



do ibc=1,nbc

    SELECT CASE (bc(ibc))


 CASE (1)

   idloop(3) = idarray(3)

 CASE (2)

   idloop(4) = idarray(4)

 CASE (3)

   idloop(3) = idarray(3)

 CASE (4)

   idloop(4) = idarray(4)

 CASE (5)

        do j=idloop(3),idloop(4) 
   
   q(1-4+0,j,1) = q(1-4+0,j,1) - param_float(5)*q2(1-4+0,j,1)
   q(1-4+0,j,2) = q(1-4+0,j,2) - param_float(5)*q2(1-4+0,j,2)
   q(1-4+0,j,3) = q(1-4+0,j,3) - param_float(5)*q2(1-4+0,j,3)
   q(1-4+0,j,4) = q(1-4+0,j,4) - param_float(5)*q2(1-4+0,j,4)
   q(1-4+0,j,5) = q(1-4+0,j,5) - param_float(5)*q2(1-4+0,j,5)
   
   q(1-4+1,j,1) = q(1-4+1,j,1) - param_float(5)*q2(1-4+1,j,1)
   q(1-4+1,j,2) = q(1-4+1,j,2) - param_float(5)*q2(1-4+1,j,2)
   q(1-4+1,j,3) = q(1-4+1,j,3) - param_float(5)*q2(1-4+1,j,3)
   q(1-4+1,j,4) = q(1-4+1,j,4) - param_float(5)*q2(1-4+1,j,4)
   q(1-4+1,j,5) = q(1-4+1,j,5) - param_float(5)*q2(1-4+1,j,5)
   
   q(1-4+2,j,1) = q(1-4+2,j,1) - param_float(5)*q2(1-4+2,j,1)
   q(1-4+2,j,2) = q(1-4+2,j,2) - param_float(5)*q2(1-4+2,j,2)
   q(1-4+2,j,3) = q(1-4+2,j,3) - param_float(5)*q2(1-4+2,j,3)
   q(1-4+2,j,4) = q(1-4+2,j,4) - param_float(5)*q2(1-4+2,j,4)
   q(1-4+2,j,5) = q(1-4+2,j,5) - param_float(5)*q2(1-4+2,j,5)
   
   q(1-4+3,j,1) = q(1-4+3,j,1) - param_float(5)*q2(1-4+3,j,1)
   q(1-4+3,j,2) = q(1-4+3,j,2) - param_float(5)*q2(1-4+3,j,2)
   q(1-4+3,j,3) = q(1-4+3,j,3) - param_float(5)*q2(1-4+3,j,3)
   q(1-4+3,j,4) = q(1-4+3,j,4) - param_float(5)*q2(1-4+3,j,4)
   q(1-4+3,j,5) = q(1-4+3,j,5) - param_float(5)*q2(1-4+3,j,5)
   
        enddo

 CASE (6)

        do j=idloop(3),idloop(4) 
   
   q(nx+4+0,j,1) = q(nx+4+0,j,1) - param_float(5)*q2(nx+4+0,j,1)
   q(nx+4+0,j,2) = q(nx+4+0,j,2) - param_float(5)*q2(nx+4+0,j,2)
   q(nx+4+0,j,3) = q(nx+4+0,j,3) - param_float(5)*q2(nx+4+0,j,3)
   q(nx+4+0,j,4) = q(nx+4+0,j,4) - param_float(5)*q2(nx+4+0,j,4)
   q(nx+4+0,j,5) = q(nx+4+0,j,5) - param_float(5)*q2(nx+4+0,j,5)
   
   q(nx+4-1,j,1) = q(nx+4-1,j,1) - param_float(5)*q2(nx+4-1,j,1)
   q(nx+4-1,j,2) = q(nx+4-1,j,2) - param_float(5)*q2(nx+4-1,j,2)
   q(nx+4-1,j,3) = q(nx+4-1,j,3) - param_float(5)*q2(nx+4-1,j,3)
   q(nx+4-1,j,4) = q(nx+4-1,j,4) - param_float(5)*q2(nx+4-1,j,4)
   q(nx+4-1,j,5) = q(nx+4-1,j,5) - param_float(5)*q2(nx+4-1,j,5)
   
   q(nx+4-2,j,1) = q(nx+4-2,j,1) - param_float(5)*q2(nx+4-2,j,1)
   q(nx+4-2,j,2) = q(nx+4-2,j,2) - param_float(5)*q2(nx+4-2,j,2)
   q(nx+4-2,j,3) = q(nx+4-2,j,3) - param_float(5)*q2(nx+4-2,j,3)
   q(nx+4-2,j,4) = q(nx+4-2,j,4) - param_float(5)*q2(nx+4-2,j,4)
   q(nx+4-2,j,5) = q(nx+4-2,j,5) - param_float(5)*q2(nx+4-2,j,5)
   
   q(nx+4-3,j,1) = q(nx+4-3,j,1) - param_float(5)*q2(nx+4-3,j,1)
   q(nx+4-3,j,2) = q(nx+4-3,j,2) - param_float(5)*q2(nx+4-3,j,2)
   q(nx+4-3,j,3) = q(nx+4-3,j,3) - param_float(5)*q2(nx+4-3,j,3)
   q(nx+4-3,j,4) = q(nx+4-3,j,4) - param_float(5)*q2(nx+4-3,j,4)
   q(nx+4-3,j,5) = q(nx+4-3,j,5) - param_float(5)*q2(nx+4-3,j,5)
   
        enddo

    END SELECT

enddo

!$OMP END PARALLEL

 end subroutine flt_x

subroutine flt_y(param_float,hlo,neq,nx,ny,nz,sizeblck,nbc,bc,q,q2)

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
! LOCAL VARIABLES

  integer, intent(in) :: nx,ny,nz,hlo,neq  
  integer, intent(in) :: sizeblck(3)

  integer :: i,j,k  
  integer :: bi,bj,bk
  integer :: biend,bjend,bkend
  integer :: size_bi,size_bj,size_bk
  integer :: ibc
  integer, intent(in) :: bc(nbc),nbc

real(wp),intent(inout) ::  q(1-hlo:nx+hlo,1-hlo:ny+hlo,neq),&
                          q2(1-hlo:nx+hlo,1-hlo:ny+hlo,neq)
  
  integer :: idloop(6),idarray(6),idflt(6)


!f2py intent(in)    :: nx,ny,nz,hlo,neq,indvars,param_float
!f2py intent(inout) :: q

!$OMP PARALLEL DEFAULT(SHARED) private(idloop)

idloop(1) = 1 
idloop(2) = nx
idloop(3) = 1 
idloop(4) = ny
idloop(5) = 1 
idloop(6) = nz
  
idflt(1) = 1 
idflt(2) = nx

idflt(3) = 1 
idflt(4) = ny

idflt(5) = 1 
idflt(6) = nz


idarray(1) = 1 -hlo
idarray(2) = nx+hlo

idarray(3) = 1 -hlo
idarray(4) = ny+hlo

idarray(5) = 1 -hlo
idarray(6) = nz+hlo


size_bk = sizeblck(3)
size_bj = sizeblck(2)
size_bi = sizeblck(1)



do ibc=1,nbc

    SELECT CASE (bc(ibc))


 CASE (1)

   idloop(1) = idarray(1)

 CASE (2)

   idloop(1) = idarray(1)

 CASE (3)

   idloop(2) = idarray(2)

 CASE (4)

   idloop(2) = idarray(2)

 CASE (5)

   idflt(1) = idarray(1)

 CASE (7)

    
         do i=idloop(1),idloop(2) 
   
   q2(i,1-4+0,1) =0.1604411764705_wp*q(i,1-4+0,1) &
               -0.2325_wp*q(i,1-4+1,1) &
               +0.0895588235295_wp*q(i,1-4+2,1) &
               -0.0175_wp*q(i,1-4+3,1) 
   
   q2(i,1-4+0,2) =0.1604411764705_wp*q(i,1-4+0,2) &
               -0.2325_wp*q(i,1-4+1,2) &
               +0.0895588235295_wp*q(i,1-4+2,2) &
               -0.0175_wp*q(i,1-4+3,2) 
   
   q2(i,1-4+0,3) =0.1604411764705_wp*q(i,1-4+0,3) &
               -0.2325_wp*q(i,1-4+1,3) &
               +0.0895588235295_wp*q(i,1-4+2,3) &
               -0.0175_wp*q(i,1-4+3,3) 
   
   q2(i,1-4+0,4) =0.1604411764705_wp*q(i,1-4+0,4) &
               -0.2325_wp*q(i,1-4+1,4) &
               +0.0895588235295_wp*q(i,1-4+2,4) &
               -0.0175_wp*q(i,1-4+3,4) 
   
   q2(i,1-4+0,5) =0.1604411764705_wp*q(i,1-4+0,5) &
               -0.2325_wp*q(i,1-4+1,5) &
               +0.0895588235295_wp*q(i,1-4+2,5) &
               -0.0175_wp*q(i,1-4+3,5) 
   
   q2(i,1-4+1,1) =-0.042888704485_wp*q(i,1-4+0,1) &
               +0.138814085762_wp*q(i,1-4+1,1) &
               -0.1784240360865_wp*q(i,1-4+2,1) &
               +0.111559546536_wp*q(i,1-4+3,1) &
               -0.0286735324325_wp*q(i,1-4+4,1) &
               -0.000373632298_wp*q(i,1-4+5,1) &
               -1.37269965e-05_wp*q(i,1-4+6,1) 
   
   q2(i,1-4+1,2) =-0.042888704485_wp*q(i,1-4+0,2) &
               +0.138814085762_wp*q(i,1-4+1,2) &
               -0.1784240360865_wp*q(i,1-4+2,2) &
               +0.111559546536_wp*q(i,1-4+3,2) &
               -0.0286735324325_wp*q(i,1-4+4,2) &
               -0.000373632298_wp*q(i,1-4+5,2) &
               -1.37269965e-05_wp*q(i,1-4+6,2) 
   
   q2(i,1-4+1,3) =-0.042888704485_wp*q(i,1-4+0,3) &
               +0.138814085762_wp*q(i,1-4+1,3) &
               -0.1784240360865_wp*q(i,1-4+2,3) &
               +0.111559546536_wp*q(i,1-4+3,3) &
               -0.0286735324325_wp*q(i,1-4+4,3) &
               -0.000373632298_wp*q(i,1-4+5,3) &
               -1.37269965e-05_wp*q(i,1-4+6,3) 
   
   q2(i,1-4+1,4) =-0.042888704485_wp*q(i,1-4+0,4) &
               +0.138814085762_wp*q(i,1-4+1,4) &
               -0.1784240360865_wp*q(i,1-4+2,4) &
               +0.111559546536_wp*q(i,1-4+3,4) &
               -0.0286735324325_wp*q(i,1-4+4,4) &
               -0.000373632298_wp*q(i,1-4+5,4) &
               -1.37269965e-05_wp*q(i,1-4+6,4) 
   
   q2(i,1-4+1,5) =-0.042888704485_wp*q(i,1-4+0,5) &
               +0.138814085762_wp*q(i,1-4+1,5) &
               -0.1784240360865_wp*q(i,1-4+2,5) &
               +0.111559546536_wp*q(i,1-4+3,5) &
               -0.0286735324325_wp*q(i,1-4+4,5) &
               -0.000373632298_wp*q(i,1-4+5,5) &
               -1.37269965e-05_wp*q(i,1-4+6,5) 
   
   q2(i,1-4+2,1) =0.052523901012_wp*q(i,1-4+0,1) &
               -0.206299133811_wp*q(i,1-4+1,1) &
               +0.35352799825_wp*q(i,1-4+2,1) &
               -0.348142394842_wp*q(i,1-4+3,1) &
               +0.181481803619_wp*q(i,1-4+4,1) &
               +0.00944080437_wp*q(i,1-4+5,1) &
               -0.077675100452_wp*q(i,1-4+6,1) &
               +0.044887364863_wp*q(i,1-4+7,1) &
               -0.009971961849_wp*q(i,1-4+8,1) &
               +0.00011335942_wp*q(i,1-4+9,1) &
               +0.00011335942_wp*q(i,1-4+10,1) 
   
   q2(i,1-4+2,2) =0.052523901012_wp*q(i,1-4+0,2) &
               -0.206299133811_wp*q(i,1-4+1,2) &
               +0.35352799825_wp*q(i,1-4+2,2) &
               -0.348142394842_wp*q(i,1-4+3,2) &
               +0.181481803619_wp*q(i,1-4+4,2) &
               +0.00944080437_wp*q(i,1-4+5,2) &
               -0.077675100452_wp*q(i,1-4+6,2) &
               +0.044887364863_wp*q(i,1-4+7,2) &
               -0.009971961849_wp*q(i,1-4+8,2) &
               +0.00011335942_wp*q(i,1-4+9,2) &
               +0.00011335942_wp*q(i,1-4+10,2) 
   
   q2(i,1-4+2,3) =0.052523901012_wp*q(i,1-4+0,3) &
               -0.206299133811_wp*q(i,1-4+1,3) &
               +0.35352799825_wp*q(i,1-4+2,3) &
               -0.348142394842_wp*q(i,1-4+3,3) &
               +0.181481803619_wp*q(i,1-4+4,3) &
               +0.00944080437_wp*q(i,1-4+5,3) &
               -0.077675100452_wp*q(i,1-4+6,3) &
               +0.044887364863_wp*q(i,1-4+7,3) &
               -0.009971961849_wp*q(i,1-4+8,3) &
               +0.00011335942_wp*q(i,1-4+9,3) &
               +0.00011335942_wp*q(i,1-4+10,3) 
   
   q2(i,1-4+2,4) =0.052523901012_wp*q(i,1-4+0,4) &
               -0.206299133811_wp*q(i,1-4+1,4) &
               +0.35352799825_wp*q(i,1-4+2,4) &
               -0.348142394842_wp*q(i,1-4+3,4) &
               +0.181481803619_wp*q(i,1-4+4,4) &
               +0.00944080437_wp*q(i,1-4+5,4) &
               -0.077675100452_wp*q(i,1-4+6,4) &
               +0.044887364863_wp*q(i,1-4+7,4) &
               -0.009971961849_wp*q(i,1-4+8,4) &
               +0.00011335942_wp*q(i,1-4+9,4) &
               +0.00011335942_wp*q(i,1-4+10,4) 
   
   q2(i,1-4+2,5) =0.052523901012_wp*q(i,1-4+0,5) &
               -0.206299133811_wp*q(i,1-4+1,5) &
               +0.35352799825_wp*q(i,1-4+2,5) &
               -0.348142394842_wp*q(i,1-4+3,5) &
               +0.181481803619_wp*q(i,1-4+4,5) &
               +0.00944080437_wp*q(i,1-4+5,5) &
               -0.077675100452_wp*q(i,1-4+6,5) &
               +0.044887364863_wp*q(i,1-4+7,5) &
               -0.009971961849_wp*q(i,1-4+8,5) &
               +0.00011335942_wp*q(i,1-4+9,5) &
               +0.00011335942_wp*q(i,1-4+10,5) 
   
   q2(i,1-4+3,1) =-5.459601e-05_wp*q(i,1-4+0,1) &
               +0.042124772446_wp*q(i,1-4+1,1) &
               -0.173103107841_wp*q(i,1-4+2,1) &
               +0.299615871352_wp*q(i,1-4+3,1) &
               -0.276543612935_wp*q(i,1-4+4,1) &
               +0.131223506571_wp*q(i,1-4+5,1) &
               -0.023424966418_wp*q(i,1-4+6,1) &
               +0.013937561779_wp*q(i,1-4+7,1) &
               -0.024565095706_wp*q(i,1-4+8,1) &
               +0.013098287852_wp*q(i,1-4+9,1) &
               -0.00230862109_wp*q(i,1-4+10,1) 
   
   q2(i,1-4+3,2) =-5.459601e-05_wp*q(i,1-4+0,2) &
               +0.042124772446_wp*q(i,1-4+1,2) &
               -0.173103107841_wp*q(i,1-4+2,2) &
               +0.299615871352_wp*q(i,1-4+3,2) &
               -0.276543612935_wp*q(i,1-4+4,2) &
               +0.131223506571_wp*q(i,1-4+5,2) &
               -0.023424966418_wp*q(i,1-4+6,2) &
               +0.013937561779_wp*q(i,1-4+7,2) &
               -0.024565095706_wp*q(i,1-4+8,2) &
               +0.013098287852_wp*q(i,1-4+9,2) &
               -0.00230862109_wp*q(i,1-4+10,2) 
   
   q2(i,1-4+3,3) =-5.459601e-05_wp*q(i,1-4+0,3) &
               +0.042124772446_wp*q(i,1-4+1,3) &
               -0.173103107841_wp*q(i,1-4+2,3) &
               +0.299615871352_wp*q(i,1-4+3,3) &
               -0.276543612935_wp*q(i,1-4+4,3) &
               +0.131223506571_wp*q(i,1-4+5,3) &
               -0.023424966418_wp*q(i,1-4+6,3) &
               +0.013937561779_wp*q(i,1-4+7,3) &
               -0.024565095706_wp*q(i,1-4+8,3) &
               +0.013098287852_wp*q(i,1-4+9,3) &
               -0.00230862109_wp*q(i,1-4+10,3) 
   
   q2(i,1-4+3,4) =-5.459601e-05_wp*q(i,1-4+0,4) &
               +0.042124772446_wp*q(i,1-4+1,4) &
               -0.173103107841_wp*q(i,1-4+2,4) &
               +0.299615871352_wp*q(i,1-4+3,4) &
               -0.276543612935_wp*q(i,1-4+4,4) &
               +0.131223506571_wp*q(i,1-4+5,4) &
               -0.023424966418_wp*q(i,1-4+6,4) &
               +0.013937561779_wp*q(i,1-4+7,4) &
               -0.024565095706_wp*q(i,1-4+8,4) &
               +0.013098287852_wp*q(i,1-4+9,4) &
               -0.00230862109_wp*q(i,1-4+10,4) 
   
   q2(i,1-4+3,5) =-5.459601e-05_wp*q(i,1-4+0,5) &
               +0.042124772446_wp*q(i,1-4+1,5) &
               -0.173103107841_wp*q(i,1-4+2,5) &
               +0.299615871352_wp*q(i,1-4+3,5) &
               -0.276543612935_wp*q(i,1-4+4,5) &
               +0.131223506571_wp*q(i,1-4+5,5) &
               -0.023424966418_wp*q(i,1-4+6,5) &
               +0.013937561779_wp*q(i,1-4+7,5) &
               -0.024565095706_wp*q(i,1-4+8,5) &
               +0.013098287852_wp*q(i,1-4+9,5) &
               -0.00230862109_wp*q(i,1-4+10,5) 
   
      enddo

 CASE (6)

   idflt(2) = idarray(2)

 CASE (8)

    
         do i=idloop(1),idloop(2) 
   
   q2(i,ny+4+0,1) =0.1604411764705_wp*q(i,ny+4+0,1) &
               -0.2325_wp*q(i,ny+4-1,1) &
               +0.0895588235295_wp*q(i,ny+4-2,1) &
               -0.0175_wp*q(i,ny+4-3,1) 
   
   q2(i,ny+4+0,2) =0.1604411764705_wp*q(i,ny+4+0,2) &
               -0.2325_wp*q(i,ny+4-1,2) &
               +0.0895588235295_wp*q(i,ny+4-2,2) &
               -0.0175_wp*q(i,ny+4-3,2) 
   
   q2(i,ny+4+0,3) =0.1604411764705_wp*q(i,ny+4+0,3) &
               -0.2325_wp*q(i,ny+4-1,3) &
               +0.0895588235295_wp*q(i,ny+4-2,3) &
               -0.0175_wp*q(i,ny+4-3,3) 
   
   q2(i,ny+4+0,4) =0.1604411764705_wp*q(i,ny+4+0,4) &
               -0.2325_wp*q(i,ny+4-1,4) &
               +0.0895588235295_wp*q(i,ny+4-2,4) &
               -0.0175_wp*q(i,ny+4-3,4) 
   
   q2(i,ny+4+0,5) =0.1604411764705_wp*q(i,ny+4+0,5) &
               -0.2325_wp*q(i,ny+4-1,5) &
               +0.0895588235295_wp*q(i,ny+4-2,5) &
               -0.0175_wp*q(i,ny+4-3,5) 
   
   q2(i,ny+4-1,1) =-0.042888704485_wp*q(i,ny+4+0,1) &
               +0.138814085762_wp*q(i,ny+4-1,1) &
               -0.1784240360865_wp*q(i,ny+4-2,1) &
               +0.111559546536_wp*q(i,ny+4-3,1) &
               -0.0286735324325_wp*q(i,ny+4-4,1) &
               -0.000373632298_wp*q(i,ny+4-5,1) &
               -1.37269965e-05_wp*q(i,ny+4-6,1) 
   
   q2(i,ny+4-1,2) =-0.042888704485_wp*q(i,ny+4+0,2) &
               +0.138814085762_wp*q(i,ny+4-1,2) &
               -0.1784240360865_wp*q(i,ny+4-2,2) &
               +0.111559546536_wp*q(i,ny+4-3,2) &
               -0.0286735324325_wp*q(i,ny+4-4,2) &
               -0.000373632298_wp*q(i,ny+4-5,2) &
               -1.37269965e-05_wp*q(i,ny+4-6,2) 
   
   q2(i,ny+4-1,3) =-0.042888704485_wp*q(i,ny+4+0,3) &
               +0.138814085762_wp*q(i,ny+4-1,3) &
               -0.1784240360865_wp*q(i,ny+4-2,3) &
               +0.111559546536_wp*q(i,ny+4-3,3) &
               -0.0286735324325_wp*q(i,ny+4-4,3) &
               -0.000373632298_wp*q(i,ny+4-5,3) &
               -1.37269965e-05_wp*q(i,ny+4-6,3) 
   
   q2(i,ny+4-1,4) =-0.042888704485_wp*q(i,ny+4+0,4) &
               +0.138814085762_wp*q(i,ny+4-1,4) &
               -0.1784240360865_wp*q(i,ny+4-2,4) &
               +0.111559546536_wp*q(i,ny+4-3,4) &
               -0.0286735324325_wp*q(i,ny+4-4,4) &
               -0.000373632298_wp*q(i,ny+4-5,4) &
               -1.37269965e-05_wp*q(i,ny+4-6,4) 
   
   q2(i,ny+4-1,5) =-0.042888704485_wp*q(i,ny+4+0,5) &
               +0.138814085762_wp*q(i,ny+4-1,5) &
               -0.1784240360865_wp*q(i,ny+4-2,5) &
               +0.111559546536_wp*q(i,ny+4-3,5) &
               -0.0286735324325_wp*q(i,ny+4-4,5) &
               -0.000373632298_wp*q(i,ny+4-5,5) &
               -1.37269965e-05_wp*q(i,ny+4-6,5) 
   
   q2(i,ny+4-2,1) =0.052523901012_wp*q(i,ny+4+0,1) &
               -0.206299133811_wp*q(i,ny+4-1,1) &
               +0.35352799825_wp*q(i,ny+4-2,1) &
               -0.348142394842_wp*q(i,ny+4-3,1) &
               +0.181481803619_wp*q(i,ny+4-4,1) &
               +0.00944080437_wp*q(i,ny+4-5,1) &
               -0.077675100452_wp*q(i,ny+4-6,1) &
               +0.044887364863_wp*q(i,ny+4-7,1) &
               -0.009971961849_wp*q(i,ny+4-8,1) &
               +0.00011335942_wp*q(i,ny+4-9,1) &
               +0.00011335942_wp*q(i,ny+4-10,1) 
   
   q2(i,ny+4-2,2) =0.052523901012_wp*q(i,ny+4+0,2) &
               -0.206299133811_wp*q(i,ny+4-1,2) &
               +0.35352799825_wp*q(i,ny+4-2,2) &
               -0.348142394842_wp*q(i,ny+4-3,2) &
               +0.181481803619_wp*q(i,ny+4-4,2) &
               +0.00944080437_wp*q(i,ny+4-5,2) &
               -0.077675100452_wp*q(i,ny+4-6,2) &
               +0.044887364863_wp*q(i,ny+4-7,2) &
               -0.009971961849_wp*q(i,ny+4-8,2) &
               +0.00011335942_wp*q(i,ny+4-9,2) &
               +0.00011335942_wp*q(i,ny+4-10,2) 
   
   q2(i,ny+4-2,3) =0.052523901012_wp*q(i,ny+4+0,3) &
               -0.206299133811_wp*q(i,ny+4-1,3) &
               +0.35352799825_wp*q(i,ny+4-2,3) &
               -0.348142394842_wp*q(i,ny+4-3,3) &
               +0.181481803619_wp*q(i,ny+4-4,3) &
               +0.00944080437_wp*q(i,ny+4-5,3) &
               -0.077675100452_wp*q(i,ny+4-6,3) &
               +0.044887364863_wp*q(i,ny+4-7,3) &
               -0.009971961849_wp*q(i,ny+4-8,3) &
               +0.00011335942_wp*q(i,ny+4-9,3) &
               +0.00011335942_wp*q(i,ny+4-10,3) 
   
   q2(i,ny+4-2,4) =0.052523901012_wp*q(i,ny+4+0,4) &
               -0.206299133811_wp*q(i,ny+4-1,4) &
               +0.35352799825_wp*q(i,ny+4-2,4) &
               -0.348142394842_wp*q(i,ny+4-3,4) &
               +0.181481803619_wp*q(i,ny+4-4,4) &
               +0.00944080437_wp*q(i,ny+4-5,4) &
               -0.077675100452_wp*q(i,ny+4-6,4) &
               +0.044887364863_wp*q(i,ny+4-7,4) &
               -0.009971961849_wp*q(i,ny+4-8,4) &
               +0.00011335942_wp*q(i,ny+4-9,4) &
               +0.00011335942_wp*q(i,ny+4-10,4) 
   
   q2(i,ny+4-2,5) =0.052523901012_wp*q(i,ny+4+0,5) &
               -0.206299133811_wp*q(i,ny+4-1,5) &
               +0.35352799825_wp*q(i,ny+4-2,5) &
               -0.348142394842_wp*q(i,ny+4-3,5) &
               +0.181481803619_wp*q(i,ny+4-4,5) &
               +0.00944080437_wp*q(i,ny+4-5,5) &
               -0.077675100452_wp*q(i,ny+4-6,5) &
               +0.044887364863_wp*q(i,ny+4-7,5) &
               -0.009971961849_wp*q(i,ny+4-8,5) &
               +0.00011335942_wp*q(i,ny+4-9,5) &
               +0.00011335942_wp*q(i,ny+4-10,5) 
   
   q2(i,ny+4-3,1) =-5.459601e-05_wp*q(i,ny+4+0,1) &
               +0.042124772446_wp*q(i,ny+4-1,1) &
               -0.173103107841_wp*q(i,ny+4-2,1) &
               +0.299615871352_wp*q(i,ny+4-3,1) &
               -0.276543612935_wp*q(i,ny+4-4,1) &
               +0.131223506571_wp*q(i,ny+4-5,1) &
               -0.023424966418_wp*q(i,ny+4-6,1) &
               +0.013937561779_wp*q(i,ny+4-7,1) &
               -0.024565095706_wp*q(i,ny+4-8,1) &
               +0.013098287852_wp*q(i,ny+4-9,1) &
               -0.00230862109_wp*q(i,ny+4-10,1) 
   
   q2(i,ny+4-3,2) =-5.459601e-05_wp*q(i,ny+4+0,2) &
               +0.042124772446_wp*q(i,ny+4-1,2) &
               -0.173103107841_wp*q(i,ny+4-2,2) &
               +0.299615871352_wp*q(i,ny+4-3,2) &
               -0.276543612935_wp*q(i,ny+4-4,2) &
               +0.131223506571_wp*q(i,ny+4-5,2) &
               -0.023424966418_wp*q(i,ny+4-6,2) &
               +0.013937561779_wp*q(i,ny+4-7,2) &
               -0.024565095706_wp*q(i,ny+4-8,2) &
               +0.013098287852_wp*q(i,ny+4-9,2) &
               -0.00230862109_wp*q(i,ny+4-10,2) 
   
   q2(i,ny+4-3,3) =-5.459601e-05_wp*q(i,ny+4+0,3) &
               +0.042124772446_wp*q(i,ny+4-1,3) &
               -0.173103107841_wp*q(i,ny+4-2,3) &
               +0.299615871352_wp*q(i,ny+4-3,3) &
               -0.276543612935_wp*q(i,ny+4-4,3) &
               +0.131223506571_wp*q(i,ny+4-5,3) &
               -0.023424966418_wp*q(i,ny+4-6,3) &
               +0.013937561779_wp*q(i,ny+4-7,3) &
               -0.024565095706_wp*q(i,ny+4-8,3) &
               +0.013098287852_wp*q(i,ny+4-9,3) &
               -0.00230862109_wp*q(i,ny+4-10,3) 
   
   q2(i,ny+4-3,4) =-5.459601e-05_wp*q(i,ny+4+0,4) &
               +0.042124772446_wp*q(i,ny+4-1,4) &
               -0.173103107841_wp*q(i,ny+4-2,4) &
               +0.299615871352_wp*q(i,ny+4-3,4) &
               -0.276543612935_wp*q(i,ny+4-4,4) &
               +0.131223506571_wp*q(i,ny+4-5,4) &
               -0.023424966418_wp*q(i,ny+4-6,4) &
               +0.013937561779_wp*q(i,ny+4-7,4) &
               -0.024565095706_wp*q(i,ny+4-8,4) &
               +0.013098287852_wp*q(i,ny+4-9,4) &
               -0.00230862109_wp*q(i,ny+4-10,4) 
   
   q2(i,ny+4-3,5) =-5.459601e-05_wp*q(i,ny+4+0,5) &
               +0.042124772446_wp*q(i,ny+4-1,5) &
               -0.173103107841_wp*q(i,ny+4-2,5) &
               +0.299615871352_wp*q(i,ny+4-3,5) &
               -0.276543612935_wp*q(i,ny+4-4,5) &
               +0.131223506571_wp*q(i,ny+4-5,5) &
               -0.023424966418_wp*q(i,ny+4-6,5) &
               +0.013937561779_wp*q(i,ny+4-7,5) &
               -0.024565095706_wp*q(i,ny+4-8,5) &
               +0.013098287852_wp*q(i,ny+4-9,5) &
               -0.00230862109_wp*q(i,ny+4-10,5) 
   
      enddo

    END SELECT

enddo

!$OMP DO SCHEDULE(GUIDED,4) COLLAPSE(2) 
  do bk=idflt(5),idflt(6),size_bk  
    do bj=idflt(3),idflt(4),size_bj 
      do bi=idflt(1),idflt(2),size_bi 
    
idloop(6) = min( bk+size_bk, idflt(6)+1)-1
idloop(4) = min( bj+size_bj, idflt(4)+1)-1
idloop(2) = min( bi+size_bi, idflt(2)+1)-1

  idloop(5) = bk
  idloop(3) = bj
  idloop(1) = bi

     do j=idloop(3),idloop(4) 
 
      do i=idloop(1),idloop(2) 

q2(i,j,1) =0.0625_wp*q(i,j-2,1) &
            -0.25_wp*q(i,j-1,1) &
            +0.375_wp*q(i,j+0,1) &
            -0.25_wp*q(i,j+1,1) &
            +0.0625_wp*q(i,j+2,1) 

q2(i,j,2) =0.0625_wp*q(i,j-2,2) &
            -0.25_wp*q(i,j-1,2) &
            +0.375_wp*q(i,j+0,2) &
            -0.25_wp*q(i,j+1,2) &
            +0.0625_wp*q(i,j+2,2) 

q2(i,j,3) =0.0625_wp*q(i,j-2,3) &
            -0.25_wp*q(i,j-1,3) &
            +0.375_wp*q(i,j+0,3) &
            -0.25_wp*q(i,j+1,3) &
            +0.0625_wp*q(i,j+2,3) 

q2(i,j,4) =0.0625_wp*q(i,j-2,4) &
            -0.25_wp*q(i,j-1,4) &
            +0.375_wp*q(i,j+0,4) &
            -0.25_wp*q(i,j+1,4) &
            +0.0625_wp*q(i,j+2,4) 

q2(i,j,5) =0.0625_wp*q(i,j-2,5) &
            -0.25_wp*q(i,j-1,5) &
            +0.375_wp*q(i,j+0,5) &
            -0.25_wp*q(i,j+1,5) &
            +0.0625_wp*q(i,j+2,5) 

     
     enddo
   enddo

     enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO

!$OMP DO SCHEDULE(GUIDED,4) COLLAPSE(2) 
  do bk=idflt(5),idflt(6),size_bk  
    do bj=idflt(3),idflt(4),size_bj 
      do bi=idflt(1),idflt(2),size_bi 
    
idloop(6) = min( bk+size_bk, idflt(6)+1)-1
idloop(4) = min( bj+size_bj, idflt(4)+1)-1
idloop(2) = min( bi+size_bi, idflt(2)+1)-1

  idloop(5) = bk
  idloop(3) = bj
  idloop(1) = bi

     do j=idloop(3),idloop(4) 
 
      do i=idloop(1),idloop(2) 
q(i,j,1) = q(i,j,1) - param_float(5)*q2(i,j,1)

q(i,j,2) = q(i,j,2) - param_float(5)*q2(i,j,2)

q(i,j,3) = q(i,j,3) - param_float(5)*q2(i,j,3)

q(i,j,4) = q(i,j,4) - param_float(5)*q2(i,j,4)

q(i,j,5) = q(i,j,5) - param_float(5)*q2(i,j,5)

     enddo
   enddo
  
     enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO

idloop(1) = 1 
idloop(2) = nx
idloop(3) = 1 
idloop(4) = ny
idloop(5) = 1 
idloop(6) = nz



do ibc=1,nbc

    SELECT CASE (bc(ibc))


 CASE (1)

   idloop(1) = idarray(1)

 CASE (2)

   idloop(1) = idarray(1)

 CASE (3)

   idloop(2) = idarray(2)

 CASE (4)

   idloop(2) = idarray(2)

 CASE (7)

    
         do i=idloop(1),idloop(2) 
   
   q(i,1-4+0,1) = q(i,1-4+0,1) - param_float(5)*q2(i,1-4+0,1)
   q(i,1-4+0,2) = q(i,1-4+0,2) - param_float(5)*q2(i,1-4+0,2)
   q(i,1-4+0,3) = q(i,1-4+0,3) - param_float(5)*q2(i,1-4+0,3)
   q(i,1-4+0,4) = q(i,1-4+0,4) - param_float(5)*q2(i,1-4+0,4)
   q(i,1-4+0,5) = q(i,1-4+0,5) - param_float(5)*q2(i,1-4+0,5)
   
   q(i,1-4+1,1) = q(i,1-4+1,1) - param_float(5)*q2(i,1-4+1,1)
   q(i,1-4+1,2) = q(i,1-4+1,2) - param_float(5)*q2(i,1-4+1,2)
   q(i,1-4+1,3) = q(i,1-4+1,3) - param_float(5)*q2(i,1-4+1,3)
   q(i,1-4+1,4) = q(i,1-4+1,4) - param_float(5)*q2(i,1-4+1,4)
   q(i,1-4+1,5) = q(i,1-4+1,5) - param_float(5)*q2(i,1-4+1,5)
   
   q(i,1-4+2,1) = q(i,1-4+2,1) - param_float(5)*q2(i,1-4+2,1)
   q(i,1-4+2,2) = q(i,1-4+2,2) - param_float(5)*q2(i,1-4+2,2)
   q(i,1-4+2,3) = q(i,1-4+2,3) - param_float(5)*q2(i,1-4+2,3)
   q(i,1-4+2,4) = q(i,1-4+2,4) - param_float(5)*q2(i,1-4+2,4)
   q(i,1-4+2,5) = q(i,1-4+2,5) - param_float(5)*q2(i,1-4+2,5)
   
   q(i,1-4+3,1) = q(i,1-4+3,1) - param_float(5)*q2(i,1-4+3,1)
   q(i,1-4+3,2) = q(i,1-4+3,2) - param_float(5)*q2(i,1-4+3,2)
   q(i,1-4+3,3) = q(i,1-4+3,3) - param_float(5)*q2(i,1-4+3,3)
   q(i,1-4+3,4) = q(i,1-4+3,4) - param_float(5)*q2(i,1-4+3,4)
   q(i,1-4+3,5) = q(i,1-4+3,5) - param_float(5)*q2(i,1-4+3,5)
   
      enddo

 CASE (8)

    
         do i=idloop(1),idloop(2) 
   
   q(i,ny+4+0,1) = q(i,ny+4+0,1) - param_float(5)*q2(i,ny+4+0,1)
   q(i,ny+4+0,2) = q(i,ny+4+0,2) - param_float(5)*q2(i,ny+4+0,2)
   q(i,ny+4+0,3) = q(i,ny+4+0,3) - param_float(5)*q2(i,ny+4+0,3)
   q(i,ny+4+0,4) = q(i,ny+4+0,4) - param_float(5)*q2(i,ny+4+0,4)
   q(i,ny+4+0,5) = q(i,ny+4+0,5) - param_float(5)*q2(i,ny+4+0,5)
   
   q(i,ny+4-1,1) = q(i,ny+4-1,1) - param_float(5)*q2(i,ny+4-1,1)
   q(i,ny+4-1,2) = q(i,ny+4-1,2) - param_float(5)*q2(i,ny+4-1,2)
   q(i,ny+4-1,3) = q(i,ny+4-1,3) - param_float(5)*q2(i,ny+4-1,3)
   q(i,ny+4-1,4) = q(i,ny+4-1,4) - param_float(5)*q2(i,ny+4-1,4)
   q(i,ny+4-1,5) = q(i,ny+4-1,5) - param_float(5)*q2(i,ny+4-1,5)
   
   q(i,ny+4-2,1) = q(i,ny+4-2,1) - param_float(5)*q2(i,ny+4-2,1)
   q(i,ny+4-2,2) = q(i,ny+4-2,2) - param_float(5)*q2(i,ny+4-2,2)
   q(i,ny+4-2,3) = q(i,ny+4-2,3) - param_float(5)*q2(i,ny+4-2,3)
   q(i,ny+4-2,4) = q(i,ny+4-2,4) - param_float(5)*q2(i,ny+4-2,4)
   q(i,ny+4-2,5) = q(i,ny+4-2,5) - param_float(5)*q2(i,ny+4-2,5)
   
   q(i,ny+4-3,1) = q(i,ny+4-3,1) - param_float(5)*q2(i,ny+4-3,1)
   q(i,ny+4-3,2) = q(i,ny+4-3,2) - param_float(5)*q2(i,ny+4-3,2)
   q(i,ny+4-3,3) = q(i,ny+4-3,3) - param_float(5)*q2(i,ny+4-3,3)
   q(i,ny+4-3,4) = q(i,ny+4-3,4) - param_float(5)*q2(i,ny+4-3,4)
   q(i,ny+4-3,5) = q(i,ny+4-3,5) - param_float(5)*q2(i,ny+4-3,5)
   
      enddo

    END SELECT

enddo

!$OMP END PARALLEL


 end subroutine flt_y


 subroutine flt_z(param_float,hlo,neq,nx,ny,nz,sizeblck,nbc,bc,q,q2)

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
! LOCAL VARIABLES

  integer, intent(in) :: nx,ny,nz,hlo,neq  
  integer, intent(in) :: sizeblck(3)

  integer :: i,j,k  
  integer :: bi,bj,bk
  integer :: biend,bjend,bkend
  integer :: size_bi,size_bj,size_bk
  integer :: ibc
  integer, intent(in) :: bc(nbc),nbc

real(wp),intent(inout) ::  q(1-hlo:nx+hlo,1-hlo:ny+hlo,neq),&
                          q2(1-hlo:nx+hlo,1-hlo:ny+hlo,neq)
  
  integer :: idloop(6),idarray(6),idflt(6)


!f2py intent(in)    :: nx,ny,nz,hlo,neq,indvars,param_float
!f2py intent(inout) :: q

!$OMP PARALLEL DEFAULT(SHARED) private(idloop)

idloop(1) = 1 
idloop(2) = nx
idloop(3) = 1 
idloop(4) = ny
idloop(5) = 1 
idloop(6) = nz
  
idflt(1) = 1 
idflt(2) = nx

idflt(3) = 1 
idflt(4) = ny

idflt(5) = 1 
idflt(6) = nz

idarray(1) = 1 -hlo
idarray(2) = nx+hlo

idarray(3) = 1 -hlo
idarray(4) = ny+hlo

idarray(5) = 1 -hlo
idarray(6) = nz+hlo

size_bk = sizeblck(3)
size_bj = sizeblck(2)
size_bi = sizeblck(1)



do ibc=1,nbc

    SELECT CASE (bc(ibc))


    END SELECT

enddo

!$OMP DO SCHEDULE(GUIDED,4) COLLAPSE(2) 
  do bk=idflt(5),idflt(6),size_bk  
    do bj=idflt(3),idflt(4),size_bj 
      do bi=idflt(1),idflt(2),size_bi 
    
idloop(6) = min( bk+size_bk, idflt(6)+1)-1
idloop(4) = min( bj+size_bj, idflt(4)+1)-1
idloop(2) = min( bi+size_bi, idflt(2)+1)-1

idloop(5) = bk
idloop(3) = bj
idloop(1) = bi 


     do j=idloop(3),idloop(4) 
 
      do i=idloop(1),idloop(2) 

     
     enddo
   enddo


     enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO

!$OMP DO SCHEDULE(GUIDED,4) COLLAPSE(2) 
  do bk=idflt(5),idflt(6),size_bk  
    do bj=idflt(3),idflt(4),size_bj 
      do bi=idflt(1),idflt(2),size_bi 
    
idloop(6) = min( bk+size_bk, idflt(6)+1)-1
idloop(4) = min( bj+size_bj, idflt(4)+1)-1
idloop(2) = min( bi+size_bi, idflt(2)+1)-1

  idloop(5) = bk
  idloop(3) = bj
  idloop(1) = bi

     do j=idloop(3),idloop(4) 
 
      do i=idloop(1),idloop(2) 
     enddo
   enddo
  
     enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO

idloop(1) = 1 
idloop(2) = nx
idloop(3) = 1 
idloop(4) = ny
idloop(5) = 1 
idloop(6) = nz



do ibc=1,nbc

    SELECT CASE (bc(ibc))


    END SELECT

enddo

!$OMP END PARALLEL

 end subroutine flt_z


subroutine bc(dir,param_int,param_float,data_float)!param_float,hlo,nrk,neq,indvars,nx,ny,nz,q,q1,q2,rhs)

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
real(wp), intent(inout) ::  data_float(*)
integer,  intent(in)    ::   param_int(*)

integer, intent(in)   :: dir

!f2py  intent(in)   :: param_int,param_float,dir
!f2py intent(inout) :: data_float    

 if     (dir == 1) then

    call periodic_x(param_float                                            ,&
                    param_int(1)                                   ,&
                    param_int(5)                                  ,&
                    param_int(2),param_int(3) , param_int(4),&
                    param_int(10)                             ,&               
                    data_float(1)                                          ,&
                    data_float(1+param_int(7)*2)                   ) ! q2 is used for filter working array
    
 elseif (dir == 2)then

    call periodic_y(param_float                                            ,&
                    param_int(1)                                   ,&
                    param_int(5)                                  ,&
                    param_int(2),param_int(3) , param_int(4),&
                    param_int(10)                             ,&                              
                    data_float(1)                                          ,&
                    data_float(1+param_int(7)*2)                   ) ! q2 is used for filter working array

 elseif (dir == 3)then


    call periodic_z(param_float                                            ,&
                    param_int(1)                                   ,&
                    param_int(5)                                  ,&
                    param_int(2),param_int(3) , param_int(4),&
                    param_int(10)                             ,&                              
                    data_float(1)                                          ,&
                    data_float(1+param_int(7)*2)                   ) ! q2 is used for filter working array

 elseif (dir == 0 )then 

  call periodic_x(param_float                                            ,&
                    param_int(1)                                   ,&
                    param_int(5)                                  ,&
                    param_int(2),param_int(3) , param_int(4),&
                    param_int(10)                             ,&               
                    data_float(1)                                          ,&
                    data_float(1+param_int(7)*2)                   ) ! q2 is used for filter working array
  call periodic_y(param_float                                            ,&
                    param_int(1)                                   ,&
                    param_int(5)                                  ,&
                    param_int(2),param_int(3) , param_int(4),&
                    param_int(10)                             ,&                              
                    data_float(1)                                          ,&
                    data_float(1+param_int(7)*2)                   ) ! q2 is used for filter working array
  call periodic_z(param_float                                            ,&
                    param_int(1)                                   ,&
                    param_int(5)                                  ,&
                    param_int(2),param_int(3) , param_int(4),&
                    param_int(10)                             ,&                              
                    data_float(1)                                          ,&
                    data_float(1+param_int(7)*2)                   ) ! q2 is used for filter working array

    
 endif



 end subroutine bc

subroutine periodic_x(param_float,hlo,neq,nx,ny,nz,sizeblck,q,q2)

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
! LOCAL VARIABLES

  integer, intent(in) :: nx,ny,nz,hlo,neq  
  integer, intent(in) :: sizeblck(3)

  integer :: i,j,k  
  integer :: bi,bj,bk
  integer :: biend,bjend,bkend
  integer :: size_bi,size_bj,size_bk

real(wp),intent(inout) ::  q(1-hlo:nx+hlo,1-hlo:ny+hlo,neq),&
                          q2(1-hlo:nx+hlo,1-hlo:ny+hlo,neq)
  
  integer :: idloop(6),idarray(6)


!f2py intent(in)    :: nx,ny,nz,hlo,neq,indvars,param_float
!f2py intent(inout) :: q

!$OMP PARALLEL DEFAULT(SHARED) private(idloop)


  idarray(1) = 1 -hlo
  idarray(2) = nx+hlo

  idarray(3) = 1 -hlo
  idarray(4) = ny+hlo

  idarray(5) = 1 -hlo
  idarray(6) = nz+hlo

size_bk = sizeblck(3)
size_bj = sizeblck(2)
size_bi = sizeblck(1)

!$OMP DO SCHEDULE(GUIDED,4) COLLAPSE(2) 
  do bk=1,nz,size_bk  
    do bj=1,ny,size_bj 
    
idloop(6) = min( bk+size_bk, nz +1)-1
idloop(4) = min( bj+size_bj, ny +1)-1
idloop(2) = 0

idloop(5) = bk
idloop(3) = bj
idloop(1) = 0 

     do j=idloop(3),idloop(4) 
 
      do i=idloop(1),idloop(2) 

q(0+0,j,1) = q(nx+0,j,1)
q(0-1,j,1) = q(nx-1,j,1)

q(nx+1+0,j,1) = q(1+0,j,1)
q(nx+1+1,j,1) = q(1+1,j,1)

q(0+0,j,2) = q(nx+0,j,2)
q(0-1,j,2) = q(nx-1,j,2)

q(nx+1+0,j,2) = q(1+0,j,2)
q(nx+1+1,j,2) = q(1+1,j,2)

q(0+0,j,3) = q(nx+0,j,3)
q(0-1,j,3) = q(nx-1,j,3)

q(nx+1+0,j,3) = q(1+0,j,3)
q(nx+1+1,j,3) = q(1+1,j,3)

q(0+0,j,4) = q(nx+0,j,4)
q(0-1,j,4) = q(nx-1,j,4)

q(nx+1+0,j,4) = q(1+0,j,4)
q(nx+1+1,j,4) = q(1+1,j,4)

q(0+0,j,5) = q(nx+0,j,5)
q(0-1,j,5) = q(nx-1,j,5)

q(nx+1+0,j,5) = q(1+0,j,5)
q(nx+1+1,j,5) = q(1+1,j,5)

     
     enddo
   enddo

  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO
!$OMP END PARALLEL
 
 end subroutine periodic_x

subroutine periodic_y(param_float,hlo,neq,nx,ny,nz,sizeblck,q,q2)

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
! LOCAL VARIABLES

  integer, intent(in) :: nx,ny,nz,hlo,neq  
  integer, intent(in) :: sizeblck(3)

  integer :: i,j,k  
  integer :: bi,bj,bk
  integer :: biend,bjend,bkend
  integer :: size_bi,size_bj,size_bk

real(wp),intent(inout) ::  q(1-hlo:nx+hlo,1-hlo:ny+hlo,neq),&
                          q2(1-hlo:nx+hlo,1-hlo:ny+hlo,neq)
  
  integer :: idloop(6),idarray(6)


!f2py intent(in)    :: nx,ny,nz,hlo,neq,indvars,param_float
!f2py intent(inout) :: q

!$OMP PARALLEL DEFAULT(SHARED) private(idloop)

 idarray(1) = 1 -hlo
 idarray(2) = nx+hlo

 idarray(3) = 1 -hlo
 idarray(4) = ny+hlo

 idarray(5) = 1 -hlo
 idarray(6) = nz+hlo


size_bk = sizeblck(3)
size_bj = sizeblck(2)
size_bi = sizeblck(1)

!$OMP DO SCHEDULE(GUIDED,4) COLLAPSE(2) 
  do bk=1,nz,size_bk      
      do bi=1-hlo,nx+hlo,size_bi 
    
  idloop(6) = min( bk+size_bk, nz+1)-1
  idloop(4) = 0
  idloop(2) = min( bi+size_bi, nx+hlo+1)-1

  idloop(5) = bk
  idloop(3) = 0
  idloop(1) = bi

     do j=idloop(3),idloop(4) 
 
      do i=idloop(1),idloop(2) 

q(i,0+0,1) = q(i,ny+0,1)
q(i,0-1,1) = q(i,ny-1,1)

q(i,ny+1+0,1) = q(i,1+0,1)
q(i,ny+1+1,1) = q(i,1+1,1)

q(i,0+0,2) = q(i,ny+0,2)
q(i,0-1,2) = q(i,ny-1,2)

q(i,ny+1+0,2) = q(i,1+0,2)
q(i,ny+1+1,2) = q(i,1+1,2)

q(i,0+0,3) = q(i,ny+0,3)
q(i,0-1,3) = q(i,ny-1,3)

q(i,ny+1+0,3) = q(i,1+0,3)
q(i,ny+1+1,3) = q(i,1+1,3)

q(i,0+0,4) = q(i,ny+0,4)
q(i,0-1,4) = q(i,ny-1,4)

q(i,ny+1+0,4) = q(i,1+0,4)
q(i,ny+1+1,4) = q(i,1+1,4)

q(i,0+0,5) = q(i,ny+0,5)
q(i,0-1,5) = q(i,ny-1,5)

q(i,ny+1+0,5) = q(i,1+0,5)
q(i,ny+1+1,5) = q(i,1+1,5)

     
     enddo
   enddo

     enddo ! END cache blocking i
enddo ! END cache blocking k
!$OMP END DO
!$OMP END PARALLEL


 end subroutine periodic_y


 subroutine periodic_z(param_float,hlo,neq,nx,ny,nz,sizeblck,q,q2)

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
! LOCAL VARIABLES

  integer, intent(in) :: nx,ny,nz,hlo,neq  
  integer, intent(in) :: sizeblck(3)

  integer :: i,j,k  
  integer :: bi,bj,bk
  integer :: biend,bjend,bkend
  integer :: size_bi,size_bj,size_bk

real(wp),intent(inout) ::  q(1-hlo:nx+hlo,1-hlo:ny+hlo,neq),&
                          q2(1-hlo:nx+hlo,1-hlo:ny+hlo,neq)
  
  integer :: idloop(6),idarray(6)


!f2py intent(in)    :: nx,ny,nz,hlo,neq,indvars,param_float
!f2py intent(inout) :: q

!$OMP PARALLEL DEFAULT(SHARED) private(idloop)


  idarray(1) = 1 -hlo
  idarray(2) = nx+hlo

  idarray(3) = 1 -hlo
  idarray(4) = ny+hlo

  idarray(5) = 1 -hlo
  idarray(6) = nz+hlo

size_bk = sizeblck(3)
size_bj = sizeblck(2)
size_bi = sizeblck(1)

!$OMP DO SCHEDULE(GUIDED,4) COLLAPSE(2)  
    do bj=1-hlo,ny+hlo,size_bj 
      do bi=1-hlo,nx+hlo,size_bi 
    
idloop(6) = 0
idloop(4) = min( bj+size_bj, ny+hlo+1)-1
idloop(2) = min( bi+size_bi, nx+hlo+1)-1

idloop(5) = 0
idloop(3) = bj
idloop(1) = bi 


     do j=idloop(3),idloop(4) 
 
      do i=idloop(1),idloop(2) 





















     
     enddo
   enddo


     enddo ! END cache blocking i
  enddo ! END cache blocking j
!$OMP END DO
!$OMP END PARALLEL

 end subroutine periodic_z

subroutine pack(buf,f,ibeg,iend,jbeg,jend,kbeg,kend,sizex,sizey,sizez,sizenv)

implicit none

integer,parameter :: wp = kind(0.0D0) ! working precision

real(wp), intent(inout) :: buf(*),f(*)
integer,intent(in)      :: ibeg,iend,jbeg,jend,kbeg,kend,sizex,sizey,sizez,sizenv

integer :: sx,sy,sz,i,j,k,n,indbuf,indf

sx = iend - ibeg 
sy = jend - jbeg 
sz = kend - kbeg 


!$OMP PARALLEL DEFAULT(SHARED) private(indbuf,indf)

!$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
do n=1,sizenv
 do k=kbeg,kend-1
    do j=jbeg,jend-1 
!$omp simd 
      do i=ibeg,iend-1 
        indbuf = i - ibeg + 1 + (j-jbeg)*sx + (k-kbeg)*sx*sy + (n-1)*sx*sy*sz
        indf   = i        + 1 + j*sizex     + k*sizex*sizey  + (n-1)*sizex*sizey*sizez
        buf(indbuf) = f(indf)
      enddo
    enddo
 enddo
enddo 
!$OMP END DO NOWAIT
!$OMP END PARALLEL

end subroutine pack  


subroutine unpack(buf,f,ibeg,iend,jbeg,jend,kbeg,kend,sizex,sizey,sizez,sizenv)

implicit none

integer,parameter :: wp = kind(0.0D0) ! working precision

real(wp), intent(inout) :: buf(*),f(*)
integer,intent(in)      :: ibeg,iend,jbeg,jend,kbeg,kend,sizex,sizey,sizez,sizenv

integer :: sx,sy,sz,i,j,k,n,indbuf,indf

sx = iend - ibeg 
sy = jend - jbeg 
sz = kend - kbeg 


!$OMP PARALLEL DEFAULT(SHARED) private(indbuf,indf)

!$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
do n=1,sizenv
 do k=kbeg,kend-1
    do j=jbeg,jend-1 
!$omp simd       
      do i=ibeg,iend-1 
        indbuf = i - ibeg + 1 + (j-jbeg)*sx + (k-kbeg)*sx*sy + (n-1)*sx*sy*sz
        indf   = i        + 1 + j*sizex     + k*sizex*sizey  + (n-1)*sizex*sizey*sizez
        f(indf) = buf(indbuf)
      enddo
    enddo
 enddo
enddo 
!$OMP END DO NOWAIT
!$OMP END PARALLEL

end subroutine unpack  

 subroutine init(param_int,param_float,data_float)

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
  real(wp), intent(inout) ::  data_float(*)
  integer,  intent(in)    :: param_int(*)

  call init_numa(param_float                                            ,&
                 param_int(1), param_int(8)                 ,&
                 param_int(5)                                    ,&
                 param_int(6)                                  ,&
                 param_int(13)                                    ,&
                 param_int(2),param_int(3) , param_int(4),&
                 param_int(10)                               ,&
                 data_float(1)                                          ,& ! q
                 data_float(1+param_int(7)  )                 ,& ! q1
                 data_float(1+param_int(7)*2)                 ,& ! q2 
                 data_float(1+param_int(7)*3)                 ,& ! rhs
                 data_float(1+param_int(7)*4)                  ) ! stored (if any)        


  end subroutine init

subroutine init_numa(param_float,hlo,nrk,neq,neqst,ind,nx,ny,nz,sizeblck,q,q1,q2,rhs,qst)

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


  
  integer :: i,j,k
  integer :: bi,bj,bk
  integer :: biend,bjend,bkend
  integer :: size_bi,size_bj,size_bk

  real(wp), intent(in)    :: param_float(*)

  integer, intent(in) :: nx,ny,nz,hlo,nrk,neq,neqst
  integer, intent(in) :: ind(1:neq+neqst)
  integer, intent(in) :: sizeblck(3)
  
real(wp),intent(inout) ::  q(1-hlo:nx+hlo,1-hlo:ny+hlo,neq),&
                          q1(1-hlo:nx+hlo,1-hlo:ny+hlo,neq),&
                          q2(1-hlo:nx+hlo,1-hlo:ny+hlo,neq),&
                         rhs(1-hlo:nx+hlo,1-hlo:ny+hlo,neq),&
                        qst(1-hlo:nx+hlo,1-hlo:ny+hlo,neqst)

  
  integer :: idloop(6),idarray(6),indvars(neq),indvarsst(neqst)


  indvars   = ind(1:neq)
  indvarsst = ind(1+neq:neq+neqst) 

!f2py intent(in)    :: q1,q2,a,nx,ny,nz,rhs,nrk
!f2py intent(inout) :: q

!$OMP PARALLEL DEFAULT(SHARED) private(idloop)

  idarray(1) = 1 -hlo
  idarray(2) = nx+hlo

  idarray(3) = 1 -hlo
  idarray(4) = ny+hlo

  idarray(5) = 1 -hlo
  idarray(6) = nz+hlo

size_bk = sizeblck(3)
size_bj = sizeblck(2)
size_bi = sizeblck(1)

!$OMP DO SCHEDULE(STATIC) 
  do bk=1,nz,size_bk  
    do bj=1,ny,size_bj 
      do bi=1,nx,size_bi 
    
idloop(6) = min( bk+size_bk, nz+1)-1
idloop(4) = min( bj+size_bj, ny+1)-1
idloop(2) = min( bi+size_bi, nx+1)-1

idloop(5) = bk
idloop(3) = bj
idloop(1) = bi 


     do j=idloop(3),idloop(4) 
 
      do i=idloop(1),idloop(2) 

q(i,j,1) = 0.0_wp
q1(i,j,1) = 0.0_wp
q2(i,j,1) = 0.0_wp
rhs(i,j,1) = 0.0_wp
q(i,j,2) = 0.0_wp
q1(i,j,2) = 0.0_wp
q2(i,j,2) = 0.0_wp
rhs(i,j,2) = 0.0_wp
q(i,j,3) = 0.0_wp
q1(i,j,3) = 0.0_wp
q2(i,j,3) = 0.0_wp
rhs(i,j,3) = 0.0_wp
q(i,j,4) = 0.0_wp
q1(i,j,4) = 0.0_wp
q2(i,j,4) = 0.0_wp
rhs(i,j,4) = 0.0_wp
q(i,j,5) = 0.0_wp
q1(i,j,5) = 0.0_wp
q2(i,j,5) = 0.0_wp
rhs(i,j,5) = 0.0_wp

     enddo
   enddo

    enddo ! END cache blocking i
  enddo ! END cache blocking j
enddo ! END cache blocking k
!$OMP END DO
!$OMP END PARALLEL

end subroutine init_numa
