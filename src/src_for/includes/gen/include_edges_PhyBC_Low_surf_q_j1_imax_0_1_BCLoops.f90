

!***********************************************************
!                                                           
! Start building layers for BC : imax j1 None **************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 1 0 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 1 0 None d(rho u)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! symm*u
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 0 None d(rho u)/dt ***********
!                                                           
!***********************************************************


q(nx+2-1,1-2+0,indvars(2)) =  qst(nx+2-1,1-2+0,indvarsst(5))*q(nx+2-1,1-2+0,indvars(2))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 0 None d(rho v)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! 0.0_wp
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 0 None d(rho v)/dt ***********
!                                                           
!***********************************************************


q(nx+2-1,1-2+0,indvars(3)) =  0.0_wp



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 0 None d(rho nut)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! symm*nut
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 0 None d(rho nut)/dt *********
!                                                           
!***********************************************************


q(nx+2-1,1-2+0,indvars(5)) =  qst(nx+2-1,1-2+0,indvarsst(5))*q(nx+2-1,1-2+0,indvars(5))

