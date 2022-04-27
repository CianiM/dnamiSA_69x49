

!***********************************************************
!                                                           
! Start building layers for BC : imax jmax None ************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 2 0 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 2 0 None M_jmax ***
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (u/c)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 2 0 None M_jmax ****************
!                                                           
!***********************************************************


qface_j(nx+4-2,1) =  (q(nx+4-2,ny+4+0,indvars(2))/(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4-2,ny+4+0,indvars(1))*((q(nx+4-2,ny+4+0,indvars(4))-&
                    0.5_wp*(q(nx+4-2,ny+4+0,indvars(2))*q(nx+4-2,ny+4+0,indvars(2))+&
                    q(nx+4-2,ny+4+0,indvars(3))*q(nx+4-2,ny+4+0,indvars(3)))))/q(nx+4-2,ny+4+0,indvars(1)))**0.5)

