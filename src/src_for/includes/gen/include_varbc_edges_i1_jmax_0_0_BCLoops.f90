

!***********************************************************
!                                                           
! Start building layers for BC : i1 jmax None **************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 0 0 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None M_jmax ***
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (u/c)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None M_jmax ****************
!                                                           
!***********************************************************


qface_j(1-2+0,1) =  (q(1-2+0,ny+2+0,indvars(2))/(param_float(23 + 5)*(param_float(3 + 5))*q(1-2+0,ny+2+0,indvars(1))*((q(1-2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+0,ny+2+0,indvars(2))*q(1-2+0,ny+2+0,indvars(2))+&
                    q(1-2+0,ny+2+0,indvars(3))*q(1-2+0,ny+2+0,indvars(3)))))/q(1-2+0,ny+2+0,indvars(1)))**0.5)

