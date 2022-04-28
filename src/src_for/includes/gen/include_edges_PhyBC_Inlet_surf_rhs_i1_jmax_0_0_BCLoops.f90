

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
! building source terms in RHS for layer 0 0 None d(rho)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (1.0_wp/c**2*((((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI)*(gamma_m1)+0.5_wp*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI+((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI)))+[rho*v]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_rhs_rho_dx_0_1m2p0p0nyp2p0k = q(1-2+0+0,ny+2+0,indvars(2))

d1_rhs_rho_dx_0_1m2p0p1nyp2p0k = q(1-2+0+1,ny+2+0,indvars(2))

d1_rhs_rho_dx_0_1m2p0p2nyp2p0k = q(1-2+0+2,ny+2+0,indvars(2))

d1_rhs_rho_dx_0_1m2p0nyp2p0k = -&
          1.5_wp*d1_rhs_rho_dx_0_1m2p0p0nyp2p0k+&
          2.0_wp*d1_rhs_rho_dx_0_1m2p0p1nyp2p0k-&
          0.5_wp*d1_rhs_rho_dx_0_1m2p0p2nyp2p0k

d1_rhs_rho_dx_0_1m2p0nyp2p0k = d1_rhs_rho_dx_0_1m2p0nyp2p0k*param_float(1)

d1_rhs_rho_dx_1_1m2p0p0nyp2p0k = (param_float(3 + 5))*q(1-2+0+0,ny+2+0,indvars(1))*((q(1-2+0+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+0+0,ny+2+0,indvars(2))*q(1-2+0+0,ny+2+0,indvars(2))+&
                    q(1-2+0+0,ny+2+0,indvars(3))*q(1-2+0+0,ny+2+0,indvars(3)))))

d1_rhs_rho_dx_1_1m2p0p1nyp2p0k = (param_float(3 + 5))*q(1-2+0+1,ny+2+0,indvars(1))*((q(1-2+0+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+0+1,ny+2+0,indvars(2))*q(1-2+0+1,ny+2+0,indvars(2))+&
                    q(1-2+0+1,ny+2+0,indvars(3))*q(1-2+0+1,ny+2+0,indvars(3)))))

d1_rhs_rho_dx_1_1m2p0p2nyp2p0k = (param_float(3 + 5))*q(1-2+0+2,ny+2+0,indvars(1))*((q(1-2+0+2,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+0+2,ny+2+0,indvars(2))*q(1-2+0+2,ny+2+0,indvars(2))+&
                    q(1-2+0+2,ny+2+0,indvars(3))*q(1-2+0+2,ny+2+0,indvars(3)))))

d1_rhs_rho_dx_1_1m2p0nyp2p0k = -&
          1.5_wp*d1_rhs_rho_dx_1_1m2p0p0nyp2p0k+&
          2.0_wp*d1_rhs_rho_dx_1_1m2p0p1nyp2p0k-&
          0.5_wp*d1_rhs_rho_dx_1_1m2p0p2nyp2p0k

d1_rhs_rho_dx_1_1m2p0nyp2p0k = d1_rhs_rho_dx_1_1m2p0nyp2p0k*param_float(1)

d1_rhs_rho_dy_0_1m2p0nyp2p0p0k = q(1-2+0,ny+2+0+0,indvars(1))*q(1-2+0,ny+2+0+0,indvars(3))

d1_rhs_rho_dy_0_1m2p0nyp2p0m1k = q(1-2+0,ny+2+0-1,indvars(1))*q(1-2+0,ny+2+0-1,indvars(3))

d1_rhs_rho_dy_0_1m2p0nyp2p0m2k = q(1-2+0,ny+2+0-2,indvars(1))*q(1-2+0,ny+2+0-2,indvars(3))

d1_rhs_rho_dy_0_1m2p0nyp2p0k = 1.5_wp*d1_rhs_rho_dy_0_1m2p0nyp2p0p0k-&
          2.0_wp*d1_rhs_rho_dy_0_1m2p0nyp2p0m1k+&
          0.5_wp*d1_rhs_rho_dy_0_1m2p0nyp2p0m2k

d1_rhs_rho_dy_0_1m2p0nyp2p0k = d1_rhs_rho_dy_0_1m2p0nyp2p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None d(rho)/dt *************
!                                                           
!***********************************************************


rhs(1-2+0,ny+2+0,indvars(1)) =   -  ( (1.0_wp/(param_float(23 + 5)*(param_float(3 + 5))*q(1-2+0,ny+2+0,indvars(1))*((q(1-2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+0,ny+2+0,indvars(2))*q(1-2+0,ny+2+0,indvars(2))+&
                    q(1-2+0,ny+2+0,indvars(3))*q(1-2+0,ny+2+0,indvars(3)))))/q(1-2+0,ny+2+0,indvars(1)))**0.5**2*(((((param_float(23 + 5)*(param_float(3 + 5))*q(1-2+0,ny+2+0,indvars(1))*((q(1-2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+0,ny+2+0,indvars(2))*q(1-2+0,ny+2+0,indvars(2))+&
                    q(1-2+0,ny+2+0,indvars(3))*q(1-2+0,ny+2+0,indvars(3)))))/q(1-2+0,ny+2+0,indvars(1)))**0.5+&
                    q(1-2+0,ny+2+0,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(1-2+0,ny+2+0,indvars(1))*((q(1-2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+0,ny+2+0,indvars(2))*q(1-2+0,ny+2+0,indvars(2))+&
                    q(1-2+0,ny+2+0,indvars(3))*q(1-2+0,ny+2+0,indvars(3)))))/q(1-2+0,ny+2+0,indvars(1)))**0.5*q(1-2+0,ny+2+0,indvars(1))*d1_rhs_rho_dx_0_1m2p0nyp2p0k+&
                    1.0_wp*d1_rhs_rho_dx_1_1m2p0nyp2p0k))*qst(1-2+0,ny+2+0,indvarsst(10)))*(param_float(3 + 5))+&
                    0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(1-2+0,ny+2+0,indvars(1))*((q(1-2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+0,ny+2+0,indvars(2))*q(1-2+0,ny+2+0,indvars(2))+&
                    q(1-2+0,ny+2+0,indvars(3))*q(1-2+0,ny+2+0,indvars(3)))))/q(1-2+0,ny+2+0,indvars(1)))**0.5+&
                    q(1-2+0,ny+2+0,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(1-2+0,ny+2+0,indvars(1))*((q(1-2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+0,ny+2+0,indvars(2))*q(1-2+0,ny+2+0,indvars(2))+&
                    q(1-2+0,ny+2+0,indvars(3))*q(1-2+0,ny+2+0,indvars(3)))))/q(1-2+0,ny+2+0,indvars(1)))**0.5*q(1-2+0,ny+2+0,indvars(1))*d1_rhs_rho_dx_0_1m2p0nyp2p0k+&
                    1.0_wp*d1_rhs_rho_dx_1_1m2p0nyp2p0k))*qst(1-2+0,ny+2+0,indvarsst(10))+&
                    (((param_float(23 + 5)*(param_float(3 + 5))*q(1-2+0,ny+2+0,indvars(1))*((q(1-2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+0,ny+2+0,indvars(2))*q(1-2+0,ny+2+0,indvars(2))+&
                    q(1-2+0,ny+2+0,indvars(3))*q(1-2+0,ny+2+0,indvars(3)))))/q(1-2+0,ny+2+0,indvars(1)))**0.5+&
                    q(1-2+0,ny+2+0,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(1-2+0,ny+2+0,indvars(1))*((q(1-2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(1-2+0,ny+2+0,indvars(2))*q(1-2+0,ny+2+0,indvars(2))+&
                    q(1-2+0,ny+2+0,indvars(3))*q(1-2+0,ny+2+0,indvars(3)))))/q(1-2+0,ny+2+0,indvars(1)))**0.5*q(1-2+0,ny+2+0,indvars(1))*d1_rhs_rho_dx_0_1m2p0nyp2p0k+&
                    1.0_wp*d1_rhs_rho_dx_1_1m2p0nyp2p0k))*qst(1-2+0,ny+2+0,indvarsst(10)))))+d1_rhs_rho_dy_0_1m2p0nyp2p0k ) 

