

!***********************************************************
!                                                           
! Start building layers for BC : i1 j1 None ****************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 0 1 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None d(rho)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (1.0_wp/c**2*((((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI)*(gamma_m1)+0.5_wp*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI+((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI)))+[rho*v]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_rhs_rho_dx_0_1m2p0p01m2p1k = q(1-2+0+0,1-2+1,indvars(2))

d1_rhs_rho_dx_0_1m2p0p11m2p1k = q(1-2+0+1,1-2+1,indvars(2))

d1_rhs_rho_dx_0_1m2p0p21m2p1k = q(1-2+0+2,1-2+1,indvars(2))

d1_rhs_rho_dx_0_1m2p01m2p1k = -&
          1.5_wp*d1_rhs_rho_dx_0_1m2p0p01m2p1k+&
          2.0_wp*d1_rhs_rho_dx_0_1m2p0p11m2p1k-&
          0.5_wp*d1_rhs_rho_dx_0_1m2p0p21m2p1k

d1_rhs_rho_dx_0_1m2p01m2p1k = d1_rhs_rho_dx_0_1m2p01m2p1k*param_float(1)

d1_rhs_rho_dx_1_1m2p0p01m2p1k = (param_float(3 + 5))*q(1-2+0+0,1-2+1,indvars(1))*((q(1-2+0+0,1-2+1,indvars(4))-&
                    0.5_wp*(q(1-2+0+0,1-2+1,indvars(2))*q(1-2+0+0,1-2+1,indvars(2))+&
                    q(1-2+0+0,1-2+1,indvars(3))*q(1-2+0+0,1-2+1,indvars(3)))))

d1_rhs_rho_dx_1_1m2p0p11m2p1k = (param_float(3 + 5))*q(1-2+0+1,1-2+1,indvars(1))*((q(1-2+0+1,1-2+1,indvars(4))-&
                    0.5_wp*(q(1-2+0+1,1-2+1,indvars(2))*q(1-2+0+1,1-2+1,indvars(2))+&
                    q(1-2+0+1,1-2+1,indvars(3))*q(1-2+0+1,1-2+1,indvars(3)))))

d1_rhs_rho_dx_1_1m2p0p21m2p1k = (param_float(3 + 5))*q(1-2+0+2,1-2+1,indvars(1))*((q(1-2+0+2,1-2+1,indvars(4))-&
                    0.5_wp*(q(1-2+0+2,1-2+1,indvars(2))*q(1-2+0+2,1-2+1,indvars(2))+&
                    q(1-2+0+2,1-2+1,indvars(3))*q(1-2+0+2,1-2+1,indvars(3)))))

d1_rhs_rho_dx_1_1m2p01m2p1k = -&
          1.5_wp*d1_rhs_rho_dx_1_1m2p0p01m2p1k+&
          2.0_wp*d1_rhs_rho_dx_1_1m2p0p11m2p1k-&
          0.5_wp*d1_rhs_rho_dx_1_1m2p0p21m2p1k

d1_rhs_rho_dx_1_1m2p01m2p1k = d1_rhs_rho_dx_1_1m2p01m2p1k*param_float(1)

d1_rhs_rho_dy_0_1m2p01m2p1m1k = q(1-2+0,1-2+1-1,indvars(1))*q(1-2+0,1-2+1-1,indvars(3))

d1_rhs_rho_dy_0_1m2p01m2p1p1k = q(1-2+0,1-2+1+1,indvars(1))*q(1-2+0,1-2+1+1,indvars(3))

d1_rhs_rho_dy_0_1m2p01m2p1k = -&
          0.5_wp*d1_rhs_rho_dy_0_1m2p01m2p1m1k+&
          0.5_wp*d1_rhs_rho_dy_0_1m2p01m2p1p1k

d1_rhs_rho_dy_0_1m2p01m2p1k = d1_rhs_rho_dy_0_1m2p01m2p1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None d(rho)/dt *************
!                                                           
!***********************************************************


rhs(1-2+0,1-2+1,indvars(1)) =   -  ( (1.0_wp/(param_float(23 + 5)*(param_float(3 + 5))*q(1-2+0,1-2+1,indvars(1))*((q(1-2+0,1-2+1,indvars(4))-&
                    0.5_wp*(q(1-2+0,1-2+1,indvars(2))*q(1-2+0,1-2+1,indvars(2))+&
                    q(1-2+0,1-2+1,indvars(3))*q(1-2+0,1-2+1,indvars(3)))))/q(1-2+0,1-2+1,indvars(1)))**0.5**2*(((((param_float(23 + 5)*(param_float(3 + 5))*q(1-2+0,1-2+1,indvars(1))*((q(1-2+0,1-2+1,indvars(4))-&
                    0.5_wp*(q(1-2+0,1-2+1,indvars(2))*q(1-2+0,1-2+1,indvars(2))+&
                    q(1-2+0,1-2+1,indvars(3))*q(1-2+0,1-2+1,indvars(3)))))/q(1-2+0,1-2+1,indvars(1)))**0.5+&
                    q(1-2+0,1-2+1,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(1-2+0,1-2+1,indvars(1))*((q(1-2+0,1-2+1,indvars(4))-&
                    0.5_wp*(q(1-2+0,1-2+1,indvars(2))*q(1-2+0,1-2+1,indvars(2))+&
                    q(1-2+0,1-2+1,indvars(3))*q(1-2+0,1-2+1,indvars(3)))))/q(1-2+0,1-2+1,indvars(1)))**0.5*q(1-2+0,1-2+1,indvars(1))*d1_rhs_rho_dx_0_1m2p01m2p1k+&
                    1.0_wp*d1_rhs_rho_dx_1_1m2p01m2p1k))*qst(1-2+0,1-2+1,indvarsst(10)))*(param_float(3 + 5))+&
                    0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(1-2+0,1-2+1,indvars(1))*((q(1-2+0,1-2+1,indvars(4))-&
                    0.5_wp*(q(1-2+0,1-2+1,indvars(2))*q(1-2+0,1-2+1,indvars(2))+&
                    q(1-2+0,1-2+1,indvars(3))*q(1-2+0,1-2+1,indvars(3)))))/q(1-2+0,1-2+1,indvars(1)))**0.5+&
                    q(1-2+0,1-2+1,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(1-2+0,1-2+1,indvars(1))*((q(1-2+0,1-2+1,indvars(4))-&
                    0.5_wp*(q(1-2+0,1-2+1,indvars(2))*q(1-2+0,1-2+1,indvars(2))+&
                    q(1-2+0,1-2+1,indvars(3))*q(1-2+0,1-2+1,indvars(3)))))/q(1-2+0,1-2+1,indvars(1)))**0.5*q(1-2+0,1-2+1,indvars(1))*d1_rhs_rho_dx_0_1m2p01m2p1k+&
                    1.0_wp*d1_rhs_rho_dx_1_1m2p01m2p1k))*qst(1-2+0,1-2+1,indvarsst(10))+&
                    (((param_float(23 + 5)*(param_float(3 + 5))*q(1-2+0,1-2+1,indvars(1))*((q(1-2+0,1-2+1,indvars(4))-&
                    0.5_wp*(q(1-2+0,1-2+1,indvars(2))*q(1-2+0,1-2+1,indvars(2))+&
                    q(1-2+0,1-2+1,indvars(3))*q(1-2+0,1-2+1,indvars(3)))))/q(1-2+0,1-2+1,indvars(1)))**0.5+&
                    q(1-2+0,1-2+1,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(1-2+0,1-2+1,indvars(1))*((q(1-2+0,1-2+1,indvars(4))-&
                    0.5_wp*(q(1-2+0,1-2+1,indvars(2))*q(1-2+0,1-2+1,indvars(2))+&
                    q(1-2+0,1-2+1,indvars(3))*q(1-2+0,1-2+1,indvars(3)))))/q(1-2+0,1-2+1,indvars(1)))**0.5*q(1-2+0,1-2+1,indvars(1))*d1_rhs_rho_dx_0_1m2p01m2p1k+&
                    1.0_wp*d1_rhs_rho_dx_1_1m2p01m2p1k))*qst(1-2+0,1-2+1,indvarsst(10)))))+d1_rhs_rho_dy_0_1m2p01m2p1k ) 

