

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
! building source terms in RHS for layer 2 0 None stemp ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (((0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))**2+(0.5_wp*(deltaxI*([v]_1x)-deltayI*([u]_1y)))**2)*2)**0.5
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_nxp4m2m1nyp4p0k = q(nx+4-2-1,ny+4+0,indvars(3))

d1_stemp_dx_0_nxp4m2p1nyp4p0k = q(nx+4-2+1,ny+4+0,indvars(3))

d1_stemp_dx_0_nxp4m2nyp4p0k = -&
          0.5_wp*d1_stemp_dx_0_nxp4m2m1nyp4p0k+&
          0.5_wp*d1_stemp_dx_0_nxp4m2p1nyp4p0k

d1_stemp_dx_0_nxp4m2nyp4p0k = d1_stemp_dx_0_nxp4m2nyp4p0k*param_float(1)

d1_stemp_dy_0_nxp4m2nyp4p0p0k = q(nx+4-2,ny+4+0+0,indvars(2))

d1_stemp_dy_0_nxp4m2nyp4p0m1k = q(nx+4-2,ny+4+0-1,indvars(2))

d1_stemp_dy_0_nxp4m2nyp4p0m2k = q(nx+4-2,ny+4+0-2,indvars(2))

d1_stemp_dy_0_nxp4m2nyp4p0k = 1.5_wp*d1_stemp_dy_0_nxp4m2nyp4p0p0k-&
          2.0_wp*d1_stemp_dy_0_nxp4m2nyp4p0m1k+&
          0.5_wp*d1_stemp_dy_0_nxp4m2nyp4p0m2k

d1_stemp_dy_0_nxp4m2nyp4p0k = d1_stemp_dy_0_nxp4m2nyp4p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 2 0 None stemp *****************
!                                                           
!***********************************************************


qst(nx+4-2,ny+4+0,indvarsst(4)) =  (((0.5_wp*(qst(nx+4-2,ny+4+0,indvarsst(11))*(d1_stemp_dy_0_nxp4m2nyp4p0k)-&
                    qst(nx+4-2,ny+4+0,indvarsst(10))*(d1_stemp_dx_0_nxp4m2nyp4p0k)))**2+&
                    (0.5_wp*(qst(nx+4-2,ny+4+0,indvarsst(10))*(d1_stemp_dx_0_nxp4m2nyp4p0k)-&
                    qst(nx+4-2,ny+4+0,indvarsst(11))*(d1_stemp_dy_0_nxp4m2nyp4p0k)))**2)*2)**0.5



!***********************************************************
!                                                           
! building source terms in RHS for layer 2 0 None SS *******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (stemp+ReI*nut/(k**2*eta**2))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 2 0 None SS ********************
!                                                           
!***********************************************************


qst(nx+4-2,ny+4+0,indvarsst(12)) =  (qst(nx+4-2,ny+4+0,indvarsst(4))+&
                    param_float(1 + 5)*q(nx+4-2,ny+4+0,indvars(5))/(param_float(9 + 5)**2*qst(nx+4-2,ny+4+0,indvarsst(2))**2))



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


qst(nx+4-2,ny+4+0,indvarsst(13)) =  (((1+&
                    param_float(21 + 5))/(((q(nx+4-2,ny+4+0,indvars(4))-&
                    0.5_wp*(q(nx+4-2,ny+4+0,indvars(2))*q(nx+4-2,ny+4+0,indvars(2))+&
                    q(nx+4-2,ny+4+0,indvars(3))*q(nx+4-2,ny+4+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4-2,ny+4+0,indvars(4))-&
                    0.5_wp*(q(nx+4-2,ny+4+0,indvars(2))*q(nx+4-2,ny+4+0,indvars(2))+&
                    q(nx+4-2,ny+4+0,indvars(3))*q(nx+4-2,ny+4+0,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+4-2,ny+4+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4-2,ny+4+0,indvars(4))-&
                    0.5_wp*(q(nx+4-2,ny+4+0,indvars(2))*q(nx+4-2,ny+4+0,indvars(2))+&
                    q(nx+4-2,ny+4+0,indvars(3))*q(nx+4-2,ny+4+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4-2,ny+4+0,indvars(4))-&
                    0.5_wp*(q(nx+4-2,ny+4+0,indvars(2))*q(nx+4-2,ny+4+0,indvars(2))+&
                    q(nx+4-2,ny+4+0,indvars(3))*q(nx+4-2,ny+4+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4-2,ny+4+0,indvars(1)))**3/((q(nx+4-2,ny+4+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4-2,ny+4+0,indvars(4))-&
                    0.5_wp*(q(nx+4-2,ny+4+0,indvars(2))*q(nx+4-2,ny+4+0,indvars(2))+&
                    q(nx+4-2,ny+4+0,indvars(3))*q(nx+4-2,ny+4+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4-2,ny+4+0,indvars(4))-&
                    0.5_wp*(q(nx+4-2,ny+4+0,indvars(2))*q(nx+4-2,ny+4+0,indvars(2))+&
                    q(nx+4-2,ny+4+0,indvars(3))*q(nx+4-2,ny+4+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4-2,ny+4+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+4-2,ny+4+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4-2,ny+4+0,indvars(4))-&
                    0.5_wp*(q(nx+4-2,ny+4+0,indvars(2))*q(nx+4-2,ny+4+0,indvars(2))+&
                    q(nx+4-2,ny+4+0,indvars(3))*q(nx+4-2,ny+4+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4-2,ny+4+0,indvars(4))-&
                    0.5_wp*(q(nx+4-2,ny+4+0,indvars(2))*q(nx+4-2,ny+4+0,indvars(2))+&
                    q(nx+4-2,ny+4+0,indvars(3))*q(nx+4-2,ny+4+0,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4-2,ny+4+0,indvars(1))))*param_float(1 + 5)*(d1_stemp_dy_0_nxp4m2nyp4p0k)*qst(nx+4-2,ny+4+0,indvarsst(11)))

