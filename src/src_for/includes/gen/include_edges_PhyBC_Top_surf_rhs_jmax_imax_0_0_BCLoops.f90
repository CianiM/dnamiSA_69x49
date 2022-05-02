

!***********************************************************
!                                                           
! Start building layers for BC : imax jmax None ************
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
! deltaxI*([rho*u]_1x)+deltayI*([rho*v]_1y)+(1.0_wp/c**2*((v*(-c**2*[rho]_1y+1.0_wp*[p]_1y))*deltayI+0.5_wp*(((c+v)*(c*rho*[v]_1y+1.0_wp*[p]_1y))*deltayI+esse*c*(1-M_jmax*M_jmax)/L_ref*(p-P0))))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_rhs_rho_dx_0_nxp2p0p0nyp2p0k = q(nx+2+0+0,ny+2+0,indvars(1))*q(nx+2+0+0,ny+2+0,indvars(2))

d1_rhs_rho_dx_0_nxp2p0m1nyp2p0k = q(nx+2+0-1,ny+2+0,indvars(1))*q(nx+2+0-1,ny+2+0,indvars(2))

d1_rhs_rho_dx_0_nxp2p0m2nyp2p0k = q(nx+2+0-2,ny+2+0,indvars(1))*q(nx+2+0-2,ny+2+0,indvars(2))

d1_rhs_rho_dx_0_nxp2p0nyp2p0k = 1.5_wp*d1_rhs_rho_dx_0_nxp2p0p0nyp2p0k-&
          2.0_wp*d1_rhs_rho_dx_0_nxp2p0m1nyp2p0k+&
          0.5_wp*d1_rhs_rho_dx_0_nxp2p0m2nyp2p0k

d1_rhs_rho_dx_0_nxp2p0nyp2p0k = d1_rhs_rho_dx_0_nxp2p0nyp2p0k*param_float(1)

d1_rhs_rho_dy_0_nxp2p0nyp2p0p0k = q(nx+2+0,ny+2+0+0,indvars(1))*q(nx+2+0,ny+2+0+0,indvars(3))

d1_rhs_rho_dy_0_nxp2p0nyp2p0m1k = q(nx+2+0,ny+2+0-1,indvars(1))*q(nx+2+0,ny+2+0-1,indvars(3))

d1_rhs_rho_dy_0_nxp2p0nyp2p0m2k = q(nx+2+0,ny+2+0-2,indvars(1))*q(nx+2+0,ny+2+0-2,indvars(3))

d1_rhs_rho_dy_0_nxp2p0nyp2p0k = 1.5_wp*d1_rhs_rho_dy_0_nxp2p0nyp2p0p0k-&
          2.0_wp*d1_rhs_rho_dy_0_nxp2p0nyp2p0m1k+&
          0.5_wp*d1_rhs_rho_dy_0_nxp2p0nyp2p0m2k

d1_rhs_rho_dy_0_nxp2p0nyp2p0k = d1_rhs_rho_dy_0_nxp2p0nyp2p0k*param_float(2)

d1_rhs_rho_dy_1_nxp2p0nyp2p0p0k = q(nx+2+0,ny+2+0+0,indvars(1))

d1_rhs_rho_dy_1_nxp2p0nyp2p0m1k = q(nx+2+0,ny+2+0-1,indvars(1))

d1_rhs_rho_dy_1_nxp2p0nyp2p0m2k = q(nx+2+0,ny+2+0-2,indvars(1))

d1_rhs_rho_dy_1_nxp2p0nyp2p0k = 1.5_wp*d1_rhs_rho_dy_1_nxp2p0nyp2p0p0k-&
          2.0_wp*d1_rhs_rho_dy_1_nxp2p0nyp2p0m1k+&
          0.5_wp*d1_rhs_rho_dy_1_nxp2p0nyp2p0m2k

d1_rhs_rho_dy_1_nxp2p0nyp2p0k = d1_rhs_rho_dy_1_nxp2p0nyp2p0k*param_float(2)

d1_rhs_rho_dy_2_nxp2p0nyp2p0p0k = (param_float(3 + 5))*q(nx+2+0,ny+2+0+0,indvars(1))*((q(nx+2+0,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0+0,indvars(2))*q(nx+2+0,ny+2+0+0,indvars(2))+&
                    q(nx+2+0,ny+2+0+0,indvars(3))*q(nx+2+0,ny+2+0+0,indvars(3)))))

d1_rhs_rho_dy_2_nxp2p0nyp2p0m1k = (param_float(3 + 5))*q(nx+2+0,ny+2+0-1,indvars(1))*((q(nx+2+0,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0-1,indvars(2))*q(nx+2+0,ny+2+0-1,indvars(2))+&
                    q(nx+2+0,ny+2+0-1,indvars(3))*q(nx+2+0,ny+2+0-1,indvars(3)))))

d1_rhs_rho_dy_2_nxp2p0nyp2p0m2k = (param_float(3 + 5))*q(nx+2+0,ny+2+0-2,indvars(1))*((q(nx+2+0,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0-2,indvars(2))*q(nx+2+0,ny+2+0-2,indvars(2))+&
                    q(nx+2+0,ny+2+0-2,indvars(3))*q(nx+2+0,ny+2+0-2,indvars(3)))))

d1_rhs_rho_dy_2_nxp2p0nyp2p0k = 1.5_wp*d1_rhs_rho_dy_2_nxp2p0nyp2p0p0k-&
          2.0_wp*d1_rhs_rho_dy_2_nxp2p0nyp2p0m1k+&
          0.5_wp*d1_rhs_rho_dy_2_nxp2p0nyp2p0m2k

d1_rhs_rho_dy_2_nxp2p0nyp2p0k = d1_rhs_rho_dy_2_nxp2p0nyp2p0k*param_float(2)

d1_rhs_rho_dy_3_nxp2p0nyp2p0p0k = q(nx+2+0,ny+2+0+0,indvars(3))

d1_rhs_rho_dy_3_nxp2p0nyp2p0m1k = q(nx+2+0,ny+2+0-1,indvars(3))

d1_rhs_rho_dy_3_nxp2p0nyp2p0m2k = q(nx+2+0,ny+2+0-2,indvars(3))

d1_rhs_rho_dy_3_nxp2p0nyp2p0k = 1.5_wp*d1_rhs_rho_dy_3_nxp2p0nyp2p0p0k-&
          2.0_wp*d1_rhs_rho_dy_3_nxp2p0nyp2p0m1k+&
          0.5_wp*d1_rhs_rho_dy_3_nxp2p0nyp2p0m2k

d1_rhs_rho_dy_3_nxp2p0nyp2p0k = d1_rhs_rho_dy_3_nxp2p0nyp2p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None d(rho)/dt *************
!                                                           
!***********************************************************


rhs(nx+2+0,ny+2+0,indvars(1)) =   -  ( qst(nx+2+0,ny+2+0,indvarsst(10))*(d1_rhs_rho_dx_0_nxp2p0nyp2p0k)+&
                    qst(nx+2+0,ny+2+0,indvarsst(11))*(d1_rhs_rho_dy_0_nxp2p0nyp2p0k)+&
                    (1.0_wp/(param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5**2*((q(nx+2+0,ny+2+0,indvars(3))*(-&
                    (param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5**2*d1_rhs_rho_dy_1_nxp2p0nyp2p0k+&
                    1.0_wp*d1_rhs_rho_dy_2_nxp2p0nyp2p0k))*qst(nx+2+0,ny+2+0,indvarsst(11))+&
                    0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5+&
                    q(nx+2+0,ny+2+0,indvars(3)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5*q(nx+2+0,ny+2+0,indvars(1))*d1_rhs_rho_dy_3_nxp2p0nyp2p0k+&
                    1.0_wp*d1_rhs_rho_dy_2_nxp2p0nyp2p0k))*qst(nx+2+0,ny+2+0,indvarsst(11))+&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5*(1-&
                    qface_j(nx+2+0,1)*qface_j(nx+2+0,1))/param_float(20 + 5)*((param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))-&
                    param_float(26 + 5))))) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None d(rho u)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! u*(1.0_wp/c**2*((v*(-c**2*[rho]_1y+1.0_wp*[p]_1y))*deltayI+0.5_wp*(((c+v)*(c*rho*[v]_1y+1.0_wp*[p]_1y))*deltayI+esse*c*(1-M_jmax*M_jmax)/L_ref*(p-P0))))+rho*(1.0_wp/(2*rho*c)*(((c+v)*(c*rho*[v]_1y+1.0_wp*[p]_1y))*deltayI-esse*c*(1-M_jmax*M_jmax)/L_ref*(p-P0)))+[rho*u*v]_1y+deltaxI*([-visc_t*(2.0_wp*deltaxI*({u}_1x)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y)))]_1x)+deltayI*([-visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_rhs_rhou_dxdx_0_0_nxp2p0p0nyp2p0k_nxp2p0p0p0nyp2p0k = q(nx+2+0+0+0,ny+2+0,indvars(2))

d2_rhs_rhou_dxdx_0_0_nxp2p0p0nyp2p0k_nxp2p0p0m1nyp2p0k = q(nx+2+0+0-1,ny+2+0,indvars(2))

d2_rhs_rhou_dxdx_0_0_nxp2p0p0nyp2p0k_nxp2p0p0m2nyp2p0k = q(nx+2+0+0-2,ny+2+0,indvars(2))

d2_rhs_rhou_dxdx_0_0_nxp2p0p0nyp2p0k = 1.5_wp*d2_rhs_rhou_dxdx_0_0_nxp2p0p0nyp2p0k_nxp2p0p0p0nyp2p0k-&
          2.0_wp*d2_rhs_rhou_dxdx_0_0_nxp2p0p0nyp2p0k_nxp2p0p0m1nyp2p0k+&
          0.5_wp*d2_rhs_rhou_dxdx_0_0_nxp2p0p0nyp2p0k_nxp2p0p0m2nyp2p0k

d2_rhs_rhou_dxdx_0_0_nxp2p0p0nyp2p0k = d2_rhs_rhou_dxdx_0_0_nxp2p0p0nyp2p0k*param_float(1)

d2_rhs_rhou_dxdx_0_0_nxp2p0m1nyp2p0k_nxp2p0m1m1nyp2p0k = q(nx+2+0-1-1,ny+2+0,indvars(2))

d2_rhs_rhou_dxdx_0_0_nxp2p0m1nyp2p0k_nxp2p0m1p1nyp2p0k = q(nx+2+0-1+1,ny+2+0,indvars(2))

d2_rhs_rhou_dxdx_0_0_nxp2p0m1nyp2p0k = -&
          0.5_wp*d2_rhs_rhou_dxdx_0_0_nxp2p0m1nyp2p0k_nxp2p0m1m1nyp2p0k+&
          0.5_wp*d2_rhs_rhou_dxdx_0_0_nxp2p0m1nyp2p0k_nxp2p0m1p1nyp2p0k

d2_rhs_rhou_dxdx_0_0_nxp2p0m1nyp2p0k = d2_rhs_rhou_dxdx_0_0_nxp2p0m1nyp2p0k*param_float(1)

d2_rhs_rhou_dxdx_0_0_nxp2p0m2nyp2p0k_nxp2p0m2m1nyp2p0k = q(nx+2+0-2-1,ny+2+0,indvars(2))

d2_rhs_rhou_dxdx_0_0_nxp2p0m2nyp2p0k_nxp2p0m2p1nyp2p0k = q(nx+2+0-2+1,ny+2+0,indvars(2))

d2_rhs_rhou_dxdx_0_0_nxp2p0m2nyp2p0k = -&
          0.5_wp*d2_rhs_rhou_dxdx_0_0_nxp2p0m2nyp2p0k_nxp2p0m2m1nyp2p0k+&
          0.5_wp*d2_rhs_rhou_dxdx_0_0_nxp2p0m2nyp2p0k_nxp2p0m2p1nyp2p0k

d2_rhs_rhou_dxdx_0_0_nxp2p0m2nyp2p0k = d2_rhs_rhou_dxdx_0_0_nxp2p0m2nyp2p0k*param_float(1)

d2_rhs_rhou_dxdy_0_0_nxp2p0p0nyp2p0k_nxp2p0p0nyp2p0p0k = q(nx+2+0+0,ny+2+0+0,indvars(3))

d2_rhs_rhou_dxdy_0_0_nxp2p0p0nyp2p0k_nxp2p0p0nyp2p0m1k = q(nx+2+0+0,ny+2+0-1,indvars(3))

d2_rhs_rhou_dxdy_0_0_nxp2p0p0nyp2p0k_nxp2p0p0nyp2p0m2k = q(nx+2+0+0,ny+2+0-2,indvars(3))

d2_rhs_rhou_dxdy_0_0_nxp2p0p0nyp2p0k = 1.5_wp*d2_rhs_rhou_dxdy_0_0_nxp2p0p0nyp2p0k_nxp2p0p0nyp2p0p0k-&
          2.0_wp*d2_rhs_rhou_dxdy_0_0_nxp2p0p0nyp2p0k_nxp2p0p0nyp2p0m1k+&
          0.5_wp*d2_rhs_rhou_dxdy_0_0_nxp2p0p0nyp2p0k_nxp2p0p0nyp2p0m2k

d2_rhs_rhou_dxdy_0_0_nxp2p0p0nyp2p0k = d2_rhs_rhou_dxdy_0_0_nxp2p0p0nyp2p0k*param_float(2)

d2_rhs_rhou_dxdy_0_0_nxp2p0m1nyp2p0k_nxp2p0m1nyp2p0p0k = q(nx+2+0-1,ny+2+0+0,indvars(3))

d2_rhs_rhou_dxdy_0_0_nxp2p0m1nyp2p0k_nxp2p0m1nyp2p0m1k = q(nx+2+0-1,ny+2+0-1,indvars(3))

d2_rhs_rhou_dxdy_0_0_nxp2p0m1nyp2p0k_nxp2p0m1nyp2p0m2k = q(nx+2+0-1,ny+2+0-2,indvars(3))

d2_rhs_rhou_dxdy_0_0_nxp2p0m1nyp2p0k = 1.5_wp*d2_rhs_rhou_dxdy_0_0_nxp2p0m1nyp2p0k_nxp2p0m1nyp2p0p0k-&
          2.0_wp*d2_rhs_rhou_dxdy_0_0_nxp2p0m1nyp2p0k_nxp2p0m1nyp2p0m1k+&
          0.5_wp*d2_rhs_rhou_dxdy_0_0_nxp2p0m1nyp2p0k_nxp2p0m1nyp2p0m2k

d2_rhs_rhou_dxdy_0_0_nxp2p0m1nyp2p0k = d2_rhs_rhou_dxdy_0_0_nxp2p0m1nyp2p0k*param_float(2)

d2_rhs_rhou_dxdy_0_0_nxp2p0m2nyp2p0k_nxp2p0m2nyp2p0p0k = q(nx+2+0-2,ny+2+0+0,indvars(3))

d2_rhs_rhou_dxdy_0_0_nxp2p0m2nyp2p0k_nxp2p0m2nyp2p0m1k = q(nx+2+0-2,ny+2+0-1,indvars(3))

d2_rhs_rhou_dxdy_0_0_nxp2p0m2nyp2p0k_nxp2p0m2nyp2p0m2k = q(nx+2+0-2,ny+2+0-2,indvars(3))

d2_rhs_rhou_dxdy_0_0_nxp2p0m2nyp2p0k = 1.5_wp*d2_rhs_rhou_dxdy_0_0_nxp2p0m2nyp2p0k_nxp2p0m2nyp2p0p0k-&
          2.0_wp*d2_rhs_rhou_dxdy_0_0_nxp2p0m2nyp2p0k_nxp2p0m2nyp2p0m1k+&
          0.5_wp*d2_rhs_rhou_dxdy_0_0_nxp2p0m2nyp2p0k_nxp2p0m2nyp2p0m2k

d2_rhs_rhou_dxdy_0_0_nxp2p0m2nyp2p0k = d2_rhs_rhou_dxdy_0_0_nxp2p0m2nyp2p0k*param_float(2)

d1_rhs_rhou_dx_0_nxp2p0p0nyp2p0k = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+2+0+0,ny+2+0,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2+0+0,ny+2+0,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+2+0+0,ny+2+0,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(nx+2+0+0,ny+2+0,indvarsst(10))*(d2_rhs_rhou_dxdx_0_0_nxp2p0p0nyp2p0k)-&
                    2.0_wp/3.0_wp*(qst(nx+2+0+0,ny+2+0,indvarsst(10))*(d2_rhs_rhou_dxdx_0_0_nxp2p0p0nyp2p0k)+&
                    qst(nx+2+0+0,ny+2+0,indvarsst(11))*(d2_rhs_rhou_dxdy_0_0_nxp2p0p0nyp2p0k)))

d1_rhs_rhou_dx_0_nxp2p0m1nyp2p0k = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+2+0-1,ny+2+0,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2+0-1,ny+2+0,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+2+0-1,ny+2+0,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(nx+2+0-1,ny+2+0,indvarsst(10))*(d2_rhs_rhou_dxdx_0_0_nxp2p0m1nyp2p0k)-&
                    2.0_wp/3.0_wp*(qst(nx+2+0-1,ny+2+0,indvarsst(10))*(d2_rhs_rhou_dxdx_0_0_nxp2p0m1nyp2p0k)+&
                    qst(nx+2+0-1,ny+2+0,indvarsst(11))*(d2_rhs_rhou_dxdy_0_0_nxp2p0m1nyp2p0k)))

d1_rhs_rhou_dx_0_nxp2p0m2nyp2p0k = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+2+0-2,ny+2+0,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2+0-2,ny+2+0,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+2+0-2,ny+2+0,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(nx+2+0-2,ny+2+0,indvarsst(10))*(d2_rhs_rhou_dxdx_0_0_nxp2p0m2nyp2p0k)-&
                    2.0_wp/3.0_wp*(qst(nx+2+0-2,ny+2+0,indvarsst(10))*(d2_rhs_rhou_dxdx_0_0_nxp2p0m2nyp2p0k)+&
                    qst(nx+2+0-2,ny+2+0,indvarsst(11))*(d2_rhs_rhou_dxdy_0_0_nxp2p0m2nyp2p0k)))

d1_rhs_rhou_dx_0_nxp2p0nyp2p0k = 1.5_wp*d1_rhs_rhou_dx_0_nxp2p0p0nyp2p0k-&
          2.0_wp*d1_rhs_rhou_dx_0_nxp2p0m1nyp2p0k+&
          0.5_wp*d1_rhs_rhou_dx_0_nxp2p0m2nyp2p0k

d1_rhs_rhou_dx_0_nxp2p0nyp2p0k = d1_rhs_rhou_dx_0_nxp2p0nyp2p0k*param_float(1)

d1_rhs_rhou_dy_6_nxp2p0nyp2p0p0k = q(nx+2+0,ny+2+0+0,indvars(1))*q(nx+2+0,ny+2+0+0,indvars(2))*q(nx+2+0,ny+2+0+0,indvars(3))

d1_rhs_rhou_dy_6_nxp2p0nyp2p0m1k = q(nx+2+0,ny+2+0-1,indvars(1))*q(nx+2+0,ny+2+0-1,indvars(2))*q(nx+2+0,ny+2+0-1,indvars(3))

d1_rhs_rhou_dy_6_nxp2p0nyp2p0m2k = q(nx+2+0,ny+2+0-2,indvars(1))*q(nx+2+0,ny+2+0-2,indvars(2))*q(nx+2+0,ny+2+0-2,indvars(3))

d1_rhs_rhou_dy_6_nxp2p0nyp2p0k = 1.5_wp*d1_rhs_rhou_dy_6_nxp2p0nyp2p0p0k-&
          2.0_wp*d1_rhs_rhou_dy_6_nxp2p0nyp2p0m1k+&
          0.5_wp*d1_rhs_rhou_dy_6_nxp2p0nyp2p0m2k

d1_rhs_rhou_dy_6_nxp2p0nyp2p0k = d1_rhs_rhou_dy_6_nxp2p0nyp2p0k*param_float(2)

d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0p0k_nxp2p0p0nyp2p0p0k = q(nx+2+0+0,ny+2+0+0,indvars(3))

d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0p0k_nxp2p0m1nyp2p0p0k = q(nx+2+0-1,ny+2+0+0,indvars(3))

d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0p0k_nxp2p0m2nyp2p0p0k = q(nx+2+0-2,ny+2+0+0,indvars(3))

d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0p0k = 1.5_wp*d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0p0k_nxp2p0p0nyp2p0p0k-&
          2.0_wp*d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0p0k_nxp2p0m1nyp2p0p0k+&
          0.5_wp*d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0p0k_nxp2p0m2nyp2p0p0k

d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0p0k = d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0p0k*param_float(1)

d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0m1k_nxp2p0p0nyp2p0m1k = q(nx+2+0+0,ny+2+0-1,indvars(3))

d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0m1k_nxp2p0m1nyp2p0m1k = q(nx+2+0-1,ny+2+0-1,indvars(3))

d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0m1k_nxp2p0m2nyp2p0m1k = q(nx+2+0-2,ny+2+0-1,indvars(3))

d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0m1k = 1.5_wp*d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0m1k_nxp2p0p0nyp2p0m1k-&
          2.0_wp*d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0m1k_nxp2p0m1nyp2p0m1k+&
          0.5_wp*d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0m1k_nxp2p0m2nyp2p0m1k

d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0m1k = d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0m1k*param_float(1)

d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0m2k_nxp2p0p0nyp2p0m2k = q(nx+2+0+0,ny+2+0-2,indvars(3))

d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0m2k_nxp2p0m1nyp2p0m2k = q(nx+2+0-1,ny+2+0-2,indvars(3))

d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0m2k_nxp2p0m2nyp2p0m2k = q(nx+2+0-2,ny+2+0-2,indvars(3))

d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0m2k = 1.5_wp*d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0m2k_nxp2p0p0nyp2p0m2k-&
          2.0_wp*d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0m2k_nxp2p0m1nyp2p0m2k+&
          0.5_wp*d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0m2k_nxp2p0m2nyp2p0m2k

d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0m2k = d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0m2k*param_float(1)

d2_rhs_rhou_dydy_7_0_nxp2p0nyp2p0p0k_nxp2p0nyp2p0p0p0k = q(nx+2+0,ny+2+0+0+0,indvars(2))

d2_rhs_rhou_dydy_7_0_nxp2p0nyp2p0p0k_nxp2p0nyp2p0p0m1k = q(nx+2+0,ny+2+0+0-1,indvars(2))

d2_rhs_rhou_dydy_7_0_nxp2p0nyp2p0p0k_nxp2p0nyp2p0p0m2k = q(nx+2+0,ny+2+0+0-2,indvars(2))

d2_rhs_rhou_dydy_7_0_nxp2p0nyp2p0p0k = 1.5_wp*d2_rhs_rhou_dydy_7_0_nxp2p0nyp2p0p0k_nxp2p0nyp2p0p0p0k-&
          2.0_wp*d2_rhs_rhou_dydy_7_0_nxp2p0nyp2p0p0k_nxp2p0nyp2p0p0m1k+&
          0.5_wp*d2_rhs_rhou_dydy_7_0_nxp2p0nyp2p0p0k_nxp2p0nyp2p0p0m2k

d2_rhs_rhou_dydy_7_0_nxp2p0nyp2p0p0k = d2_rhs_rhou_dydy_7_0_nxp2p0nyp2p0p0k*param_float(2)

d2_rhs_rhou_dydy_7_0_nxp2p0nyp2p0m1k_nxp2p0nyp2p0m1m1k = q(nx+2+0,ny+2+0-1-1,indvars(2))

d2_rhs_rhou_dydy_7_0_nxp2p0nyp2p0m1k_nxp2p0nyp2p0m1p1k = q(nx+2+0,ny+2+0-1+1,indvars(2))

d2_rhs_rhou_dydy_7_0_nxp2p0nyp2p0m1k = -&
          0.5_wp*d2_rhs_rhou_dydy_7_0_nxp2p0nyp2p0m1k_nxp2p0nyp2p0m1m1k+&
          0.5_wp*d2_rhs_rhou_dydy_7_0_nxp2p0nyp2p0m1k_nxp2p0nyp2p0m1p1k

d2_rhs_rhou_dydy_7_0_nxp2p0nyp2p0m1k = d2_rhs_rhou_dydy_7_0_nxp2p0nyp2p0m1k*param_float(2)

d2_rhs_rhou_dydy_7_0_nxp2p0nyp2p0m2k_nxp2p0nyp2p0m2m1k = q(nx+2+0,ny+2+0-2-1,indvars(2))

d2_rhs_rhou_dydy_7_0_nxp2p0nyp2p0m2k_nxp2p0nyp2p0m2p1k = q(nx+2+0,ny+2+0-2+1,indvars(2))

d2_rhs_rhou_dydy_7_0_nxp2p0nyp2p0m2k = -&
          0.5_wp*d2_rhs_rhou_dydy_7_0_nxp2p0nyp2p0m2k_nxp2p0nyp2p0m2m1k+&
          0.5_wp*d2_rhs_rhou_dydy_7_0_nxp2p0nyp2p0m2k_nxp2p0nyp2p0m2p1k

d2_rhs_rhou_dydy_7_0_nxp2p0nyp2p0m2k = d2_rhs_rhou_dydy_7_0_nxp2p0nyp2p0m2k*param_float(2)

d1_rhs_rhou_dy_7_nxp2p0nyp2p0p0k = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+2+0,ny+2+0+0,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2+0,ny+2+0+0,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+2+0,ny+2+0+0,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(nx+2+0,ny+2+0+0,indvarsst(11))*(d2_rhs_rhou_dydy_7_0_nxp2p0nyp2p0p0k)+&
                    qst(nx+2+0,ny+2+0+0,indvarsst(10))*(d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0p0k))

d1_rhs_rhou_dy_7_nxp2p0nyp2p0m1k = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+2+0,ny+2+0-1,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2+0,ny+2+0-1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+2+0,ny+2+0-1,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(nx+2+0,ny+2+0-1,indvarsst(11))*(d2_rhs_rhou_dydy_7_0_nxp2p0nyp2p0m1k)+&
                    qst(nx+2+0,ny+2+0-1,indvarsst(10))*(d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0m1k))

d1_rhs_rhou_dy_7_nxp2p0nyp2p0m2k = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+2+0,ny+2+0-2,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2+0,ny+2+0-2,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+2+0,ny+2+0-2,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(nx+2+0,ny+2+0-2,indvarsst(11))*(d2_rhs_rhou_dydy_7_0_nxp2p0nyp2p0m2k)+&
                    qst(nx+2+0,ny+2+0-2,indvarsst(10))*(d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0m2k))

d1_rhs_rhou_dy_7_nxp2p0nyp2p0k = 1.5_wp*d1_rhs_rhou_dy_7_nxp2p0nyp2p0p0k-&
          2.0_wp*d1_rhs_rhou_dy_7_nxp2p0nyp2p0m1k+&
          0.5_wp*d1_rhs_rhou_dy_7_nxp2p0nyp2p0m2k

d1_rhs_rhou_dy_7_nxp2p0nyp2p0k = d1_rhs_rhou_dy_7_nxp2p0nyp2p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None d(rho u)/dt ***********
!                                                           
!***********************************************************


rhs(nx+2+0,ny+2+0,indvars(2)) =   -  ( q(nx+2+0,ny+2+0,indvars(2))*(1.0_wp/(param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5**2*((q(nx+2+0,ny+2+0,indvars(3))*(-&
                    (param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5**2*d1_rhs_rho_dy_1_nxp2p0nyp2p0k+&
                    1.0_wp*d1_rhs_rho_dy_2_nxp2p0nyp2p0k))*qst(nx+2+0,ny+2+0,indvarsst(11))+&
                    0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5+&
                    q(nx+2+0,ny+2+0,indvars(3)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5*q(nx+2+0,ny+2+0,indvars(1))*d1_rhs_rho_dy_3_nxp2p0nyp2p0k+&
                    1.0_wp*d1_rhs_rho_dy_2_nxp2p0nyp2p0k))*qst(nx+2+0,ny+2+0,indvarsst(11))+&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5*(1-&
                    qface_j(nx+2+0,1)*qface_j(nx+2+0,1))/param_float(20 + 5)*((param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))-&
                    param_float(26 + 5)))))+&
                    q(nx+2+0,ny+2+0,indvars(1))*(1.0_wp/(2*q(nx+2+0,ny+2+0,indvars(1))*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5)*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5+&
                    q(nx+2+0,ny+2+0,indvars(3)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5*q(nx+2+0,ny+2+0,indvars(1))*d1_rhs_rho_dy_3_nxp2p0nyp2p0k+&
                    1.0_wp*d1_rhs_rho_dy_2_nxp2p0nyp2p0k))*qst(nx+2+0,ny+2+0,indvarsst(11))-&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5*(1-&
                    qface_j(nx+2+0,1)*qface_j(nx+2+0,1))/param_float(20 + 5)*((param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))-&
                    param_float(26 + 5))))+&
                    d1_rhs_rhou_dy_6_nxp2p0nyp2p0k+&
                    qst(nx+2+0,ny+2+0,indvarsst(10))*(d1_rhs_rhou_dx_0_nxp2p0nyp2p0k)+&
                    qst(nx+2+0,ny+2+0,indvarsst(11))*(d1_rhs_rhou_dy_7_nxp2p0nyp2p0k) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None d(rho v)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! v*(1.0_wp/c**2*((v*(-c**2*[rho]_1y+1.0_wp*[p]_1y))*deltayI+0.5_wp*(((c+v)*(c*rho*[v]_1y+1.0_wp*[p]_1y))*deltayI+esse*c*(1-M_jmax*M_jmax)/L_ref*(p-P0))))+rho*((1.0_wp*v*[u]_1y)*deltayI)+[rho*v*v]_1y+deltaxI*([-visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x))]_1x)+deltayI*([-visc_t*(2.0_wp*deltayI*({v}_1y)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y)))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_rhs_rhov_dxdx_0_0_nxp2p0p0nyp2p0k_nxp2p0p0p0nyp2p0k = q(nx+2+0+0+0,ny+2+0,indvars(3))

d2_rhs_rhov_dxdx_0_0_nxp2p0p0nyp2p0k_nxp2p0p0m1nyp2p0k = q(nx+2+0+0-1,ny+2+0,indvars(3))

d2_rhs_rhov_dxdx_0_0_nxp2p0p0nyp2p0k_nxp2p0p0m2nyp2p0k = q(nx+2+0+0-2,ny+2+0,indvars(3))

d2_rhs_rhov_dxdx_0_0_nxp2p0p0nyp2p0k = 1.5_wp*d2_rhs_rhov_dxdx_0_0_nxp2p0p0nyp2p0k_nxp2p0p0p0nyp2p0k-&
          2.0_wp*d2_rhs_rhov_dxdx_0_0_nxp2p0p0nyp2p0k_nxp2p0p0m1nyp2p0k+&
          0.5_wp*d2_rhs_rhov_dxdx_0_0_nxp2p0p0nyp2p0k_nxp2p0p0m2nyp2p0k

d2_rhs_rhov_dxdx_0_0_nxp2p0p0nyp2p0k = d2_rhs_rhov_dxdx_0_0_nxp2p0p0nyp2p0k*param_float(1)

d2_rhs_rhov_dxdx_0_0_nxp2p0m1nyp2p0k_nxp2p0m1m1nyp2p0k = q(nx+2+0-1-1,ny+2+0,indvars(3))

d2_rhs_rhov_dxdx_0_0_nxp2p0m1nyp2p0k_nxp2p0m1p1nyp2p0k = q(nx+2+0-1+1,ny+2+0,indvars(3))

d2_rhs_rhov_dxdx_0_0_nxp2p0m1nyp2p0k = -&
          0.5_wp*d2_rhs_rhov_dxdx_0_0_nxp2p0m1nyp2p0k_nxp2p0m1m1nyp2p0k+&
          0.5_wp*d2_rhs_rhov_dxdx_0_0_nxp2p0m1nyp2p0k_nxp2p0m1p1nyp2p0k

d2_rhs_rhov_dxdx_0_0_nxp2p0m1nyp2p0k = d2_rhs_rhov_dxdx_0_0_nxp2p0m1nyp2p0k*param_float(1)

d2_rhs_rhov_dxdx_0_0_nxp2p0m2nyp2p0k_nxp2p0m2m1nyp2p0k = q(nx+2+0-2-1,ny+2+0,indvars(3))

d2_rhs_rhov_dxdx_0_0_nxp2p0m2nyp2p0k_nxp2p0m2p1nyp2p0k = q(nx+2+0-2+1,ny+2+0,indvars(3))

d2_rhs_rhov_dxdx_0_0_nxp2p0m2nyp2p0k = -&
          0.5_wp*d2_rhs_rhov_dxdx_0_0_nxp2p0m2nyp2p0k_nxp2p0m2m1nyp2p0k+&
          0.5_wp*d2_rhs_rhov_dxdx_0_0_nxp2p0m2nyp2p0k_nxp2p0m2p1nyp2p0k

d2_rhs_rhov_dxdx_0_0_nxp2p0m2nyp2p0k = d2_rhs_rhov_dxdx_0_0_nxp2p0m2nyp2p0k*param_float(1)

d2_rhs_rhov_dxdy_0_0_nxp2p0p0nyp2p0k_nxp2p0p0nyp2p0p0k = q(nx+2+0+0,ny+2+0+0,indvars(2))

d2_rhs_rhov_dxdy_0_0_nxp2p0p0nyp2p0k_nxp2p0p0nyp2p0m1k = q(nx+2+0+0,ny+2+0-1,indvars(2))

d2_rhs_rhov_dxdy_0_0_nxp2p0p0nyp2p0k_nxp2p0p0nyp2p0m2k = q(nx+2+0+0,ny+2+0-2,indvars(2))

d2_rhs_rhov_dxdy_0_0_nxp2p0p0nyp2p0k = 1.5_wp*d2_rhs_rhov_dxdy_0_0_nxp2p0p0nyp2p0k_nxp2p0p0nyp2p0p0k-&
          2.0_wp*d2_rhs_rhov_dxdy_0_0_nxp2p0p0nyp2p0k_nxp2p0p0nyp2p0m1k+&
          0.5_wp*d2_rhs_rhov_dxdy_0_0_nxp2p0p0nyp2p0k_nxp2p0p0nyp2p0m2k

d2_rhs_rhov_dxdy_0_0_nxp2p0p0nyp2p0k = d2_rhs_rhov_dxdy_0_0_nxp2p0p0nyp2p0k*param_float(2)

d2_rhs_rhov_dxdy_0_0_nxp2p0m1nyp2p0k_nxp2p0m1nyp2p0p0k = q(nx+2+0-1,ny+2+0+0,indvars(2))

d2_rhs_rhov_dxdy_0_0_nxp2p0m1nyp2p0k_nxp2p0m1nyp2p0m1k = q(nx+2+0-1,ny+2+0-1,indvars(2))

d2_rhs_rhov_dxdy_0_0_nxp2p0m1nyp2p0k_nxp2p0m1nyp2p0m2k = q(nx+2+0-1,ny+2+0-2,indvars(2))

d2_rhs_rhov_dxdy_0_0_nxp2p0m1nyp2p0k = 1.5_wp*d2_rhs_rhov_dxdy_0_0_nxp2p0m1nyp2p0k_nxp2p0m1nyp2p0p0k-&
          2.0_wp*d2_rhs_rhov_dxdy_0_0_nxp2p0m1nyp2p0k_nxp2p0m1nyp2p0m1k+&
          0.5_wp*d2_rhs_rhov_dxdy_0_0_nxp2p0m1nyp2p0k_nxp2p0m1nyp2p0m2k

d2_rhs_rhov_dxdy_0_0_nxp2p0m1nyp2p0k = d2_rhs_rhov_dxdy_0_0_nxp2p0m1nyp2p0k*param_float(2)

d2_rhs_rhov_dxdy_0_0_nxp2p0m2nyp2p0k_nxp2p0m2nyp2p0p0k = q(nx+2+0-2,ny+2+0+0,indvars(2))

d2_rhs_rhov_dxdy_0_0_nxp2p0m2nyp2p0k_nxp2p0m2nyp2p0m1k = q(nx+2+0-2,ny+2+0-1,indvars(2))

d2_rhs_rhov_dxdy_0_0_nxp2p0m2nyp2p0k_nxp2p0m2nyp2p0m2k = q(nx+2+0-2,ny+2+0-2,indvars(2))

d2_rhs_rhov_dxdy_0_0_nxp2p0m2nyp2p0k = 1.5_wp*d2_rhs_rhov_dxdy_0_0_nxp2p0m2nyp2p0k_nxp2p0m2nyp2p0p0k-&
          2.0_wp*d2_rhs_rhov_dxdy_0_0_nxp2p0m2nyp2p0k_nxp2p0m2nyp2p0m1k+&
          0.5_wp*d2_rhs_rhov_dxdy_0_0_nxp2p0m2nyp2p0k_nxp2p0m2nyp2p0m2k

d2_rhs_rhov_dxdy_0_0_nxp2p0m2nyp2p0k = d2_rhs_rhov_dxdy_0_0_nxp2p0m2nyp2p0k*param_float(2)

d1_rhs_rhov_dx_0_nxp2p0p0nyp2p0k = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+2+0+0,ny+2+0,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2+0+0,ny+2+0,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+2+0+0,ny+2+0,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(nx+2+0+0,ny+2+0,indvarsst(11))*(d2_rhs_rhov_dxdy_0_0_nxp2p0p0nyp2p0k)+&
                    qst(nx+2+0+0,ny+2+0,indvarsst(10))*(d2_rhs_rhov_dxdx_0_0_nxp2p0p0nyp2p0k))

d1_rhs_rhov_dx_0_nxp2p0m1nyp2p0k = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+2+0-1,ny+2+0,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2+0-1,ny+2+0,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+2+0-1,ny+2+0,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(nx+2+0-1,ny+2+0,indvarsst(11))*(d2_rhs_rhov_dxdy_0_0_nxp2p0m1nyp2p0k)+&
                    qst(nx+2+0-1,ny+2+0,indvarsst(10))*(d2_rhs_rhov_dxdx_0_0_nxp2p0m1nyp2p0k))

d1_rhs_rhov_dx_0_nxp2p0m2nyp2p0k = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+2+0-2,ny+2+0,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2+0-2,ny+2+0,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+2+0-2,ny+2+0,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(nx+2+0-2,ny+2+0,indvarsst(11))*(d2_rhs_rhov_dxdy_0_0_nxp2p0m2nyp2p0k)+&
                    qst(nx+2+0-2,ny+2+0,indvarsst(10))*(d2_rhs_rhov_dxdx_0_0_nxp2p0m2nyp2p0k))

d1_rhs_rhov_dx_0_nxp2p0nyp2p0k = 1.5_wp*d1_rhs_rhov_dx_0_nxp2p0p0nyp2p0k-&
          2.0_wp*d1_rhs_rhov_dx_0_nxp2p0m1nyp2p0k+&
          0.5_wp*d1_rhs_rhov_dx_0_nxp2p0m2nyp2p0k

d1_rhs_rhov_dx_0_nxp2p0nyp2p0k = d1_rhs_rhov_dx_0_nxp2p0nyp2p0k*param_float(1)

d1_rhs_rhov_dy_4_nxp2p0nyp2p0p0k = q(nx+2+0,ny+2+0+0,indvars(2))

d1_rhs_rhov_dy_4_nxp2p0nyp2p0m1k = q(nx+2+0,ny+2+0-1,indvars(2))

d1_rhs_rhov_dy_4_nxp2p0nyp2p0m2k = q(nx+2+0,ny+2+0-2,indvars(2))

d1_rhs_rhov_dy_4_nxp2p0nyp2p0k = 1.5_wp*d1_rhs_rhov_dy_4_nxp2p0nyp2p0p0k-&
          2.0_wp*d1_rhs_rhov_dy_4_nxp2p0nyp2p0m1k+&
          0.5_wp*d1_rhs_rhov_dy_4_nxp2p0nyp2p0m2k

d1_rhs_rhov_dy_4_nxp2p0nyp2p0k = d1_rhs_rhov_dy_4_nxp2p0nyp2p0k*param_float(2)

d1_rhs_rhov_dy_5_nxp2p0nyp2p0p0k = q(nx+2+0,ny+2+0+0,indvars(1))*q(nx+2+0,ny+2+0+0,indvars(3))*q(nx+2+0,ny+2+0+0,indvars(3))

d1_rhs_rhov_dy_5_nxp2p0nyp2p0m1k = q(nx+2+0,ny+2+0-1,indvars(1))*q(nx+2+0,ny+2+0-1,indvars(3))*q(nx+2+0,ny+2+0-1,indvars(3))

d1_rhs_rhov_dy_5_nxp2p0nyp2p0m2k = q(nx+2+0,ny+2+0-2,indvars(1))*q(nx+2+0,ny+2+0-2,indvars(3))*q(nx+2+0,ny+2+0-2,indvars(3))

d1_rhs_rhov_dy_5_nxp2p0nyp2p0k = 1.5_wp*d1_rhs_rhov_dy_5_nxp2p0nyp2p0p0k-&
          2.0_wp*d1_rhs_rhov_dy_5_nxp2p0nyp2p0m1k+&
          0.5_wp*d1_rhs_rhov_dy_5_nxp2p0nyp2p0m2k

d1_rhs_rhov_dy_5_nxp2p0nyp2p0k = d1_rhs_rhov_dy_5_nxp2p0nyp2p0k*param_float(2)

d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0p0k_nxp2p0p0nyp2p0p0k = q(nx+2+0+0,ny+2+0+0,indvars(2))

d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0p0k_nxp2p0m1nyp2p0p0k = q(nx+2+0-1,ny+2+0+0,indvars(2))

d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0p0k_nxp2p0m2nyp2p0p0k = q(nx+2+0-2,ny+2+0+0,indvars(2))

d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0p0k = 1.5_wp*d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0p0k_nxp2p0p0nyp2p0p0k-&
          2.0_wp*d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0p0k_nxp2p0m1nyp2p0p0k+&
          0.5_wp*d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0p0k_nxp2p0m2nyp2p0p0k

d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0p0k = d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0p0k*param_float(1)

d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0m1k_nxp2p0p0nyp2p0m1k = q(nx+2+0+0,ny+2+0-1,indvars(2))

d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0m1k_nxp2p0m1nyp2p0m1k = q(nx+2+0-1,ny+2+0-1,indvars(2))

d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0m1k_nxp2p0m2nyp2p0m1k = q(nx+2+0-2,ny+2+0-1,indvars(2))

d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0m1k = 1.5_wp*d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0m1k_nxp2p0p0nyp2p0m1k-&
          2.0_wp*d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0m1k_nxp2p0m1nyp2p0m1k+&
          0.5_wp*d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0m1k_nxp2p0m2nyp2p0m1k

d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0m1k = d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0m1k*param_float(1)

d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0m2k_nxp2p0p0nyp2p0m2k = q(nx+2+0+0,ny+2+0-2,indvars(2))

d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0m2k_nxp2p0m1nyp2p0m2k = q(nx+2+0-1,ny+2+0-2,indvars(2))

d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0m2k_nxp2p0m2nyp2p0m2k = q(nx+2+0-2,ny+2+0-2,indvars(2))

d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0m2k = 1.5_wp*d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0m2k_nxp2p0p0nyp2p0m2k-&
          2.0_wp*d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0m2k_nxp2p0m1nyp2p0m2k+&
          0.5_wp*d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0m2k_nxp2p0m2nyp2p0m2k

d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0m2k = d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0m2k*param_float(1)

d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0p0k_nxp2p0nyp2p0p0p0k = q(nx+2+0,ny+2+0+0+0,indvars(3))

d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0p0k_nxp2p0nyp2p0p0m1k = q(nx+2+0,ny+2+0+0-1,indvars(3))

d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0p0k_nxp2p0nyp2p0p0m2k = q(nx+2+0,ny+2+0+0-2,indvars(3))

d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0p0k = 1.5_wp*d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0p0k_nxp2p0nyp2p0p0p0k-&
          2.0_wp*d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0p0k_nxp2p0nyp2p0p0m1k+&
          0.5_wp*d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0p0k_nxp2p0nyp2p0p0m2k

d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0p0k = d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0p0k*param_float(2)

d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0m1k_nxp2p0nyp2p0m1m1k = q(nx+2+0,ny+2+0-1-1,indvars(3))

d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0m1k_nxp2p0nyp2p0m1p1k = q(nx+2+0,ny+2+0-1+1,indvars(3))

d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0m1k = -&
          0.5_wp*d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0m1k_nxp2p0nyp2p0m1m1k+&
          0.5_wp*d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0m1k_nxp2p0nyp2p0m1p1k

d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0m1k = d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0m1k*param_float(2)

d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0m2k_nxp2p0nyp2p0m2m1k = q(nx+2+0,ny+2+0-2-1,indvars(3))

d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0m2k_nxp2p0nyp2p0m2p1k = q(nx+2+0,ny+2+0-2+1,indvars(3))

d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0m2k = -&
          0.5_wp*d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0m2k_nxp2p0nyp2p0m2m1k+&
          0.5_wp*d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0m2k_nxp2p0nyp2p0m2p1k

d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0m2k = d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0m2k*param_float(2)

d1_rhs_rhov_dy_6_nxp2p0nyp2p0p0k = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+2+0,ny+2+0+0,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2+0,ny+2+0+0,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+2+0,ny+2+0+0,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(nx+2+0,ny+2+0+0,indvarsst(11))*(d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0p0k)-&
                    2.0_wp/3.0_wp*(qst(nx+2+0,ny+2+0+0,indvarsst(10))*(d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0p0k)+&
                    qst(nx+2+0,ny+2+0+0,indvarsst(11))*(d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0p0k)))

d1_rhs_rhov_dy_6_nxp2p0nyp2p0m1k = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+2+0,ny+2+0-1,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2+0,ny+2+0-1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+2+0,ny+2+0-1,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(nx+2+0,ny+2+0-1,indvarsst(11))*(d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0m1k)-&
                    2.0_wp/3.0_wp*(qst(nx+2+0,ny+2+0-1,indvarsst(10))*(d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0m1k)+&
                    qst(nx+2+0,ny+2+0-1,indvarsst(11))*(d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0m1k)))

d1_rhs_rhov_dy_6_nxp2p0nyp2p0m2k = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+2+0,ny+2+0-2,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2+0,ny+2+0-2,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+2+0,ny+2+0-2,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(nx+2+0,ny+2+0-2,indvarsst(11))*(d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0m2k)-&
                    2.0_wp/3.0_wp*(qst(nx+2+0,ny+2+0-2,indvarsst(10))*(d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0m2k)+&
                    qst(nx+2+0,ny+2+0-2,indvarsst(11))*(d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0m2k)))

d1_rhs_rhov_dy_6_nxp2p0nyp2p0k = 1.5_wp*d1_rhs_rhov_dy_6_nxp2p0nyp2p0p0k-&
          2.0_wp*d1_rhs_rhov_dy_6_nxp2p0nyp2p0m1k+&
          0.5_wp*d1_rhs_rhov_dy_6_nxp2p0nyp2p0m2k

d1_rhs_rhov_dy_6_nxp2p0nyp2p0k = d1_rhs_rhov_dy_6_nxp2p0nyp2p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None d(rho v)/dt ***********
!                                                           
!***********************************************************


rhs(nx+2+0,ny+2+0,indvars(3)) =   -  ( q(nx+2+0,ny+2+0,indvars(3))*(1.0_wp/(param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5**2*((q(nx+2+0,ny+2+0,indvars(3))*(-&
                    (param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5**2*d1_rhs_rho_dy_1_nxp2p0nyp2p0k+&
                    1.0_wp*d1_rhs_rho_dy_2_nxp2p0nyp2p0k))*qst(nx+2+0,ny+2+0,indvarsst(11))+&
                    0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5+&
                    q(nx+2+0,ny+2+0,indvars(3)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5*q(nx+2+0,ny+2+0,indvars(1))*d1_rhs_rho_dy_3_nxp2p0nyp2p0k+&
                    1.0_wp*d1_rhs_rho_dy_2_nxp2p0nyp2p0k))*qst(nx+2+0,ny+2+0,indvarsst(11))+&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5*(1-&
                    qface_j(nx+2+0,1)*qface_j(nx+2+0,1))/param_float(20 + 5)*((param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))-&
                    param_float(26 + 5)))))+&
                    q(nx+2+0,ny+2+0,indvars(1))*((1.0_wp*q(nx+2+0,ny+2+0,indvars(3))*d1_rhs_rhov_dy_4_nxp2p0nyp2p0k)*qst(nx+2+0,ny+2+0,indvarsst(11)))+&
                    d1_rhs_rhov_dy_5_nxp2p0nyp2p0k+&
                    qst(nx+2+0,ny+2+0,indvarsst(10))*(d1_rhs_rhov_dx_0_nxp2p0nyp2p0k)+&
                    qst(nx+2+0,ny+2+0,indvarsst(11))*(d1_rhs_rhov_dy_6_nxp2p0nyp2p0k) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None d(rho et)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! 0.5_wp*(u**2+v**2)*(1.0_wp/c**2*((v*(-c**2*[rho]_1y+1.0_wp*[p]_1y))*deltayI+0.5_wp*(((c+v)*(c*rho*[v]_1y+1.0_wp*[p]_1y))*deltayI+esse*c*(1-M_jmax*M_jmax)/L_ref*(p-P0))))+1.0_wp/gamma_m1*(0.5_wp*(((c+v)*(c*rho*[v]_1y+1.0_wp*[p]_1y))*deltayI+esse*c*(1-M_jmax*M_jmax)/L_ref*(p-P0)))+rho*u*(1.0_wp/(2*rho*c)*(((c+v)*(c*rho*[v]_1y+1.0_wp*[p]_1y))*deltayI-esse*c*(1-M_jmax*M_jmax)/L_ref*(p-P0)))+rho*v*((1.0_wp*v*[u]_1y)*deltayI)+[(rho*et+p)*v]_1y+deltaxI*([-kappa*deltaxI*({T}_1x)-u*(visc_t*(2.0_wp*deltaxI*({u}_1x)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y))))-v*(visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x)))]_1x)+deltayI*([-kappa*deltayI*({T}_1y)-u*(visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x)))-v*(visc_t*(2.0_wp*deltayI*({v}_1y)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y))))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_rhs_et_dxdx_0_0_nxp2p0p0nyp2p0k_nxp2p0p0p0nyp2p0k = ((q(nx+2+0+0+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0+0,ny+2+0,indvars(2))*q(nx+2+0+0+0,ny+2+0,indvars(2))+&
                    q(nx+2+0+0+0,ny+2+0,indvars(3))*q(nx+2+0+0+0,ny+2+0,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dxdx_0_0_nxp2p0p0nyp2p0k_nxp2p0p0m1nyp2p0k = ((q(nx+2+0+0-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0-1,ny+2+0,indvars(2))*q(nx+2+0+0-1,ny+2+0,indvars(2))+&
                    q(nx+2+0+0-1,ny+2+0,indvars(3))*q(nx+2+0+0-1,ny+2+0,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dxdx_0_0_nxp2p0p0nyp2p0k_nxp2p0p0m2nyp2p0k = ((q(nx+2+0+0-2,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0+0-2,ny+2+0,indvars(2))*q(nx+2+0+0-2,ny+2+0,indvars(2))+&
                    q(nx+2+0+0-2,ny+2+0,indvars(3))*q(nx+2+0+0-2,ny+2+0,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dxdx_0_0_nxp2p0p0nyp2p0k = 1.5_wp*d2_rhs_et_dxdx_0_0_nxp2p0p0nyp2p0k_nxp2p0p0p0nyp2p0k-&
          2.0_wp*d2_rhs_et_dxdx_0_0_nxp2p0p0nyp2p0k_nxp2p0p0m1nyp2p0k+&
          0.5_wp*d2_rhs_et_dxdx_0_0_nxp2p0p0nyp2p0k_nxp2p0p0m2nyp2p0k

d2_rhs_et_dxdx_0_0_nxp2p0p0nyp2p0k = d2_rhs_et_dxdx_0_0_nxp2p0p0nyp2p0k*param_float(1)

d2_rhs_et_dxdx_0_0_nxp2p0m1nyp2p0k_nxp2p0m1m1nyp2p0k = ((q(nx+2+0-1-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1-1,ny+2+0,indvars(2))*q(nx+2+0-1-1,ny+2+0,indvars(2))+&
                    q(nx+2+0-1-1,ny+2+0,indvars(3))*q(nx+2+0-1-1,ny+2+0,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dxdx_0_0_nxp2p0m1nyp2p0k_nxp2p0m1p1nyp2p0k = ((q(nx+2+0-1+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-1+1,ny+2+0,indvars(2))*q(nx+2+0-1+1,ny+2+0,indvars(2))+&
                    q(nx+2+0-1+1,ny+2+0,indvars(3))*q(nx+2+0-1+1,ny+2+0,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dxdx_0_0_nxp2p0m1nyp2p0k = -&
          0.5_wp*d2_rhs_et_dxdx_0_0_nxp2p0m1nyp2p0k_nxp2p0m1m1nyp2p0k+&
          0.5_wp*d2_rhs_et_dxdx_0_0_nxp2p0m1nyp2p0k_nxp2p0m1p1nyp2p0k

d2_rhs_et_dxdx_0_0_nxp2p0m1nyp2p0k = d2_rhs_et_dxdx_0_0_nxp2p0m1nyp2p0k*param_float(1)

d2_rhs_et_dxdx_0_0_nxp2p0m2nyp2p0k_nxp2p0m2m1nyp2p0k = ((q(nx+2+0-2-1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2-1,ny+2+0,indvars(2))*q(nx+2+0-2-1,ny+2+0,indvars(2))+&
                    q(nx+2+0-2-1,ny+2+0,indvars(3))*q(nx+2+0-2-1,ny+2+0,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dxdx_0_0_nxp2p0m2nyp2p0k_nxp2p0m2p1nyp2p0k = ((q(nx+2+0-2+1,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0-2+1,ny+2+0,indvars(2))*q(nx+2+0-2+1,ny+2+0,indvars(2))+&
                    q(nx+2+0-2+1,ny+2+0,indvars(3))*q(nx+2+0-2+1,ny+2+0,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dxdx_0_0_nxp2p0m2nyp2p0k = -&
          0.5_wp*d2_rhs_et_dxdx_0_0_nxp2p0m2nyp2p0k_nxp2p0m2m1nyp2p0k+&
          0.5_wp*d2_rhs_et_dxdx_0_0_nxp2p0m2nyp2p0k_nxp2p0m2p1nyp2p0k

d2_rhs_et_dxdx_0_0_nxp2p0m2nyp2p0k = d2_rhs_et_dxdx_0_0_nxp2p0m2nyp2p0k*param_float(1)

d1_rhs_et_dx_0_nxp2p0p0nyp2p0k = -param_float(2 + 5)*qst(nx+2+0+0,ny+2+0,indvarsst(10))*(d2_rhs_et_dxdx_0_0_nxp2p0p0nyp2p0k)-&
                    q(nx+2+0+0,ny+2+0,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+2+0+0,ny+2+0,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2+0+0,ny+2+0,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+2+0+0,ny+2+0,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(nx+2+0+0,ny+2+0,indvarsst(10))*(d2_rhs_rhou_dxdx_0_0_nxp2p0p0nyp2p0k)-&
                    2.0_wp/3.0_wp*(qst(nx+2+0+0,ny+2+0,indvarsst(10))*(d2_rhs_rhou_dxdx_0_0_nxp2p0p0nyp2p0k)+&
                    qst(nx+2+0+0,ny+2+0,indvarsst(11))*(d2_rhs_rhou_dxdy_0_0_nxp2p0p0nyp2p0k))))-&
                    q(nx+2+0+0,ny+2+0,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+2+0+0,ny+2+0,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2+0+0,ny+2+0,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+2+0+0,ny+2+0,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(nx+2+0+0,ny+2+0,indvarsst(11))*(d2_rhs_rhov_dxdy_0_0_nxp2p0p0nyp2p0k)+&
                    qst(nx+2+0+0,ny+2+0,indvarsst(10))*(d2_rhs_rhov_dxdx_0_0_nxp2p0p0nyp2p0k)))

d1_rhs_et_dx_0_nxp2p0m1nyp2p0k = -param_float(2 + 5)*qst(nx+2+0-1,ny+2+0,indvarsst(10))*(d2_rhs_et_dxdx_0_0_nxp2p0m1nyp2p0k)-&
                    q(nx+2+0-1,ny+2+0,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+2+0-1,ny+2+0,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2+0-1,ny+2+0,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+2+0-1,ny+2+0,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(nx+2+0-1,ny+2+0,indvarsst(10))*(d2_rhs_rhou_dxdx_0_0_nxp2p0m1nyp2p0k)-&
                    2.0_wp/3.0_wp*(qst(nx+2+0-1,ny+2+0,indvarsst(10))*(d2_rhs_rhou_dxdx_0_0_nxp2p0m1nyp2p0k)+&
                    qst(nx+2+0-1,ny+2+0,indvarsst(11))*(d2_rhs_rhou_dxdy_0_0_nxp2p0m1nyp2p0k))))-&
                    q(nx+2+0-1,ny+2+0,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+2+0-1,ny+2+0,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2+0-1,ny+2+0,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+2+0-1,ny+2+0,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(nx+2+0-1,ny+2+0,indvarsst(11))*(d2_rhs_rhov_dxdy_0_0_nxp2p0m1nyp2p0k)+&
                    qst(nx+2+0-1,ny+2+0,indvarsst(10))*(d2_rhs_rhov_dxdx_0_0_nxp2p0m1nyp2p0k)))

d1_rhs_et_dx_0_nxp2p0m2nyp2p0k = -param_float(2 + 5)*qst(nx+2+0-2,ny+2+0,indvarsst(10))*(d2_rhs_et_dxdx_0_0_nxp2p0m2nyp2p0k)-&
                    q(nx+2+0-2,ny+2+0,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+2+0-2,ny+2+0,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2+0-2,ny+2+0,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+2+0-2,ny+2+0,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(nx+2+0-2,ny+2+0,indvarsst(10))*(d2_rhs_rhou_dxdx_0_0_nxp2p0m2nyp2p0k)-&
                    2.0_wp/3.0_wp*(qst(nx+2+0-2,ny+2+0,indvarsst(10))*(d2_rhs_rhou_dxdx_0_0_nxp2p0m2nyp2p0k)+&
                    qst(nx+2+0-2,ny+2+0,indvarsst(11))*(d2_rhs_rhou_dxdy_0_0_nxp2p0m2nyp2p0k))))-&
                    q(nx+2+0-2,ny+2+0,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+2+0-2,ny+2+0,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2+0-2,ny+2+0,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+2+0-2,ny+2+0,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(nx+2+0-2,ny+2+0,indvarsst(11))*(d2_rhs_rhov_dxdy_0_0_nxp2p0m2nyp2p0k)+&
                    qst(nx+2+0-2,ny+2+0,indvarsst(10))*(d2_rhs_rhov_dxdx_0_0_nxp2p0m2nyp2p0k)))

d1_rhs_et_dx_0_nxp2p0nyp2p0k = 1.5_wp*d1_rhs_et_dx_0_nxp2p0p0nyp2p0k-&
          2.0_wp*d1_rhs_et_dx_0_nxp2p0m1nyp2p0k+&
          0.5_wp*d1_rhs_et_dx_0_nxp2p0m2nyp2p0k

d1_rhs_et_dx_0_nxp2p0nyp2p0k = d1_rhs_et_dx_0_nxp2p0nyp2p0k*param_float(1)

d1_rhs_et_dy_9_nxp2p0nyp2p0p0k = (q(nx+2+0,ny+2+0+0,indvars(1))*q(nx+2+0,ny+2+0+0,indvars(4))+&
                    (param_float(3 + 5))*q(nx+2+0,ny+2+0+0,indvars(1))*((q(nx+2+0,ny+2+0+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0+0,indvars(2))*q(nx+2+0,ny+2+0+0,indvars(2))+&
                    q(nx+2+0,ny+2+0+0,indvars(3))*q(nx+2+0,ny+2+0+0,indvars(3))))))*q(nx+2+0,ny+2+0+0,indvars(3))

d1_rhs_et_dy_9_nxp2p0nyp2p0m1k = (q(nx+2+0,ny+2+0-1,indvars(1))*q(nx+2+0,ny+2+0-1,indvars(4))+&
                    (param_float(3 + 5))*q(nx+2+0,ny+2+0-1,indvars(1))*((q(nx+2+0,ny+2+0-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0-1,indvars(2))*q(nx+2+0,ny+2+0-1,indvars(2))+&
                    q(nx+2+0,ny+2+0-1,indvars(3))*q(nx+2+0,ny+2+0-1,indvars(3))))))*q(nx+2+0,ny+2+0-1,indvars(3))

d1_rhs_et_dy_9_nxp2p0nyp2p0m2k = (q(nx+2+0,ny+2+0-2,indvars(1))*q(nx+2+0,ny+2+0-2,indvars(4))+&
                    (param_float(3 + 5))*q(nx+2+0,ny+2+0-2,indvars(1))*((q(nx+2+0,ny+2+0-2,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0-2,indvars(2))*q(nx+2+0,ny+2+0-2,indvars(2))+&
                    q(nx+2+0,ny+2+0-2,indvars(3))*q(nx+2+0,ny+2+0-2,indvars(3))))))*q(nx+2+0,ny+2+0-2,indvars(3))

d1_rhs_et_dy_9_nxp2p0nyp2p0k = 1.5_wp*d1_rhs_et_dy_9_nxp2p0nyp2p0p0k-&
          2.0_wp*d1_rhs_et_dy_9_nxp2p0nyp2p0m1k+&
          0.5_wp*d1_rhs_et_dy_9_nxp2p0nyp2p0m2k

d1_rhs_et_dy_9_nxp2p0nyp2p0k = d1_rhs_et_dy_9_nxp2p0nyp2p0k*param_float(2)

d2_rhs_et_dydy_10_0_nxp2p0nyp2p0p0k_nxp2p0nyp2p0p0p0k = ((q(nx+2+0,ny+2+0+0+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0+0+0,indvars(2))*q(nx+2+0,ny+2+0+0+0,indvars(2))+&
                    q(nx+2+0,ny+2+0+0+0,indvars(3))*q(nx+2+0,ny+2+0+0+0,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_10_0_nxp2p0nyp2p0p0k_nxp2p0nyp2p0p0m1k = ((q(nx+2+0,ny+2+0+0-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0+0-1,indvars(2))*q(nx+2+0,ny+2+0+0-1,indvars(2))+&
                    q(nx+2+0,ny+2+0+0-1,indvars(3))*q(nx+2+0,ny+2+0+0-1,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_10_0_nxp2p0nyp2p0p0k_nxp2p0nyp2p0p0m2k = ((q(nx+2+0,ny+2+0+0-2,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0+0-2,indvars(2))*q(nx+2+0,ny+2+0+0-2,indvars(2))+&
                    q(nx+2+0,ny+2+0+0-2,indvars(3))*q(nx+2+0,ny+2+0+0-2,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_10_0_nxp2p0nyp2p0p0k = 1.5_wp*d2_rhs_et_dydy_10_0_nxp2p0nyp2p0p0k_nxp2p0nyp2p0p0p0k-&
          2.0_wp*d2_rhs_et_dydy_10_0_nxp2p0nyp2p0p0k_nxp2p0nyp2p0p0m1k+&
          0.5_wp*d2_rhs_et_dydy_10_0_nxp2p0nyp2p0p0k_nxp2p0nyp2p0p0m2k

d2_rhs_et_dydy_10_0_nxp2p0nyp2p0p0k = d2_rhs_et_dydy_10_0_nxp2p0nyp2p0p0k*param_float(2)

d2_rhs_et_dydy_10_0_nxp2p0nyp2p0m1k_nxp2p0nyp2p0m1m1k = ((q(nx+2+0,ny+2+0-1-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0-1-1,indvars(2))*q(nx+2+0,ny+2+0-1-1,indvars(2))+&
                    q(nx+2+0,ny+2+0-1-1,indvars(3))*q(nx+2+0,ny+2+0-1-1,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_10_0_nxp2p0nyp2p0m1k_nxp2p0nyp2p0m1p1k = ((q(nx+2+0,ny+2+0-1+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0-1+1,indvars(2))*q(nx+2+0,ny+2+0-1+1,indvars(2))+&
                    q(nx+2+0,ny+2+0-1+1,indvars(3))*q(nx+2+0,ny+2+0-1+1,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_10_0_nxp2p0nyp2p0m1k = -&
          0.5_wp*d2_rhs_et_dydy_10_0_nxp2p0nyp2p0m1k_nxp2p0nyp2p0m1m1k+&
          0.5_wp*d2_rhs_et_dydy_10_0_nxp2p0nyp2p0m1k_nxp2p0nyp2p0m1p1k

d2_rhs_et_dydy_10_0_nxp2p0nyp2p0m1k = d2_rhs_et_dydy_10_0_nxp2p0nyp2p0m1k*param_float(2)

d2_rhs_et_dydy_10_0_nxp2p0nyp2p0m2k_nxp2p0nyp2p0m2m1k = ((q(nx+2+0,ny+2+0-2-1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0-2-1,indvars(2))*q(nx+2+0,ny+2+0-2-1,indvars(2))+&
                    q(nx+2+0,ny+2+0-2-1,indvars(3))*q(nx+2+0,ny+2+0-2-1,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_10_0_nxp2p0nyp2p0m2k_nxp2p0nyp2p0m2p1k = ((q(nx+2+0,ny+2+0-2+1,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0-2+1,indvars(2))*q(nx+2+0,ny+2+0-2+1,indvars(2))+&
                    q(nx+2+0,ny+2+0-2+1,indvars(3))*q(nx+2+0,ny+2+0-2+1,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_10_0_nxp2p0nyp2p0m2k = -&
          0.5_wp*d2_rhs_et_dydy_10_0_nxp2p0nyp2p0m2k_nxp2p0nyp2p0m2m1k+&
          0.5_wp*d2_rhs_et_dydy_10_0_nxp2p0nyp2p0m2k_nxp2p0nyp2p0m2p1k

d2_rhs_et_dydy_10_0_nxp2p0nyp2p0m2k = d2_rhs_et_dydy_10_0_nxp2p0nyp2p0m2k*param_float(2)

d1_rhs_et_dy_10_nxp2p0nyp2p0p0k = -param_float(2 + 5)*qst(nx+2+0,ny+2+0+0,indvarsst(11))*(d2_rhs_et_dydy_10_0_nxp2p0nyp2p0p0k)-&
                    q(nx+2+0,ny+2+0+0,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+2+0,ny+2+0+0,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2+0,ny+2+0+0,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+2+0,ny+2+0+0,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(nx+2+0,ny+2+0+0,indvarsst(11))*(d2_rhs_rhou_dydy_7_0_nxp2p0nyp2p0p0k)+&
                    qst(nx+2+0,ny+2+0+0,indvarsst(10))*(d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0p0k)))-&
                    q(nx+2+0,ny+2+0+0,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+2+0,ny+2+0+0,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2+0,ny+2+0+0,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+2+0,ny+2+0+0,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(nx+2+0,ny+2+0+0,indvarsst(11))*(d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0p0k)-&
                    2.0_wp/3.0_wp*(qst(nx+2+0,ny+2+0+0,indvarsst(10))*(d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0p0k)+&
                    qst(nx+2+0,ny+2+0+0,indvarsst(11))*(d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0p0k))))

d1_rhs_et_dy_10_nxp2p0nyp2p0m1k = -param_float(2 + 5)*qst(nx+2+0,ny+2+0-1,indvarsst(11))*(d2_rhs_et_dydy_10_0_nxp2p0nyp2p0m1k)-&
                    q(nx+2+0,ny+2+0-1,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+2+0,ny+2+0-1,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2+0,ny+2+0-1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+2+0,ny+2+0-1,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(nx+2+0,ny+2+0-1,indvarsst(11))*(d2_rhs_rhou_dydy_7_0_nxp2p0nyp2p0m1k)+&
                    qst(nx+2+0,ny+2+0-1,indvarsst(10))*(d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0m1k)))-&
                    q(nx+2+0,ny+2+0-1,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+2+0,ny+2+0-1,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2+0,ny+2+0-1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+2+0,ny+2+0-1,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(nx+2+0,ny+2+0-1,indvarsst(11))*(d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0m1k)-&
                    2.0_wp/3.0_wp*(qst(nx+2+0,ny+2+0-1,indvarsst(10))*(d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0m1k)+&
                    qst(nx+2+0,ny+2+0-1,indvarsst(11))*(d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0m1k))))

d1_rhs_et_dy_10_nxp2p0nyp2p0m2k = -param_float(2 + 5)*qst(nx+2+0,ny+2+0-2,indvarsst(11))*(d2_rhs_et_dydy_10_0_nxp2p0nyp2p0m2k)-&
                    q(nx+2+0,ny+2+0-2,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+2+0,ny+2+0-2,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2+0,ny+2+0-2,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+2+0,ny+2+0-2,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(nx+2+0,ny+2+0-2,indvarsst(11))*(d2_rhs_rhou_dydy_7_0_nxp2p0nyp2p0m2k)+&
                    qst(nx+2+0,ny+2+0-2,indvarsst(10))*(d2_rhs_rhou_dydx_7_0_nxp2p0nyp2p0m2k)))-&
                    q(nx+2+0,ny+2+0-2,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+2+0,ny+2+0-2,indvars(5))/1.0_wp)**3.0_wp/((q(nx+2+0,ny+2+0-2,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+2+0,ny+2+0-2,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(nx+2+0,ny+2+0-2,indvarsst(11))*(d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0m2k)-&
                    2.0_wp/3.0_wp*(qst(nx+2+0,ny+2+0-2,indvarsst(10))*(d2_rhs_rhov_dydx_6_0_nxp2p0nyp2p0m2k)+&
                    qst(nx+2+0,ny+2+0-2,indvarsst(11))*(d2_rhs_rhov_dydy_6_0_nxp2p0nyp2p0m2k))))

d1_rhs_et_dy_10_nxp2p0nyp2p0k = 1.5_wp*d1_rhs_et_dy_10_nxp2p0nyp2p0p0k-&
          2.0_wp*d1_rhs_et_dy_10_nxp2p0nyp2p0m1k+&
          0.5_wp*d1_rhs_et_dy_10_nxp2p0nyp2p0m2k

d1_rhs_et_dy_10_nxp2p0nyp2p0k = d1_rhs_et_dy_10_nxp2p0nyp2p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None d(rho et)/dt **********
!                                                           
!***********************************************************


rhs(nx+2+0,ny+2+0,indvars(4)) =   -  ( 0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))**2+&
                    q(nx+2+0,ny+2+0,indvars(3))**2)*(1.0_wp/(param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5**2*((q(nx+2+0,ny+2+0,indvars(3))*(-&
                    (param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5**2*d1_rhs_rho_dy_1_nxp2p0nyp2p0k+&
                    1.0_wp*d1_rhs_rho_dy_2_nxp2p0nyp2p0k))*qst(nx+2+0,ny+2+0,indvarsst(11))+&
                    0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5+&
                    q(nx+2+0,ny+2+0,indvars(3)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5*q(nx+2+0,ny+2+0,indvars(1))*d1_rhs_rho_dy_3_nxp2p0nyp2p0k+&
                    1.0_wp*d1_rhs_rho_dy_2_nxp2p0nyp2p0k))*qst(nx+2+0,ny+2+0,indvarsst(11))+&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5*(1-&
                    qface_j(nx+2+0,1)*qface_j(nx+2+0,1))/param_float(20 + 5)*((param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))-&
                    param_float(26 + 5)))))+&
                    1.0_wp/param_float(3 + 5)*(0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5+&
                    q(nx+2+0,ny+2+0,indvars(3)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5*q(nx+2+0,ny+2+0,indvars(1))*d1_rhs_rho_dy_3_nxp2p0nyp2p0k+&
                    1.0_wp*d1_rhs_rho_dy_2_nxp2p0nyp2p0k))*qst(nx+2+0,ny+2+0,indvarsst(11))+&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5*(1-&
                    qface_j(nx+2+0,1)*qface_j(nx+2+0,1))/param_float(20 + 5)*((param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))-&
                    param_float(26 + 5))))+&
                    q(nx+2+0,ny+2+0,indvars(1))*q(nx+2+0,ny+2+0,indvars(2))*(1.0_wp/(2*q(nx+2+0,ny+2+0,indvars(1))*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5)*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5+&
                    q(nx+2+0,ny+2+0,indvars(3)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5*q(nx+2+0,ny+2+0,indvars(1))*d1_rhs_rho_dy_3_nxp2p0nyp2p0k+&
                    1.0_wp*d1_rhs_rho_dy_2_nxp2p0nyp2p0k))*qst(nx+2+0,ny+2+0,indvarsst(11))-&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))/q(nx+2+0,ny+2+0,indvars(1)))**0.5*(1-&
                    qface_j(nx+2+0,1)*qface_j(nx+2+0,1))/param_float(20 + 5)*((param_float(3 + 5))*q(nx+2+0,ny+2+0,indvars(1))*((q(nx+2+0,ny+2+0,indvars(4))-&
                    0.5_wp*(q(nx+2+0,ny+2+0,indvars(2))*q(nx+2+0,ny+2+0,indvars(2))+&
                    q(nx+2+0,ny+2+0,indvars(3))*q(nx+2+0,ny+2+0,indvars(3)))))-&
                    param_float(26 + 5))))+&
                    q(nx+2+0,ny+2+0,indvars(1))*q(nx+2+0,ny+2+0,indvars(3))*((1.0_wp*q(nx+2+0,ny+2+0,indvars(3))*d1_rhs_rhov_dy_4_nxp2p0nyp2p0k)*qst(nx+2+0,ny+2+0,indvarsst(11)))+&
                    d1_rhs_et_dy_9_nxp2p0nyp2p0k+&
                    qst(nx+2+0,ny+2+0,indvarsst(10))*(d1_rhs_et_dx_0_nxp2p0nyp2p0k)+&
                    qst(nx+2+0,ny+2+0,indvarsst(11))*(d1_rhs_et_dy_10_nxp2p0nyp2p0k) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None d(rho nut)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([(rho*u*nut)]_1x)+-deltaxI*([(ReI*(1.0_wp+chi)*sigmaI*deltaxI*({nut}_1x))]_1x)-ReI*Cb2*sigmaI*((deltaxI)**2*([rho*nut]_1x)*([nut]_1x))-Cb1*(1-ft2)*SS*rho*nut+ReI*(Cw1*fw-Cb1/k**2*ft2)*rho*(nut/eta)**2.0_wp
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_rhs_nut_dx_0_nxp2p0p0nyp2p0k = (q(nx+2+0+0,ny+2+0,indvars(1))*q(nx+2+0+0,ny+2+0,indvars(2))*q(nx+2+0+0,ny+2+0,indvars(5)))

d1_rhs_nut_dx_0_nxp2p0m1nyp2p0k = (q(nx+2+0-1,ny+2+0,indvars(1))*q(nx+2+0-1,ny+2+0,indvars(2))*q(nx+2+0-1,ny+2+0,indvars(5)))

d1_rhs_nut_dx_0_nxp2p0m2nyp2p0k = (q(nx+2+0-2,ny+2+0,indvars(1))*q(nx+2+0-2,ny+2+0,indvars(2))*q(nx+2+0-2,ny+2+0,indvars(5)))

d1_rhs_nut_dx_0_nxp2p0nyp2p0k = 1.5_wp*d1_rhs_nut_dx_0_nxp2p0p0nyp2p0k-&
          2.0_wp*d1_rhs_nut_dx_0_nxp2p0m1nyp2p0k+&
          0.5_wp*d1_rhs_nut_dx_0_nxp2p0m2nyp2p0k

d1_rhs_nut_dx_0_nxp2p0nyp2p0k = d1_rhs_nut_dx_0_nxp2p0nyp2p0k*param_float(1)

d2_rhs_nut_dxdx_1_0_nxp2p0p0nyp2p0k_nxp2p0p0p0nyp2p0k = q(nx+2+0+0+0,ny+2+0,indvars(5))

d2_rhs_nut_dxdx_1_0_nxp2p0p0nyp2p0k_nxp2p0p0m1nyp2p0k = q(nx+2+0+0-1,ny+2+0,indvars(5))

d2_rhs_nut_dxdx_1_0_nxp2p0p0nyp2p0k_nxp2p0p0m2nyp2p0k = q(nx+2+0+0-2,ny+2+0,indvars(5))

d2_rhs_nut_dxdx_1_0_nxp2p0p0nyp2p0k = 1.5_wp*d2_rhs_nut_dxdx_1_0_nxp2p0p0nyp2p0k_nxp2p0p0p0nyp2p0k-&
          2.0_wp*d2_rhs_nut_dxdx_1_0_nxp2p0p0nyp2p0k_nxp2p0p0m1nyp2p0k+&
          0.5_wp*d2_rhs_nut_dxdx_1_0_nxp2p0p0nyp2p0k_nxp2p0p0m2nyp2p0k

d2_rhs_nut_dxdx_1_0_nxp2p0p0nyp2p0k = d2_rhs_nut_dxdx_1_0_nxp2p0p0nyp2p0k*param_float(1)

d2_rhs_nut_dxdx_1_0_nxp2p0m1nyp2p0k_nxp2p0m1m1nyp2p0k = q(nx+2+0-1-1,ny+2+0,indvars(5))

d2_rhs_nut_dxdx_1_0_nxp2p0m1nyp2p0k_nxp2p0m1p1nyp2p0k = q(nx+2+0-1+1,ny+2+0,indvars(5))

d2_rhs_nut_dxdx_1_0_nxp2p0m1nyp2p0k = -&
          0.5_wp*d2_rhs_nut_dxdx_1_0_nxp2p0m1nyp2p0k_nxp2p0m1m1nyp2p0k+&
          0.5_wp*d2_rhs_nut_dxdx_1_0_nxp2p0m1nyp2p0k_nxp2p0m1p1nyp2p0k

d2_rhs_nut_dxdx_1_0_nxp2p0m1nyp2p0k = d2_rhs_nut_dxdx_1_0_nxp2p0m1nyp2p0k*param_float(1)

d2_rhs_nut_dxdx_1_0_nxp2p0m2nyp2p0k_nxp2p0m2m1nyp2p0k = q(nx+2+0-2-1,ny+2+0,indvars(5))

d2_rhs_nut_dxdx_1_0_nxp2p0m2nyp2p0k_nxp2p0m2p1nyp2p0k = q(nx+2+0-2+1,ny+2+0,indvars(5))

d2_rhs_nut_dxdx_1_0_nxp2p0m2nyp2p0k = -&
          0.5_wp*d2_rhs_nut_dxdx_1_0_nxp2p0m2nyp2p0k_nxp2p0m2m1nyp2p0k+&
          0.5_wp*d2_rhs_nut_dxdx_1_0_nxp2p0m2nyp2p0k_nxp2p0m2p1nyp2p0k

d2_rhs_nut_dxdx_1_0_nxp2p0m2nyp2p0k = d2_rhs_nut_dxdx_1_0_nxp2p0m2nyp2p0k*param_float(1)

d1_rhs_nut_dx_2_nxp2p0p0nyp2p0k = q(nx+2+0+0,ny+2+0,indvars(1))*q(nx+2+0+0,ny+2+0,indvars(5))

d1_rhs_nut_dx_2_nxp2p0m1nyp2p0k = q(nx+2+0-1,ny+2+0,indvars(1))*q(nx+2+0-1,ny+2+0,indvars(5))

d1_rhs_nut_dx_2_nxp2p0m2nyp2p0k = q(nx+2+0-2,ny+2+0,indvars(1))*q(nx+2+0-2,ny+2+0,indvars(5))

d1_rhs_nut_dx_2_nxp2p0nyp2p0k = 1.5_wp*d1_rhs_nut_dx_2_nxp2p0p0nyp2p0k-&
          2.0_wp*d1_rhs_nut_dx_2_nxp2p0m1nyp2p0k+&
          0.5_wp*d1_rhs_nut_dx_2_nxp2p0m2nyp2p0k

d1_rhs_nut_dx_2_nxp2p0nyp2p0k = d1_rhs_nut_dx_2_nxp2p0nyp2p0k*param_float(1)

d1_rhs_nut_dx_3_nxp2p0p0nyp2p0k = q(nx+2+0+0,ny+2+0,indvars(5))

d1_rhs_nut_dx_3_nxp2p0m1nyp2p0k = q(nx+2+0-1,ny+2+0,indvars(5))

d1_rhs_nut_dx_3_nxp2p0m2nyp2p0k = q(nx+2+0-2,ny+2+0,indvars(5))

d1_rhs_nut_dx_3_nxp2p0nyp2p0k = 1.5_wp*d1_rhs_nut_dx_3_nxp2p0p0nyp2p0k-&
          2.0_wp*d1_rhs_nut_dx_3_nxp2p0m1nyp2p0k+&
          0.5_wp*d1_rhs_nut_dx_3_nxp2p0m2nyp2p0k

d1_rhs_nut_dx_3_nxp2p0nyp2p0k = d1_rhs_nut_dx_3_nxp2p0nyp2p0k*param_float(1)

d1_rhs_nut_dx_1_nxp2p0p0nyp2p0k = (param_float(1 + 5)*(1.0_wp+&
                    (q(nx+2+0+0,ny+2+0,indvars(5))/1.0_wp))*param_float(18 + 5)*qst(nx+2+0+0,ny+2+0,indvarsst(10))*(d2_rhs_nut_dxdx_1_0_nxp2p0p0nyp2p0k))

d1_rhs_nut_dx_1_nxp2p0m1nyp2p0k = (param_float(1 + 5)*(1.0_wp+&
                    (q(nx+2+0-1,ny+2+0,indvars(5))/1.0_wp))*param_float(18 + 5)*qst(nx+2+0-1,ny+2+0,indvarsst(10))*(d2_rhs_nut_dxdx_1_0_nxp2p0m1nyp2p0k))

d1_rhs_nut_dx_1_nxp2p0m2nyp2p0k = (param_float(1 + 5)*(1.0_wp+&
                    (q(nx+2+0-2,ny+2+0,indvars(5))/1.0_wp))*param_float(18 + 5)*qst(nx+2+0-2,ny+2+0,indvarsst(10))*(d2_rhs_nut_dxdx_1_0_nxp2p0m2nyp2p0k))

d1_rhs_nut_dx_1_nxp2p0nyp2p0k = 1.5_wp*d1_rhs_nut_dx_1_nxp2p0p0nyp2p0k-&
          2.0_wp*d1_rhs_nut_dx_1_nxp2p0m1nyp2p0k+&
          0.5_wp*d1_rhs_nut_dx_1_nxp2p0m2nyp2p0k

d1_rhs_nut_dx_1_nxp2p0nyp2p0k = d1_rhs_nut_dx_1_nxp2p0nyp2p0k*param_float(1)



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None d(rho nut)/dt *********
!                                                           
!***********************************************************


rhs(nx+2+0,ny+2+0,indvars(5)) =   -  ( qst(nx+2+0,ny+2+0,indvarsst(10))*(d1_rhs_nut_dx_0_nxp2p0nyp2p0k)+&
                    -&
                    qst(nx+2+0,ny+2+0,indvarsst(10))*(d1_rhs_nut_dx_1_nxp2p0nyp2p0k)-&
                    param_float(1 + 5)*param_float(7 + 5)*param_float(18 + 5)*((qst(nx+2+0,ny+2+0,indvarsst(10)))**2*(d1_rhs_nut_dx_2_nxp2p0nyp2p0k)*(d1_rhs_nut_dx_3_nxp2p0nyp2p0k))-&
                    param_float(6 + 5)*(1-&
                    (param_float(16 + 5)*exp(-&
                    param_float(17 + 5)*(q(nx+2+0,ny+2+0,indvars(5))/1.0_wp)**2.0_wp)))*qst(nx+2+0,ny+2+0,indvarsst(13))*q(nx+2+0,ny+2+0,indvars(1))*q(nx+2+0,ny+2+0,indvars(5))+&
                    param_float(1 + 5)*(param_float(10 + 5)*qst(nx+2+0,ny+2+0,indvarsst(16))-&
                    param_float(6 + 5)/param_float(9 + 5)**2*(param_float(16 + 5)*exp(-&
                    param_float(17 + 5)*(q(nx+2+0,ny+2+0,indvars(5))/1.0_wp)**2.0_wp)))*q(nx+2+0,ny+2+0,indvars(1))*(q(nx+2+0,ny+2+0,indvars(5))/qst(nx+2+0,ny+2+0,indvarsst(2)))**2.0_wp ) 

