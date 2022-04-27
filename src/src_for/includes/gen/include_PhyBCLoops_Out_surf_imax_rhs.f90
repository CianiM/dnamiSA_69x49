

!***********************************************************
!                                                           
! Start building layers for BC : imax None None ************
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
! deltaxI*([rho*u]_1x)+deltayI*([rho*v]_1y)+(1.0_wp/c**2*((u*(-c**2*[rho]_1x+1.0_wp*[p]_1x))*deltaxI+0.5_wp*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI+esse*c*(1-M_imax*M_imax)/L_ref*(p-P0))))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_rhs_rho_dx_0_nxp4p0p0jk = q(nx+4+0+0,j,indvars(1))*q(nx+4+0+0,j,indvars(2))

d1_rhs_rho_dx_0_nxp4p0m1jk = q(nx+4+0-1,j,indvars(1))*q(nx+4+0-1,j,indvars(2))

d1_rhs_rho_dx_0_nxp4p0m2jk = q(nx+4+0-2,j,indvars(1))*q(nx+4+0-2,j,indvars(2))

d1_rhs_rho_dx_0_nxp4p0jk = 1.5_wp*d1_rhs_rho_dx_0_nxp4p0p0jk-&
          2.0_wp*d1_rhs_rho_dx_0_nxp4p0m1jk+&
          0.5_wp*d1_rhs_rho_dx_0_nxp4p0m2jk

d1_rhs_rho_dx_0_nxp4p0jk = d1_rhs_rho_dx_0_nxp4p0jk*param_float(1)

d1_rhs_rho_dx_1_nxp4p0p0jk = q(nx+4+0+0,j,indvars(1))

d1_rhs_rho_dx_1_nxp4p0m1jk = q(nx+4+0-1,j,indvars(1))

d1_rhs_rho_dx_1_nxp4p0m2jk = q(nx+4+0-2,j,indvars(1))

d1_rhs_rho_dx_1_nxp4p0jk = 1.5_wp*d1_rhs_rho_dx_1_nxp4p0p0jk-&
          2.0_wp*d1_rhs_rho_dx_1_nxp4p0m1jk+&
          0.5_wp*d1_rhs_rho_dx_1_nxp4p0m2jk

d1_rhs_rho_dx_1_nxp4p0jk = d1_rhs_rho_dx_1_nxp4p0jk*param_float(1)

d1_rhs_rho_dx_2_nxp4p0p0jk = (param_float(3 + 5))*q(nx+4+0+0,j,indvars(1))*((q(nx+4+0+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,j,indvars(2))*q(nx+4+0+0,j,indvars(2))+&
                    q(nx+4+0+0,j,indvars(3))*q(nx+4+0+0,j,indvars(3)))))

d1_rhs_rho_dx_2_nxp4p0m1jk = (param_float(3 + 5))*q(nx+4+0-1,j,indvars(1))*((q(nx+4+0-1,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,j,indvars(2))*q(nx+4+0-1,j,indvars(2))+&
                    q(nx+4+0-1,j,indvars(3))*q(nx+4+0-1,j,indvars(3)))))

d1_rhs_rho_dx_2_nxp4p0m2jk = (param_float(3 + 5))*q(nx+4+0-2,j,indvars(1))*((q(nx+4+0-2,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,j,indvars(2))*q(nx+4+0-2,j,indvars(2))+&
                    q(nx+4+0-2,j,indvars(3))*q(nx+4+0-2,j,indvars(3)))))

d1_rhs_rho_dx_2_nxp4p0jk = 1.5_wp*d1_rhs_rho_dx_2_nxp4p0p0jk-&
          2.0_wp*d1_rhs_rho_dx_2_nxp4p0m1jk+&
          0.5_wp*d1_rhs_rho_dx_2_nxp4p0m2jk

d1_rhs_rho_dx_2_nxp4p0jk = d1_rhs_rho_dx_2_nxp4p0jk*param_float(1)

d1_rhs_rho_dx_3_nxp4p0p0jk = q(nx+4+0+0,j,indvars(2))

d1_rhs_rho_dx_3_nxp4p0m1jk = q(nx+4+0-1,j,indvars(2))

d1_rhs_rho_dx_3_nxp4p0m2jk = q(nx+4+0-2,j,indvars(2))

d1_rhs_rho_dx_3_nxp4p0jk = 1.5_wp*d1_rhs_rho_dx_3_nxp4p0p0jk-&
          2.0_wp*d1_rhs_rho_dx_3_nxp4p0m1jk+&
          0.5_wp*d1_rhs_rho_dx_3_nxp4p0m2jk

d1_rhs_rho_dx_3_nxp4p0jk = d1_rhs_rho_dx_3_nxp4p0jk*param_float(1)

d1_rhs_rho_dy_0_nxp4p0jm1k = q(nx+4+0,j-1,indvars(1))*q(nx+4+0,j-1,indvars(3))

d1_rhs_rho_dy_0_nxp4p0jp1k = q(nx+4+0,j+1,indvars(1))*q(nx+4+0,j+1,indvars(3))

d1_rhs_rho_dy_0_nxp4p0jk = -&
          0.5_wp*d1_rhs_rho_dy_0_nxp4p0jm1k+&
          0.5_wp*d1_rhs_rho_dy_0_nxp4p0jp1k

d1_rhs_rho_dy_0_nxp4p0jk = d1_rhs_rho_dy_0_nxp4p0jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None d(rho)/dt **********
!                                                           
!***********************************************************


rhs(nx+4+0,j,indvars(1)) =   -  ( qst(nx+4+0,j,indvarsst(10))*(d1_rhs_rho_dx_0_nxp4p0jk)+&
                    qst(nx+4+0,j,indvarsst(11))*(d1_rhs_rho_dy_0_nxp4p0jk)+&
                    (1.0_wp/(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5**2*((q(nx+4+0,j,indvars(2))*(-&
                    (param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5**2*d1_rhs_rho_dx_1_nxp4p0jk+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p0jk))*qst(nx+4+0,j,indvarsst(10))+&
                    0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5+&
                    q(nx+4+0,j,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5*q(nx+4+0,j,indvars(1))*d1_rhs_rho_dx_3_nxp4p0jk+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p0jk))*qst(nx+4+0,j,indvarsst(10))+&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5*(1-&
                    qface_i(j,2)*qface_i(j,2))/param_float(20 + 5)*((param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))-&
                    param_float(26 + 5))))) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None d(rho u)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! u*(1.0_wp/c**2*((u*(-c**2*[rho]_1x+1.0_wp*[p]_1x))*deltaxI+0.5_wp*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI+esse*c*(1-M_imax*M_imax)/L_ref*(p-P0))))+rho*(1.0_wp/(2*rho*c)*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI-esse*c*(1-M_imax*M_imax)/L_ref*(p-P0)))+[rho*u*v]_1y+deltaxI*([-visc_t*(2.0_wp*deltaxI*({u}_1x)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y)))]_1x)+deltayI*([-visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_rhs_u_dxdx_6_0_nxp4p0p0jk_nxp4p0p0p0jk = q(nx+4+0+0+0,j,indvars(2))

d2_rhs_u_dxdx_6_0_nxp4p0p0jk_nxp4p0p0m1jk = q(nx+4+0+0-1,j,indvars(2))

d2_rhs_u_dxdx_6_0_nxp4p0p0jk_nxp4p0p0m2jk = q(nx+4+0+0-2,j,indvars(2))

d2_rhs_u_dxdx_6_0_nxp4p0p0jk = 1.5_wp*d2_rhs_u_dxdx_6_0_nxp4p0p0jk_nxp4p0p0p0jk-&
          2.0_wp*d2_rhs_u_dxdx_6_0_nxp4p0p0jk_nxp4p0p0m1jk+&
          0.5_wp*d2_rhs_u_dxdx_6_0_nxp4p0p0jk_nxp4p0p0m2jk

d2_rhs_u_dxdx_6_0_nxp4p0p0jk = d2_rhs_u_dxdx_6_0_nxp4p0p0jk*param_float(1)

d2_rhs_u_dxdx_6_0_nxp4p0m1jk_nxp4p0m1m1jk = q(nx+4+0-1-1,j,indvars(2))

d2_rhs_u_dxdx_6_0_nxp4p0m1jk_nxp4p0m1p1jk = q(nx+4+0-1+1,j,indvars(2))

d2_rhs_u_dxdx_6_0_nxp4p0m1jk = -&
          0.5_wp*d2_rhs_u_dxdx_6_0_nxp4p0m1jk_nxp4p0m1m1jk+&
          0.5_wp*d2_rhs_u_dxdx_6_0_nxp4p0m1jk_nxp4p0m1p1jk

d2_rhs_u_dxdx_6_0_nxp4p0m1jk = d2_rhs_u_dxdx_6_0_nxp4p0m1jk*param_float(1)

d2_rhs_u_dxdx_6_0_nxp4p0m2jk_nxp4p0m2m1jk = q(nx+4+0-2-1,j,indvars(2))

d2_rhs_u_dxdx_6_0_nxp4p0m2jk_nxp4p0m2p1jk = q(nx+4+0-2+1,j,indvars(2))

d2_rhs_u_dxdx_6_0_nxp4p0m2jk = -&
          0.5_wp*d2_rhs_u_dxdx_6_0_nxp4p0m2jk_nxp4p0m2m1jk+&
          0.5_wp*d2_rhs_u_dxdx_6_0_nxp4p0m2jk_nxp4p0m2p1jk

d2_rhs_u_dxdx_6_0_nxp4p0m2jk = d2_rhs_u_dxdx_6_0_nxp4p0m2jk*param_float(1)

d2_rhs_u_dxdy_6_0_nxp4p0p0jk_nxp4p0p0jm1k = q(nx+4+0+0,j-1,indvars(3))

d2_rhs_u_dxdy_6_0_nxp4p0p0jk_nxp4p0p0jp1k = q(nx+4+0+0,j+1,indvars(3))

d2_rhs_u_dxdy_6_0_nxp4p0p0jk = -&
          0.5_wp*d2_rhs_u_dxdy_6_0_nxp4p0p0jk_nxp4p0p0jm1k+&
          0.5_wp*d2_rhs_u_dxdy_6_0_nxp4p0p0jk_nxp4p0p0jp1k

d2_rhs_u_dxdy_6_0_nxp4p0p0jk = d2_rhs_u_dxdy_6_0_nxp4p0p0jk*param_float(2)

d2_rhs_u_dxdy_6_0_nxp4p0m1jk_nxp4p0m1jm1k = q(nx+4+0-1,j-1,indvars(3))

d2_rhs_u_dxdy_6_0_nxp4p0m1jk_nxp4p0m1jp1k = q(nx+4+0-1,j+1,indvars(3))

d2_rhs_u_dxdy_6_0_nxp4p0m1jk = -&
          0.5_wp*d2_rhs_u_dxdy_6_0_nxp4p0m1jk_nxp4p0m1jm1k+&
          0.5_wp*d2_rhs_u_dxdy_6_0_nxp4p0m1jk_nxp4p0m1jp1k

d2_rhs_u_dxdy_6_0_nxp4p0m1jk = d2_rhs_u_dxdy_6_0_nxp4p0m1jk*param_float(2)

d2_rhs_u_dxdy_6_0_nxp4p0m2jk_nxp4p0m2jm1k = q(nx+4+0-2,j-1,indvars(3))

d2_rhs_u_dxdy_6_0_nxp4p0m2jk_nxp4p0m2jp1k = q(nx+4+0-2,j+1,indvars(3))

d2_rhs_u_dxdy_6_0_nxp4p0m2jk = -&
          0.5_wp*d2_rhs_u_dxdy_6_0_nxp4p0m2jk_nxp4p0m2jm1k+&
          0.5_wp*d2_rhs_u_dxdy_6_0_nxp4p0m2jk_nxp4p0m2jp1k

d2_rhs_u_dxdy_6_0_nxp4p0m2jk = d2_rhs_u_dxdy_6_0_nxp4p0m2jk*param_float(2)

d1_rhs_u_dx_6_nxp4p0p0jk = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0+0,j,indvars(5))/1.0_wp*q(nx+4+0+0,j,indvars(1)))**3.0_wp/((q(nx+4+0+0,j,indvars(5))/1.0_wp*q(nx+4+0+0,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0+0,j,indvars(5))/1.0_wp*q(nx+4+0+0,j,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0+0,j,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0p0jk)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0+0,j,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0p0jk)+&
                    qst(nx+4+0+0,j,indvarsst(11))*(d2_rhs_u_dxdy_6_0_nxp4p0p0jk)))

d1_rhs_u_dx_6_nxp4p0m1jk = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0-1,j,indvars(5))/1.0_wp*q(nx+4+0-1,j,indvars(1)))**3.0_wp/((q(nx+4+0-1,j,indvars(5))/1.0_wp*q(nx+4+0-1,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0-1,j,indvars(5))/1.0_wp*q(nx+4+0-1,j,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0-1,j,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0m1jk)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0-1,j,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0m1jk)+&
                    qst(nx+4+0-1,j,indvarsst(11))*(d2_rhs_u_dxdy_6_0_nxp4p0m1jk)))

d1_rhs_u_dx_6_nxp4p0m2jk = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0-2,j,indvars(5))/1.0_wp*q(nx+4+0-2,j,indvars(1)))**3.0_wp/((q(nx+4+0-2,j,indvars(5))/1.0_wp*q(nx+4+0-2,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0-2,j,indvars(5))/1.0_wp*q(nx+4+0-2,j,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0-2,j,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0m2jk)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0-2,j,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0m2jk)+&
                    qst(nx+4+0-2,j,indvarsst(11))*(d2_rhs_u_dxdy_6_0_nxp4p0m2jk)))

d1_rhs_u_dx_6_nxp4p0jk = 1.5_wp*d1_rhs_u_dx_6_nxp4p0p0jk-&
          2.0_wp*d1_rhs_u_dx_6_nxp4p0m1jk+&
          0.5_wp*d1_rhs_u_dx_6_nxp4p0m2jk

d1_rhs_u_dx_6_nxp4p0jk = d1_rhs_u_dx_6_nxp4p0jk*param_float(1)

d1_rhs_u_dy_0_nxp4p0jm1k = q(nx+4+0,j-1,indvars(1))*q(nx+4+0,j-1,indvars(2))*q(nx+4+0,j-1,indvars(3))

d1_rhs_u_dy_0_nxp4p0jp1k = q(nx+4+0,j+1,indvars(1))*q(nx+4+0,j+1,indvars(2))*q(nx+4+0,j+1,indvars(3))

d1_rhs_u_dy_0_nxp4p0jk = -&
          0.5_wp*d1_rhs_u_dy_0_nxp4p0jm1k+&
          0.5_wp*d1_rhs_u_dy_0_nxp4p0jp1k

d1_rhs_u_dy_0_nxp4p0jk = d1_rhs_u_dy_0_nxp4p0jk*param_float(2)

d2_rhs_u_dydx_1_0_nxp4p0jm1k_nxp4p0p0jm1k = q(nx+4+0+0,j-1,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p0jm1k_nxp4p0m1jm1k = q(nx+4+0-1,j-1,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p0jm1k_nxp4p0m2jm1k = q(nx+4+0-2,j-1,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p0jm1k = 1.5_wp*d2_rhs_u_dydx_1_0_nxp4p0jm1k_nxp4p0p0jm1k-&
          2.0_wp*d2_rhs_u_dydx_1_0_nxp4p0jm1k_nxp4p0m1jm1k+&
          0.5_wp*d2_rhs_u_dydx_1_0_nxp4p0jm1k_nxp4p0m2jm1k

d2_rhs_u_dydx_1_0_nxp4p0jm1k = d2_rhs_u_dydx_1_0_nxp4p0jm1k*param_float(1)

d2_rhs_u_dydx_1_0_nxp4p0jp1k_nxp4p0p0jp1k = q(nx+4+0+0,j+1,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p0jp1k_nxp4p0m1jp1k = q(nx+4+0-1,j+1,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p0jp1k_nxp4p0m2jp1k = q(nx+4+0-2,j+1,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p0jp1k = 1.5_wp*d2_rhs_u_dydx_1_0_nxp4p0jp1k_nxp4p0p0jp1k-&
          2.0_wp*d2_rhs_u_dydx_1_0_nxp4p0jp1k_nxp4p0m1jp1k+&
          0.5_wp*d2_rhs_u_dydx_1_0_nxp4p0jp1k_nxp4p0m2jp1k

d2_rhs_u_dydx_1_0_nxp4p0jp1k = d2_rhs_u_dydx_1_0_nxp4p0jp1k*param_float(1)

d2_rhs_u_dydy_1_0_nxp4p0jm1k_nxp4p0jm1m1k = q(nx+4+0,j-1-1,indvars(2))

d2_rhs_u_dydy_1_0_nxp4p0jm1k_nxp4p0jm1p1k = q(nx+4+0,j-1+1,indvars(2))

d2_rhs_u_dydy_1_0_nxp4p0jm1k = -&
          0.5_wp*d2_rhs_u_dydy_1_0_nxp4p0jm1k_nxp4p0jm1m1k+&
          0.5_wp*d2_rhs_u_dydy_1_0_nxp4p0jm1k_nxp4p0jm1p1k

d2_rhs_u_dydy_1_0_nxp4p0jm1k = d2_rhs_u_dydy_1_0_nxp4p0jm1k*param_float(2)

d2_rhs_u_dydy_1_0_nxp4p0jp1k_nxp4p0jp1m1k = q(nx+4+0,j+1-1,indvars(2))

d2_rhs_u_dydy_1_0_nxp4p0jp1k_nxp4p0jp1p1k = q(nx+4+0,j+1+1,indvars(2))

d2_rhs_u_dydy_1_0_nxp4p0jp1k = -&
          0.5_wp*d2_rhs_u_dydy_1_0_nxp4p0jp1k_nxp4p0jp1m1k+&
          0.5_wp*d2_rhs_u_dydy_1_0_nxp4p0jp1k_nxp4p0jp1p1k

d2_rhs_u_dydy_1_0_nxp4p0jp1k = d2_rhs_u_dydy_1_0_nxp4p0jp1k*param_float(2)

d1_rhs_u_dy_1_nxp4p0jm1k = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0,j-1,indvars(5))/1.0_wp*q(nx+4+0,j-1,indvars(1)))**3.0_wp/((q(nx+4+0,j-1,indvars(5))/1.0_wp*q(nx+4+0,j-1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0,j-1,indvars(5))/1.0_wp*q(nx+4+0,j-1,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0,j-1,indvarsst(11))*(d2_rhs_u_dydy_1_0_nxp4p0jm1k)+&
                    qst(nx+4+0,j-1,indvarsst(10))*(d2_rhs_u_dydx_1_0_nxp4p0jm1k))

d1_rhs_u_dy_1_nxp4p0jp1k = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0,j+1,indvars(5))/1.0_wp*q(nx+4+0,j+1,indvars(1)))**3.0_wp/((q(nx+4+0,j+1,indvars(5))/1.0_wp*q(nx+4+0,j+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0,j+1,indvars(5))/1.0_wp*q(nx+4+0,j+1,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0,j+1,indvarsst(11))*(d2_rhs_u_dydy_1_0_nxp4p0jp1k)+&
                    qst(nx+4+0,j+1,indvarsst(10))*(d2_rhs_u_dydx_1_0_nxp4p0jp1k))

d1_rhs_u_dy_1_nxp4p0jk = -&
          0.5_wp*d1_rhs_u_dy_1_nxp4p0jm1k+&
          0.5_wp*d1_rhs_u_dy_1_nxp4p0jp1k

d1_rhs_u_dy_1_nxp4p0jk = d1_rhs_u_dy_1_nxp4p0jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None d(rho u)/dt ********
!                                                           
!***********************************************************


rhs(nx+4+0,j,indvars(2)) =   -  ( q(nx+4+0,j,indvars(2))*(1.0_wp/(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5**2*((q(nx+4+0,j,indvars(2))*(-&
                    (param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5**2*d1_rhs_rho_dx_1_nxp4p0jk+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p0jk))*qst(nx+4+0,j,indvarsst(10))+&
                    0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5+&
                    q(nx+4+0,j,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5*q(nx+4+0,j,indvars(1))*d1_rhs_rho_dx_3_nxp4p0jk+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p0jk))*qst(nx+4+0,j,indvarsst(10))+&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5*(1-&
                    qface_i(j,2)*qface_i(j,2))/param_float(20 + 5)*((param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))-&
                    param_float(26 + 5)))))+&
                    q(nx+4+0,j,indvars(1))*(1.0_wp/(2*q(nx+4+0,j,indvars(1))*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5)*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5+&
                    q(nx+4+0,j,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5*q(nx+4+0,j,indvars(1))*d1_rhs_rho_dx_3_nxp4p0jk+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p0jk))*qst(nx+4+0,j,indvarsst(10))-&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5*(1-&
                    qface_i(j,2)*qface_i(j,2))/param_float(20 + 5)*((param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))-&
                    param_float(26 + 5))))+&
                    d1_rhs_u_dy_0_nxp4p0jk+&
                    qst(nx+4+0,j,indvarsst(10))*(d1_rhs_u_dx_6_nxp4p0jk)+&
                    qst(nx+4+0,j,indvarsst(11))*(d1_rhs_u_dy_1_nxp4p0jk) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None d(rho v)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! v*(1.0_wp/c**2*((u*(-c**2*[rho]_1x+1.0_wp*[p]_1x))*deltaxI+0.5_wp*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI+esse*c*(1-M_imax*M_imax)/L_ref*(p-P0))))+rho*((1.0_wp*u*[v]_1x)*deltaxI)+[rho*v*v]_1y+deltaxI*([-visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x))]_1x)+deltayI*([-visc_t*(2.0_wp*deltayI*({v}_1y)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y)))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_rhs_v_dx_4_nxp4p0p0jk = q(nx+4+0+0,j,indvars(3))

d1_rhs_v_dx_4_nxp4p0m1jk = q(nx+4+0-1,j,indvars(3))

d1_rhs_v_dx_4_nxp4p0m2jk = q(nx+4+0-2,j,indvars(3))

d1_rhs_v_dx_4_nxp4p0jk = 1.5_wp*d1_rhs_v_dx_4_nxp4p0p0jk-&
          2.0_wp*d1_rhs_v_dx_4_nxp4p0m1jk+&
          0.5_wp*d1_rhs_v_dx_4_nxp4p0m2jk

d1_rhs_v_dx_4_nxp4p0jk = d1_rhs_v_dx_4_nxp4p0jk*param_float(1)

d2_rhs_v_dxdx_5_0_nxp4p0p0jk_nxp4p0p0p0jk = q(nx+4+0+0+0,j,indvars(3))

d2_rhs_v_dxdx_5_0_nxp4p0p0jk_nxp4p0p0m1jk = q(nx+4+0+0-1,j,indvars(3))

d2_rhs_v_dxdx_5_0_nxp4p0p0jk_nxp4p0p0m2jk = q(nx+4+0+0-2,j,indvars(3))

d2_rhs_v_dxdx_5_0_nxp4p0p0jk = 1.5_wp*d2_rhs_v_dxdx_5_0_nxp4p0p0jk_nxp4p0p0p0jk-&
          2.0_wp*d2_rhs_v_dxdx_5_0_nxp4p0p0jk_nxp4p0p0m1jk+&
          0.5_wp*d2_rhs_v_dxdx_5_0_nxp4p0p0jk_nxp4p0p0m2jk

d2_rhs_v_dxdx_5_0_nxp4p0p0jk = d2_rhs_v_dxdx_5_0_nxp4p0p0jk*param_float(1)

d2_rhs_v_dxdx_5_0_nxp4p0m1jk_nxp4p0m1m1jk = q(nx+4+0-1-1,j,indvars(3))

d2_rhs_v_dxdx_5_0_nxp4p0m1jk_nxp4p0m1p1jk = q(nx+4+0-1+1,j,indvars(3))

d2_rhs_v_dxdx_5_0_nxp4p0m1jk = -&
          0.5_wp*d2_rhs_v_dxdx_5_0_nxp4p0m1jk_nxp4p0m1m1jk+&
          0.5_wp*d2_rhs_v_dxdx_5_0_nxp4p0m1jk_nxp4p0m1p1jk

d2_rhs_v_dxdx_5_0_nxp4p0m1jk = d2_rhs_v_dxdx_5_0_nxp4p0m1jk*param_float(1)

d2_rhs_v_dxdx_5_0_nxp4p0m2jk_nxp4p0m2m1jk = q(nx+4+0-2-1,j,indvars(3))

d2_rhs_v_dxdx_5_0_nxp4p0m2jk_nxp4p0m2p1jk = q(nx+4+0-2+1,j,indvars(3))

d2_rhs_v_dxdx_5_0_nxp4p0m2jk = -&
          0.5_wp*d2_rhs_v_dxdx_5_0_nxp4p0m2jk_nxp4p0m2m1jk+&
          0.5_wp*d2_rhs_v_dxdx_5_0_nxp4p0m2jk_nxp4p0m2p1jk

d2_rhs_v_dxdx_5_0_nxp4p0m2jk = d2_rhs_v_dxdx_5_0_nxp4p0m2jk*param_float(1)

d2_rhs_v_dxdy_5_0_nxp4p0p0jk_nxp4p0p0jm1k = q(nx+4+0+0,j-1,indvars(2))

d2_rhs_v_dxdy_5_0_nxp4p0p0jk_nxp4p0p0jp1k = q(nx+4+0+0,j+1,indvars(2))

d2_rhs_v_dxdy_5_0_nxp4p0p0jk = -&
          0.5_wp*d2_rhs_v_dxdy_5_0_nxp4p0p0jk_nxp4p0p0jm1k+&
          0.5_wp*d2_rhs_v_dxdy_5_0_nxp4p0p0jk_nxp4p0p0jp1k

d2_rhs_v_dxdy_5_0_nxp4p0p0jk = d2_rhs_v_dxdy_5_0_nxp4p0p0jk*param_float(2)

d2_rhs_v_dxdy_5_0_nxp4p0m1jk_nxp4p0m1jm1k = q(nx+4+0-1,j-1,indvars(2))

d2_rhs_v_dxdy_5_0_nxp4p0m1jk_nxp4p0m1jp1k = q(nx+4+0-1,j+1,indvars(2))

d2_rhs_v_dxdy_5_0_nxp4p0m1jk = -&
          0.5_wp*d2_rhs_v_dxdy_5_0_nxp4p0m1jk_nxp4p0m1jm1k+&
          0.5_wp*d2_rhs_v_dxdy_5_0_nxp4p0m1jk_nxp4p0m1jp1k

d2_rhs_v_dxdy_5_0_nxp4p0m1jk = d2_rhs_v_dxdy_5_0_nxp4p0m1jk*param_float(2)

d2_rhs_v_dxdy_5_0_nxp4p0m2jk_nxp4p0m2jm1k = q(nx+4+0-2,j-1,indvars(2))

d2_rhs_v_dxdy_5_0_nxp4p0m2jk_nxp4p0m2jp1k = q(nx+4+0-2,j+1,indvars(2))

d2_rhs_v_dxdy_5_0_nxp4p0m2jk = -&
          0.5_wp*d2_rhs_v_dxdy_5_0_nxp4p0m2jk_nxp4p0m2jm1k+&
          0.5_wp*d2_rhs_v_dxdy_5_0_nxp4p0m2jk_nxp4p0m2jp1k

d2_rhs_v_dxdy_5_0_nxp4p0m2jk = d2_rhs_v_dxdy_5_0_nxp4p0m2jk*param_float(2)

d1_rhs_v_dx_5_nxp4p0p0jk = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0+0,j,indvars(5))/1.0_wp*q(nx+4+0+0,j,indvars(1)))**3.0_wp/((q(nx+4+0+0,j,indvars(5))/1.0_wp*q(nx+4+0+0,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0+0,j,indvars(5))/1.0_wp*q(nx+4+0+0,j,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0+0,j,indvarsst(11))*(d2_rhs_v_dxdy_5_0_nxp4p0p0jk)+&
                    qst(nx+4+0+0,j,indvarsst(10))*(d2_rhs_v_dxdx_5_0_nxp4p0p0jk))

d1_rhs_v_dx_5_nxp4p0m1jk = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0-1,j,indvars(5))/1.0_wp*q(nx+4+0-1,j,indvars(1)))**3.0_wp/((q(nx+4+0-1,j,indvars(5))/1.0_wp*q(nx+4+0-1,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0-1,j,indvars(5))/1.0_wp*q(nx+4+0-1,j,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0-1,j,indvarsst(11))*(d2_rhs_v_dxdy_5_0_nxp4p0m1jk)+&
                    qst(nx+4+0-1,j,indvarsst(10))*(d2_rhs_v_dxdx_5_0_nxp4p0m1jk))

d1_rhs_v_dx_5_nxp4p0m2jk = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0-2,j,indvars(5))/1.0_wp*q(nx+4+0-2,j,indvars(1)))**3.0_wp/((q(nx+4+0-2,j,indvars(5))/1.0_wp*q(nx+4+0-2,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0-2,j,indvars(5))/1.0_wp*q(nx+4+0-2,j,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0-2,j,indvarsst(11))*(d2_rhs_v_dxdy_5_0_nxp4p0m2jk)+&
                    qst(nx+4+0-2,j,indvarsst(10))*(d2_rhs_v_dxdx_5_0_nxp4p0m2jk))

d1_rhs_v_dx_5_nxp4p0jk = 1.5_wp*d1_rhs_v_dx_5_nxp4p0p0jk-&
          2.0_wp*d1_rhs_v_dx_5_nxp4p0m1jk+&
          0.5_wp*d1_rhs_v_dx_5_nxp4p0m2jk

d1_rhs_v_dx_5_nxp4p0jk = d1_rhs_v_dx_5_nxp4p0jk*param_float(1)

d1_rhs_v_dy_0_nxp4p0jm1k = q(nx+4+0,j-1,indvars(1))*q(nx+4+0,j-1,indvars(3))*q(nx+4+0,j-1,indvars(3))

d1_rhs_v_dy_0_nxp4p0jp1k = q(nx+4+0,j+1,indvars(1))*q(nx+4+0,j+1,indvars(3))*q(nx+4+0,j+1,indvars(3))

d1_rhs_v_dy_0_nxp4p0jk = -&
          0.5_wp*d1_rhs_v_dy_0_nxp4p0jm1k+&
          0.5_wp*d1_rhs_v_dy_0_nxp4p0jp1k

d1_rhs_v_dy_0_nxp4p0jk = d1_rhs_v_dy_0_nxp4p0jk*param_float(2)

d2_rhs_v_dydx_1_0_nxp4p0jm1k_nxp4p0p0jm1k = q(nx+4+0+0,j-1,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p0jm1k_nxp4p0m1jm1k = q(nx+4+0-1,j-1,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p0jm1k_nxp4p0m2jm1k = q(nx+4+0-2,j-1,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p0jm1k = 1.5_wp*d2_rhs_v_dydx_1_0_nxp4p0jm1k_nxp4p0p0jm1k-&
          2.0_wp*d2_rhs_v_dydx_1_0_nxp4p0jm1k_nxp4p0m1jm1k+&
          0.5_wp*d2_rhs_v_dydx_1_0_nxp4p0jm1k_nxp4p0m2jm1k

d2_rhs_v_dydx_1_0_nxp4p0jm1k = d2_rhs_v_dydx_1_0_nxp4p0jm1k*param_float(1)

d2_rhs_v_dydx_1_0_nxp4p0jp1k_nxp4p0p0jp1k = q(nx+4+0+0,j+1,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p0jp1k_nxp4p0m1jp1k = q(nx+4+0-1,j+1,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p0jp1k_nxp4p0m2jp1k = q(nx+4+0-2,j+1,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p0jp1k = 1.5_wp*d2_rhs_v_dydx_1_0_nxp4p0jp1k_nxp4p0p0jp1k-&
          2.0_wp*d2_rhs_v_dydx_1_0_nxp4p0jp1k_nxp4p0m1jp1k+&
          0.5_wp*d2_rhs_v_dydx_1_0_nxp4p0jp1k_nxp4p0m2jp1k

d2_rhs_v_dydx_1_0_nxp4p0jp1k = d2_rhs_v_dydx_1_0_nxp4p0jp1k*param_float(1)

d2_rhs_v_dydy_1_0_nxp4p0jm1k_nxp4p0jm1m1k = q(nx+4+0,j-1-1,indvars(3))

d2_rhs_v_dydy_1_0_nxp4p0jm1k_nxp4p0jm1p1k = q(nx+4+0,j-1+1,indvars(3))

d2_rhs_v_dydy_1_0_nxp4p0jm1k = -&
          0.5_wp*d2_rhs_v_dydy_1_0_nxp4p0jm1k_nxp4p0jm1m1k+&
          0.5_wp*d2_rhs_v_dydy_1_0_nxp4p0jm1k_nxp4p0jm1p1k

d2_rhs_v_dydy_1_0_nxp4p0jm1k = d2_rhs_v_dydy_1_0_nxp4p0jm1k*param_float(2)

d2_rhs_v_dydy_1_0_nxp4p0jp1k_nxp4p0jp1m1k = q(nx+4+0,j+1-1,indvars(3))

d2_rhs_v_dydy_1_0_nxp4p0jp1k_nxp4p0jp1p1k = q(nx+4+0,j+1+1,indvars(3))

d2_rhs_v_dydy_1_0_nxp4p0jp1k = -&
          0.5_wp*d2_rhs_v_dydy_1_0_nxp4p0jp1k_nxp4p0jp1m1k+&
          0.5_wp*d2_rhs_v_dydy_1_0_nxp4p0jp1k_nxp4p0jp1p1k

d2_rhs_v_dydy_1_0_nxp4p0jp1k = d2_rhs_v_dydy_1_0_nxp4p0jp1k*param_float(2)

d1_rhs_v_dy_1_nxp4p0jm1k = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0,j-1,indvars(5))/1.0_wp*q(nx+4+0,j-1,indvars(1)))**3.0_wp/((q(nx+4+0,j-1,indvars(5))/1.0_wp*q(nx+4+0,j-1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0,j-1,indvars(5))/1.0_wp*q(nx+4+0,j-1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0,j-1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp4p0jm1k)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0,j-1,indvarsst(10))*(d2_rhs_v_dydx_1_0_nxp4p0jm1k)+&
                    qst(nx+4+0,j-1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp4p0jm1k)))

d1_rhs_v_dy_1_nxp4p0jp1k = -(1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0,j+1,indvars(5))/1.0_wp*q(nx+4+0,j+1,indvars(1)))**3.0_wp/((q(nx+4+0,j+1,indvars(5))/1.0_wp*q(nx+4+0,j+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0,j+1,indvars(5))/1.0_wp*q(nx+4+0,j+1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0,j+1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp4p0jp1k)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0,j+1,indvarsst(10))*(d2_rhs_v_dydx_1_0_nxp4p0jp1k)+&
                    qst(nx+4+0,j+1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp4p0jp1k)))

d1_rhs_v_dy_1_nxp4p0jk = -&
          0.5_wp*d1_rhs_v_dy_1_nxp4p0jm1k+&
          0.5_wp*d1_rhs_v_dy_1_nxp4p0jp1k

d1_rhs_v_dy_1_nxp4p0jk = d1_rhs_v_dy_1_nxp4p0jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None d(rho v)/dt ********
!                                                           
!***********************************************************


rhs(nx+4+0,j,indvars(3)) =   -  ( q(nx+4+0,j,indvars(3))*(1.0_wp/(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5**2*((q(nx+4+0,j,indvars(2))*(-&
                    (param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5**2*d1_rhs_rho_dx_1_nxp4p0jk+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p0jk))*qst(nx+4+0,j,indvarsst(10))+&
                    0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5+&
                    q(nx+4+0,j,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5*q(nx+4+0,j,indvars(1))*d1_rhs_rho_dx_3_nxp4p0jk+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p0jk))*qst(nx+4+0,j,indvarsst(10))+&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5*(1-&
                    qface_i(j,2)*qface_i(j,2))/param_float(20 + 5)*((param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))-&
                    param_float(26 + 5)))))+&
                    q(nx+4+0,j,indvars(1))*((1.0_wp*q(nx+4+0,j,indvars(2))*d1_rhs_v_dx_4_nxp4p0jk)*qst(nx+4+0,j,indvarsst(10)))+&
                    d1_rhs_v_dy_0_nxp4p0jk+&
                    qst(nx+4+0,j,indvarsst(10))*(d1_rhs_v_dx_5_nxp4p0jk)+&
                    qst(nx+4+0,j,indvarsst(11))*(d1_rhs_v_dy_1_nxp4p0jk) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None d(rho et)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! 0.5_wp*(u**2+v**2)*(1.0_wp/c**2*((u*(-c**2*[rho]_1x+1.0_wp*[p]_1x))*deltaxI+0.5_wp*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI+esse*c*(1-M_imax*M_imax)/L_ref*(p-P0))))+1.0_wp/gamma_m1*(0.5_wp*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI+esse*c*(1-M_imax*M_imax)/L_ref*(p-P0)))+rho*u*(1.0_wp/(2*rho*c)*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI-esse*c*(1-M_imax*M_imax)/L_ref*(p-P0)))+rho*v*((1.0_wp*u*[v]_1x)*deltaxI)+[(rho*et+p)*v]_1y+deltaxI*([-kappa*deltaxI*({T}_1x)-u*(visc_t*(2.0_wp*deltaxI*({u}_1x)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y))))-v*(visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x)))]_1x)+deltayI*([-kappa*deltayI*({T}_1y)-u*(visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x)))-v*(visc_t*(2.0_wp*deltayI*({v}_1y)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y))))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_rhs_et_dxdx_9_0_nxp4p0p0jk_nxp4p0p0p0jk = ((q(nx+4+0+0+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0+0,j,indvars(2))*q(nx+4+0+0+0,j,indvars(2))+&
                    q(nx+4+0+0+0,j,indvars(3))*q(nx+4+0+0+0,j,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dxdx_9_0_nxp4p0p0jk_nxp4p0p0m1jk = ((q(nx+4+0+0-1,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0-1,j,indvars(2))*q(nx+4+0+0-1,j,indvars(2))+&
                    q(nx+4+0+0-1,j,indvars(3))*q(nx+4+0+0-1,j,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dxdx_9_0_nxp4p0p0jk_nxp4p0p0m2jk = ((q(nx+4+0+0-2,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0-2,j,indvars(2))*q(nx+4+0+0-2,j,indvars(2))+&
                    q(nx+4+0+0-2,j,indvars(3))*q(nx+4+0+0-2,j,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dxdx_9_0_nxp4p0p0jk = 1.5_wp*d2_rhs_et_dxdx_9_0_nxp4p0p0jk_nxp4p0p0p0jk-&
          2.0_wp*d2_rhs_et_dxdx_9_0_nxp4p0p0jk_nxp4p0p0m1jk+&
          0.5_wp*d2_rhs_et_dxdx_9_0_nxp4p0p0jk_nxp4p0p0m2jk

d2_rhs_et_dxdx_9_0_nxp4p0p0jk = d2_rhs_et_dxdx_9_0_nxp4p0p0jk*param_float(1)

d2_rhs_et_dxdx_9_0_nxp4p0m1jk_nxp4p0m1m1jk = ((q(nx+4+0-1-1,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1-1,j,indvars(2))*q(nx+4+0-1-1,j,indvars(2))+&
                    q(nx+4+0-1-1,j,indvars(3))*q(nx+4+0-1-1,j,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dxdx_9_0_nxp4p0m1jk_nxp4p0m1p1jk = ((q(nx+4+0-1+1,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1+1,j,indvars(2))*q(nx+4+0-1+1,j,indvars(2))+&
                    q(nx+4+0-1+1,j,indvars(3))*q(nx+4+0-1+1,j,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dxdx_9_0_nxp4p0m1jk = -&
          0.5_wp*d2_rhs_et_dxdx_9_0_nxp4p0m1jk_nxp4p0m1m1jk+&
          0.5_wp*d2_rhs_et_dxdx_9_0_nxp4p0m1jk_nxp4p0m1p1jk

d2_rhs_et_dxdx_9_0_nxp4p0m1jk = d2_rhs_et_dxdx_9_0_nxp4p0m1jk*param_float(1)

d2_rhs_et_dxdx_9_0_nxp4p0m2jk_nxp4p0m2m1jk = ((q(nx+4+0-2-1,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2-1,j,indvars(2))*q(nx+4+0-2-1,j,indvars(2))+&
                    q(nx+4+0-2-1,j,indvars(3))*q(nx+4+0-2-1,j,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dxdx_9_0_nxp4p0m2jk_nxp4p0m2p1jk = ((q(nx+4+0-2+1,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2+1,j,indvars(2))*q(nx+4+0-2+1,j,indvars(2))+&
                    q(nx+4+0-2+1,j,indvars(3))*q(nx+4+0-2+1,j,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dxdx_9_0_nxp4p0m2jk = -&
          0.5_wp*d2_rhs_et_dxdx_9_0_nxp4p0m2jk_nxp4p0m2m1jk+&
          0.5_wp*d2_rhs_et_dxdx_9_0_nxp4p0m2jk_nxp4p0m2p1jk

d2_rhs_et_dxdx_9_0_nxp4p0m2jk = d2_rhs_et_dxdx_9_0_nxp4p0m2jk*param_float(1)

d1_rhs_et_dx_9_nxp4p0p0jk = -param_float(2 + 5)*qst(nx+4+0+0,j,indvarsst(10))*(d2_rhs_et_dxdx_9_0_nxp4p0p0jk)-&
                    q(nx+4+0+0,j,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0+0,j,indvars(5))/1.0_wp*q(nx+4+0+0,j,indvars(1)))**3.0_wp/((q(nx+4+0+0,j,indvars(5))/1.0_wp*q(nx+4+0+0,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0+0,j,indvars(5))/1.0_wp*q(nx+4+0+0,j,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0+0,j,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0p0jk)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0+0,j,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0p0jk)+&
                    qst(nx+4+0+0,j,indvarsst(11))*(d2_rhs_u_dxdy_6_0_nxp4p0p0jk))))-&
                    q(nx+4+0+0,j,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0+0,j,indvars(5))/1.0_wp*q(nx+4+0+0,j,indvars(1)))**3.0_wp/((q(nx+4+0+0,j,indvars(5))/1.0_wp*q(nx+4+0+0,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0+0,j,indvars(5))/1.0_wp*q(nx+4+0+0,j,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0+0,j,indvarsst(11))*(d2_rhs_v_dxdy_5_0_nxp4p0p0jk)+&
                    qst(nx+4+0+0,j,indvarsst(10))*(d2_rhs_v_dxdx_5_0_nxp4p0p0jk)))

d1_rhs_et_dx_9_nxp4p0m1jk = -param_float(2 + 5)*qst(nx+4+0-1,j,indvarsst(10))*(d2_rhs_et_dxdx_9_0_nxp4p0m1jk)-&
                    q(nx+4+0-1,j,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0-1,j,indvars(5))/1.0_wp*q(nx+4+0-1,j,indvars(1)))**3.0_wp/((q(nx+4+0-1,j,indvars(5))/1.0_wp*q(nx+4+0-1,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0-1,j,indvars(5))/1.0_wp*q(nx+4+0-1,j,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0-1,j,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0m1jk)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0-1,j,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0m1jk)+&
                    qst(nx+4+0-1,j,indvarsst(11))*(d2_rhs_u_dxdy_6_0_nxp4p0m1jk))))-&
                    q(nx+4+0-1,j,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0-1,j,indvars(5))/1.0_wp*q(nx+4+0-1,j,indvars(1)))**3.0_wp/((q(nx+4+0-1,j,indvars(5))/1.0_wp*q(nx+4+0-1,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0-1,j,indvars(5))/1.0_wp*q(nx+4+0-1,j,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0-1,j,indvarsst(11))*(d2_rhs_v_dxdy_5_0_nxp4p0m1jk)+&
                    qst(nx+4+0-1,j,indvarsst(10))*(d2_rhs_v_dxdx_5_0_nxp4p0m1jk)))

d1_rhs_et_dx_9_nxp4p0m2jk = -param_float(2 + 5)*qst(nx+4+0-2,j,indvarsst(10))*(d2_rhs_et_dxdx_9_0_nxp4p0m2jk)-&
                    q(nx+4+0-2,j,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0-2,j,indvars(5))/1.0_wp*q(nx+4+0-2,j,indvars(1)))**3.0_wp/((q(nx+4+0-2,j,indvars(5))/1.0_wp*q(nx+4+0-2,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0-2,j,indvars(5))/1.0_wp*q(nx+4+0-2,j,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0-2,j,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0m2jk)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0-2,j,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0m2jk)+&
                    qst(nx+4+0-2,j,indvarsst(11))*(d2_rhs_u_dxdy_6_0_nxp4p0m2jk))))-&
                    q(nx+4+0-2,j,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0-2,j,indvars(5))/1.0_wp*q(nx+4+0-2,j,indvars(1)))**3.0_wp/((q(nx+4+0-2,j,indvars(5))/1.0_wp*q(nx+4+0-2,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0-2,j,indvars(5))/1.0_wp*q(nx+4+0-2,j,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0-2,j,indvarsst(11))*(d2_rhs_v_dxdy_5_0_nxp4p0m2jk)+&
                    qst(nx+4+0-2,j,indvarsst(10))*(d2_rhs_v_dxdx_5_0_nxp4p0m2jk)))

d1_rhs_et_dx_9_nxp4p0jk = 1.5_wp*d1_rhs_et_dx_9_nxp4p0p0jk-&
          2.0_wp*d1_rhs_et_dx_9_nxp4p0m1jk+&
          0.5_wp*d1_rhs_et_dx_9_nxp4p0m2jk

d1_rhs_et_dx_9_nxp4p0jk = d1_rhs_et_dx_9_nxp4p0jk*param_float(1)

d1_rhs_et_dy_0_nxp4p0jm1k = (q(nx+4+0,j-1,indvars(1))*q(nx+4+0,j-1,indvars(4))+&
                    (param_float(3 + 5))*q(nx+4+0,j-1,indvars(1))*((q(nx+4+0,j-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j-1,indvars(2))*q(nx+4+0,j-1,indvars(2))+&
                    q(nx+4+0,j-1,indvars(3))*q(nx+4+0,j-1,indvars(3))))))*q(nx+4+0,j-1,indvars(3))

d1_rhs_et_dy_0_nxp4p0jp1k = (q(nx+4+0,j+1,indvars(1))*q(nx+4+0,j+1,indvars(4))+&
                    (param_float(3 + 5))*q(nx+4+0,j+1,indvars(1))*((q(nx+4+0,j+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j+1,indvars(2))*q(nx+4+0,j+1,indvars(2))+&
                    q(nx+4+0,j+1,indvars(3))*q(nx+4+0,j+1,indvars(3))))))*q(nx+4+0,j+1,indvars(3))

d1_rhs_et_dy_0_nxp4p0jk = -&
          0.5_wp*d1_rhs_et_dy_0_nxp4p0jm1k+&
          0.5_wp*d1_rhs_et_dy_0_nxp4p0jp1k

d1_rhs_et_dy_0_nxp4p0jk = d1_rhs_et_dy_0_nxp4p0jk*param_float(2)

d2_rhs_et_dydy_1_0_nxp4p0jm1k_nxp4p0jm1m1k = ((q(nx+4+0,j-1-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j-1-1,indvars(2))*q(nx+4+0,j-1-1,indvars(2))+&
                    q(nx+4+0,j-1-1,indvars(3))*q(nx+4+0,j-1-1,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_1_0_nxp4p0jm1k_nxp4p0jm1p1k = ((q(nx+4+0,j-1+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j-1+1,indvars(2))*q(nx+4+0,j-1+1,indvars(2))+&
                    q(nx+4+0,j-1+1,indvars(3))*q(nx+4+0,j-1+1,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_1_0_nxp4p0jm1k = -&
          0.5_wp*d2_rhs_et_dydy_1_0_nxp4p0jm1k_nxp4p0jm1m1k+&
          0.5_wp*d2_rhs_et_dydy_1_0_nxp4p0jm1k_nxp4p0jm1p1k

d2_rhs_et_dydy_1_0_nxp4p0jm1k = d2_rhs_et_dydy_1_0_nxp4p0jm1k*param_float(2)

d2_rhs_et_dydy_1_0_nxp4p0jp1k_nxp4p0jp1m1k = ((q(nx+4+0,j+1-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j+1-1,indvars(2))*q(nx+4+0,j+1-1,indvars(2))+&
                    q(nx+4+0,j+1-1,indvars(3))*q(nx+4+0,j+1-1,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_1_0_nxp4p0jp1k_nxp4p0jp1p1k = ((q(nx+4+0,j+1+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j+1+1,indvars(2))*q(nx+4+0,j+1+1,indvars(2))+&
                    q(nx+4+0,j+1+1,indvars(3))*q(nx+4+0,j+1+1,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_1_0_nxp4p0jp1k = -&
          0.5_wp*d2_rhs_et_dydy_1_0_nxp4p0jp1k_nxp4p0jp1m1k+&
          0.5_wp*d2_rhs_et_dydy_1_0_nxp4p0jp1k_nxp4p0jp1p1k

d2_rhs_et_dydy_1_0_nxp4p0jp1k = d2_rhs_et_dydy_1_0_nxp4p0jp1k*param_float(2)

d1_rhs_et_dy_1_nxp4p0jm1k = -param_float(2 + 5)*qst(nx+4+0,j-1,indvarsst(11))*(d2_rhs_et_dydy_1_0_nxp4p0jm1k)-&
                    q(nx+4+0,j-1,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0,j-1,indvars(5))/1.0_wp*q(nx+4+0,j-1,indvars(1)))**3.0_wp/((q(nx+4+0,j-1,indvars(5))/1.0_wp*q(nx+4+0,j-1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0,j-1,indvars(5))/1.0_wp*q(nx+4+0,j-1,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0,j-1,indvarsst(11))*(d2_rhs_u_dydy_1_0_nxp4p0jm1k)+&
                    qst(nx+4+0,j-1,indvarsst(10))*(d2_rhs_u_dydx_1_0_nxp4p0jm1k)))-&
                    q(nx+4+0,j-1,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0,j-1,indvars(5))/1.0_wp*q(nx+4+0,j-1,indvars(1)))**3.0_wp/((q(nx+4+0,j-1,indvars(5))/1.0_wp*q(nx+4+0,j-1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0,j-1,indvars(5))/1.0_wp*q(nx+4+0,j-1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0,j-1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp4p0jm1k)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0,j-1,indvarsst(10))*(d2_rhs_v_dydx_1_0_nxp4p0jm1k)+&
                    qst(nx+4+0,j-1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp4p0jm1k))))

d1_rhs_et_dy_1_nxp4p0jp1k = -param_float(2 + 5)*qst(nx+4+0,j+1,indvarsst(11))*(d2_rhs_et_dydy_1_0_nxp4p0jp1k)-&
                    q(nx+4+0,j+1,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0,j+1,indvars(5))/1.0_wp*q(nx+4+0,j+1,indvars(1)))**3.0_wp/((q(nx+4+0,j+1,indvars(5))/1.0_wp*q(nx+4+0,j+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0,j+1,indvars(5))/1.0_wp*q(nx+4+0,j+1,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0,j+1,indvarsst(11))*(d2_rhs_u_dydy_1_0_nxp4p0jp1k)+&
                    qst(nx+4+0,j+1,indvarsst(10))*(d2_rhs_u_dydx_1_0_nxp4p0jp1k)))-&
                    q(nx+4+0,j+1,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0,j+1,indvars(5))/1.0_wp*q(nx+4+0,j+1,indvars(1)))**3.0_wp/((q(nx+4+0,j+1,indvars(5))/1.0_wp*q(nx+4+0,j+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0,j+1,indvars(5))/1.0_wp*q(nx+4+0,j+1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0,j+1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp4p0jp1k)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0,j+1,indvarsst(10))*(d2_rhs_v_dydx_1_0_nxp4p0jp1k)+&
                    qst(nx+4+0,j+1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp4p0jp1k))))

d1_rhs_et_dy_1_nxp4p0jk = -&
          0.5_wp*d1_rhs_et_dy_1_nxp4p0jm1k+&
          0.5_wp*d1_rhs_et_dy_1_nxp4p0jp1k

d1_rhs_et_dy_1_nxp4p0jk = d1_rhs_et_dy_1_nxp4p0jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None d(rho et)/dt *******
!                                                           
!***********************************************************


rhs(nx+4+0,j,indvars(4)) =   -  ( 0.5_wp*(q(nx+4+0,j,indvars(2))**2+&
                    q(nx+4+0,j,indvars(3))**2)*(1.0_wp/(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5**2*((q(nx+4+0,j,indvars(2))*(-&
                    (param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5**2*d1_rhs_rho_dx_1_nxp4p0jk+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p0jk))*qst(nx+4+0,j,indvarsst(10))+&
                    0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5+&
                    q(nx+4+0,j,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5*q(nx+4+0,j,indvars(1))*d1_rhs_rho_dx_3_nxp4p0jk+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p0jk))*qst(nx+4+0,j,indvarsst(10))+&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5*(1-&
                    qface_i(j,2)*qface_i(j,2))/param_float(20 + 5)*((param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))-&
                    param_float(26 + 5)))))+&
                    1.0_wp/param_float(3 + 5)*(0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5+&
                    q(nx+4+0,j,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5*q(nx+4+0,j,indvars(1))*d1_rhs_rho_dx_3_nxp4p0jk+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p0jk))*qst(nx+4+0,j,indvarsst(10))+&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5*(1-&
                    qface_i(j,2)*qface_i(j,2))/param_float(20 + 5)*((param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))-&
                    param_float(26 + 5))))+&
                    q(nx+4+0,j,indvars(1))*q(nx+4+0,j,indvars(2))*(1.0_wp/(2*q(nx+4+0,j,indvars(1))*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5)*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5+&
                    q(nx+4+0,j,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5*q(nx+4+0,j,indvars(1))*d1_rhs_rho_dx_3_nxp4p0jk+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p0jk))*qst(nx+4+0,j,indvarsst(10))-&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/q(nx+4+0,j,indvars(1)))**0.5*(1-&
                    qface_i(j,2)*qface_i(j,2))/param_float(20 + 5)*((param_float(3 + 5))*q(nx+4+0,j,indvars(1))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))-&
                    param_float(26 + 5))))+&
                    q(nx+4+0,j,indvars(1))*q(nx+4+0,j,indvars(3))*((1.0_wp*q(nx+4+0,j,indvars(2))*d1_rhs_v_dx_4_nxp4p0jk)*qst(nx+4+0,j,indvarsst(10)))+&
                    d1_rhs_et_dy_0_nxp4p0jk+&
                    qst(nx+4+0,j,indvarsst(10))*(d1_rhs_et_dx_9_nxp4p0jk)+&
                    qst(nx+4+0,j,indvarsst(11))*(d1_rhs_et_dy_1_nxp4p0jk) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None d(rho nut)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ([rho*u*nut-visc_t*sigmaI*deltayI*({nut}_1y)]_1y)*deltayI+-ReI*Cb2*sigmaI*((deltayI)**2*([rho*nut]_1y)*([nut]_1y))-Cb1*(1-ft2)*SS*rho*nut+ReI*(Cw1*fw-Cb1/k**2*ft2)*rho*nut**2/eta**2
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_rhs_nut_dydy_0_0_nxp4p0jm1k_nxp4p0jm1m1k = q(nx+4+0,j-1-1,indvars(5))

d2_rhs_nut_dydy_0_0_nxp4p0jm1k_nxp4p0jm1p1k = q(nx+4+0,j-1+1,indvars(5))

d2_rhs_nut_dydy_0_0_nxp4p0jm1k = -&
          0.5_wp*d2_rhs_nut_dydy_0_0_nxp4p0jm1k_nxp4p0jm1m1k+&
          0.5_wp*d2_rhs_nut_dydy_0_0_nxp4p0jm1k_nxp4p0jm1p1k

d2_rhs_nut_dydy_0_0_nxp4p0jm1k = d2_rhs_nut_dydy_0_0_nxp4p0jm1k*param_float(2)

d2_rhs_nut_dydy_0_0_nxp4p0jp1k_nxp4p0jp1m1k = q(nx+4+0,j+1-1,indvars(5))

d2_rhs_nut_dydy_0_0_nxp4p0jp1k_nxp4p0jp1p1k = q(nx+4+0,j+1+1,indvars(5))

d2_rhs_nut_dydy_0_0_nxp4p0jp1k = -&
          0.5_wp*d2_rhs_nut_dydy_0_0_nxp4p0jp1k_nxp4p0jp1m1k+&
          0.5_wp*d2_rhs_nut_dydy_0_0_nxp4p0jp1k_nxp4p0jp1p1k

d2_rhs_nut_dydy_0_0_nxp4p0jp1k = d2_rhs_nut_dydy_0_0_nxp4p0jp1k*param_float(2)

d1_rhs_nut_dy_1_nxp4p0jm1k = q(nx+4+0,j-1,indvars(1))*q(nx+4+0,j-1,indvars(5))

d1_rhs_nut_dy_1_nxp4p0jp1k = q(nx+4+0,j+1,indvars(1))*q(nx+4+0,j+1,indvars(5))

d1_rhs_nut_dy_1_nxp4p0jk = -&
          0.5_wp*d1_rhs_nut_dy_1_nxp4p0jm1k+&
          0.5_wp*d1_rhs_nut_dy_1_nxp4p0jp1k

d1_rhs_nut_dy_1_nxp4p0jk = d1_rhs_nut_dy_1_nxp4p0jk*param_float(2)

d1_rhs_nut_dy_2_nxp4p0jm1k = q(nx+4+0,j-1,indvars(5))

d1_rhs_nut_dy_2_nxp4p0jp1k = q(nx+4+0,j+1,indvars(5))

d1_rhs_nut_dy_2_nxp4p0jk = -&
          0.5_wp*d1_rhs_nut_dy_2_nxp4p0jm1k+&
          0.5_wp*d1_rhs_nut_dy_2_nxp4p0jp1k

d1_rhs_nut_dy_2_nxp4p0jk = d1_rhs_nut_dy_2_nxp4p0jk*param_float(2)

d1_rhs_nut_dy_0_nxp4p0jm1k = q(nx+4+0,j-1,indvars(1))*q(nx+4+0,j-1,indvars(2))*q(nx+4+0,j-1,indvars(5))-&
                    (1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0,j-1,indvars(5))/1.0_wp*q(nx+4+0,j-1,indvars(1)))**3.0_wp/((q(nx+4+0,j-1,indvars(5))/1.0_wp*q(nx+4+0,j-1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0,j-1,indvars(5))/1.0_wp*q(nx+4+0,j-1,indvars(1))))*param_float(1 + 5)*param_float(18 + 5)*qst(nx+4+0,j-1,indvarsst(11))*(d2_rhs_nut_dydy_0_0_nxp4p0jm1k)

d1_rhs_nut_dy_0_nxp4p0jp1k = q(nx+4+0,j+1,indvars(1))*q(nx+4+0,j+1,indvars(2))*q(nx+4+0,j+1,indvars(5))-&
                    (1.0_wp)*(1.0_wp+&
                    ((q(nx+4+0,j+1,indvars(5))/1.0_wp*q(nx+4+0,j+1,indvars(1)))**3.0_wp/((q(nx+4+0,j+1,indvars(5))/1.0_wp*q(nx+4+0,j+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4+0,j+1,indvars(5))/1.0_wp*q(nx+4+0,j+1,indvars(1))))*param_float(1 + 5)*param_float(18 + 5)*qst(nx+4+0,j+1,indvarsst(11))*(d2_rhs_nut_dydy_0_0_nxp4p0jp1k)

d1_rhs_nut_dy_0_nxp4p0jk = -&
          0.5_wp*d1_rhs_nut_dy_0_nxp4p0jm1k+&
          0.5_wp*d1_rhs_nut_dy_0_nxp4p0jp1k

d1_rhs_nut_dy_0_nxp4p0jk = d1_rhs_nut_dy_0_nxp4p0jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None d(rho nut)/dt ******
!                                                           
!***********************************************************


rhs(nx+4+0,j,indvars(5)) =   -  ( (d1_rhs_nut_dy_0_nxp4p0jk)*qst(nx+4+0,j,indvarsst(11))+&
                    -&
                    param_float(1 + 5)*param_float(7 + 5)*param_float(18 + 5)*((qst(nx+4+0,j,indvarsst(11)))**2*(d1_rhs_nut_dy_1_nxp4p0jk)*(d1_rhs_nut_dy_2_nxp4p0jk))-&
                    param_float(6 + 5)*(1-&
                    (param_float(16 + 5)*exp(-&
                    param_float(17 + 5)*(q(nx+4+0,j,indvars(5))/1.0_wp*q(nx+4+0,j,indvars(1)))**2.0_wp)))*qst(nx+4+0,j,indvarsst(12))*q(nx+4+0,j,indvars(1))*q(nx+4+0,j,indvars(5))+&
                    param_float(1 + 5)*(param_float(10 + 5)*(((min(param_float(1 + 5)*(q(nx+4+0,j,indvars(5))/(qst(nx+4+0,j,indvarsst(12))*param_float(9 + 5)**2.0_wp*qst(nx+4+0,j,indvarsst(2))**2.0_wp)),10.0_wp))+&
                    param_float(11 + 5)*((min(param_float(1 + 5)*(q(nx+4+0,j,indvars(5))/(qst(nx+4+0,j,indvarsst(12))*param_float(9 + 5)**2.0_wp*qst(nx+4+0,j,indvarsst(2))**2.0_wp)),10.0_wp))**6.0_wp-&
                    (min(param_float(1 + 5)*(q(nx+4+0,j,indvars(5))/(qst(nx+4+0,j,indvarsst(12))*param_float(9 + 5)**2.0_wp*qst(nx+4+0,j,indvarsst(2))**2.0_wp)),10.0_wp))))*((1.0_wp+&
                    param_float(12 + 5)**6.0_wp)/(((min(param_float(1 + 5)*(q(nx+4+0,j,indvars(5))/(qst(nx+4+0,j,indvarsst(12))*param_float(9 + 5)**2.0_wp*qst(nx+4+0,j,indvarsst(2))**2.0_wp)),10.0_wp))+&
                    param_float(11 + 5)*((min(param_float(1 + 5)*(q(nx+4+0,j,indvars(5))/(qst(nx+4+0,j,indvarsst(12))*param_float(9 + 5)**2.0_wp*qst(nx+4+0,j,indvarsst(2))**2.0_wp)),10.0_wp))**6.0_wp-&
                    (min(param_float(1 + 5)*(q(nx+4+0,j,indvars(5))/(qst(nx+4+0,j,indvarsst(12))*param_float(9 + 5)**2.0_wp*qst(nx+4+0,j,indvarsst(2))**2.0_wp)),10.0_wp))))**6.0_wp+&
                    param_float(12 + 5)**6.0_wp))**(1.0_wp/6.0_wp))-&
                    param_float(6 + 5)/param_float(9 + 5)**2*(param_float(16 + 5)*exp(-&
                    param_float(17 + 5)*(q(nx+4+0,j,indvars(5))/1.0_wp*q(nx+4+0,j,indvars(1)))**2.0_wp)))*q(nx+4+0,j,indvars(1))*q(nx+4+0,j,indvars(5))**2/qst(nx+4+0,j,indvarsst(2))**2 ) 

     enddo
