

!***********************************************************
!                                                           
! Start building layers for BC : imax j1 None **************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 0 3 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 0 3 None d(rho)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([rho*u]_1x)+deltayI*([rho*v]_1y)+(1.0_wp/c**2*((u*(-c**2*[rho]_1x+1.0_wp*[p]_1x))*deltaxI+0.5_wp*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI+esse*c*(1-M_imax*M_imax)/L_ref*(p-P0))))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_rhs_rho_dx_0_nxp4p0p01m4p3k = q(nx+4+0+0,1-4+3,indvars(1))*q(nx+4+0+0,1-4+3,indvars(2))

d1_rhs_rho_dx_0_nxp4p0m11m4p3k = q(nx+4+0-1,1-4+3,indvars(1))*q(nx+4+0-1,1-4+3,indvars(2))

d1_rhs_rho_dx_0_nxp4p0m21m4p3k = q(nx+4+0-2,1-4+3,indvars(1))*q(nx+4+0-2,1-4+3,indvars(2))

d1_rhs_rho_dx_0_nxp4p01m4p3k = 1.5_wp*d1_rhs_rho_dx_0_nxp4p0p01m4p3k-&
          2.0_wp*d1_rhs_rho_dx_0_nxp4p0m11m4p3k+&
          0.5_wp*d1_rhs_rho_dx_0_nxp4p0m21m4p3k

d1_rhs_rho_dx_0_nxp4p01m4p3k = d1_rhs_rho_dx_0_nxp4p01m4p3k*param_float(1)

d1_rhs_rho_dx_1_nxp4p0p01m4p3k = q(nx+4+0+0,1-4+3,indvars(1))

d1_rhs_rho_dx_1_nxp4p0m11m4p3k = q(nx+4+0-1,1-4+3,indvars(1))

d1_rhs_rho_dx_1_nxp4p0m21m4p3k = q(nx+4+0-2,1-4+3,indvars(1))

d1_rhs_rho_dx_1_nxp4p01m4p3k = 1.5_wp*d1_rhs_rho_dx_1_nxp4p0p01m4p3k-&
          2.0_wp*d1_rhs_rho_dx_1_nxp4p0m11m4p3k+&
          0.5_wp*d1_rhs_rho_dx_1_nxp4p0m21m4p3k

d1_rhs_rho_dx_1_nxp4p01m4p3k = d1_rhs_rho_dx_1_nxp4p01m4p3k*param_float(1)

d1_rhs_rho_dx_2_nxp4p0p01m4p3k = (param_float(3 + 5))*q(nx+4+0+0,1-4+3,indvars(1))*((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))

d1_rhs_rho_dx_2_nxp4p0m11m4p3k = (param_float(3 + 5))*q(nx+4+0-1,1-4+3,indvars(1))*((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))

d1_rhs_rho_dx_2_nxp4p0m21m4p3k = (param_float(3 + 5))*q(nx+4+0-2,1-4+3,indvars(1))*((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))

d1_rhs_rho_dx_2_nxp4p01m4p3k = 1.5_wp*d1_rhs_rho_dx_2_nxp4p0p01m4p3k-&
          2.0_wp*d1_rhs_rho_dx_2_nxp4p0m11m4p3k+&
          0.5_wp*d1_rhs_rho_dx_2_nxp4p0m21m4p3k

d1_rhs_rho_dx_2_nxp4p01m4p3k = d1_rhs_rho_dx_2_nxp4p01m4p3k*param_float(1)

d1_rhs_rho_dx_3_nxp4p0p01m4p3k = q(nx+4+0+0,1-4+3,indvars(2))

d1_rhs_rho_dx_3_nxp4p0m11m4p3k = q(nx+4+0-1,1-4+3,indvars(2))

d1_rhs_rho_dx_3_nxp4p0m21m4p3k = q(nx+4+0-2,1-4+3,indvars(2))

d1_rhs_rho_dx_3_nxp4p01m4p3k = 1.5_wp*d1_rhs_rho_dx_3_nxp4p0p01m4p3k-&
          2.0_wp*d1_rhs_rho_dx_3_nxp4p0m11m4p3k+&
          0.5_wp*d1_rhs_rho_dx_3_nxp4p0m21m4p3k

d1_rhs_rho_dx_3_nxp4p01m4p3k = d1_rhs_rho_dx_3_nxp4p01m4p3k*param_float(1)

d1_rhs_rho_dy_0_nxp4p01m4p3m1k = q(nx+4+0,1-4+3-1,indvars(1))*q(nx+4+0,1-4+3-1,indvars(3))

d1_rhs_rho_dy_0_nxp4p01m4p3p1k = q(nx+4+0,1-4+3+1,indvars(1))*q(nx+4+0,1-4+3+1,indvars(3))

d1_rhs_rho_dy_0_nxp4p01m4p3k = -&
          0.5_wp*d1_rhs_rho_dy_0_nxp4p01m4p3m1k+&
          0.5_wp*d1_rhs_rho_dy_0_nxp4p01m4p3p1k

d1_rhs_rho_dy_0_nxp4p01m4p3k = d1_rhs_rho_dy_0_nxp4p01m4p3k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 3 None d(rho)/dt *************
!                                                           
!***********************************************************


rhs(nx+4+0,1-4+3,indvars(1)) =   -  ( qst(nx+4+0,1-4+3,indvarsst(10))*(d1_rhs_rho_dx_0_nxp4p01m4p3k)+&
                    qst(nx+4+0,1-4+3,indvarsst(11))*(d1_rhs_rho_dy_0_nxp4p01m4p3k)+&
                    (1.0_wp/(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1)))**2*((q(nx+4+0,1-4+3,indvars(2))*(-&
                    (param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1)))**2*d1_rhs_rho_dx_1_nxp4p01m4p3k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p01m4p3k))*qst(nx+4+0,1-4+3,indvarsst(10))+&
                    0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1)))+&
                    q(nx+4+0,1-4+3,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1)))*q(nx+4+0,1-4+3,indvars(1))*d1_rhs_rho_dx_3_nxp4p01m4p3k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p01m4p3k))*qst(nx+4+0,1-4+3,indvarsst(10))+&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1)))*(1-&
                    qface_i(1-4+3,2)*qface_i(1-4+3,2))/param_float(20 + 5)*((param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))-&
                    param_float(26 + 5))))) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 3 None d(rho u)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! u*(1.0_wp/c**2*((u*(-c**2*[rho]_1x+1.0_wp*[p]_1x))*deltaxI+0.5_wp*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI+esse*c*(1-M_imax*M_imax)/L_ref*(p-P0))))+rho*(1.0_wp/(2*rho*c)*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI-esse*c*(1-M_imax*M_imax)/L_ref*(p-P0)))+[rho*u*v]_1y+deltaxI*([-visc_t*(2.0_wp*deltaxI*({u}_1x)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y)))]_1x)+deltayI*([-visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_rhs_u_dxdx_6_0_nxp4p0p01m4p3k_nxp4p0p0p01m4p3k = q(nx+4+0+0+0,1-4+3,indvars(2))

d2_rhs_u_dxdx_6_0_nxp4p0p01m4p3k_nxp4p0p0m11m4p3k = q(nx+4+0+0-1,1-4+3,indvars(2))

d2_rhs_u_dxdx_6_0_nxp4p0p01m4p3k_nxp4p0p0m21m4p3k = q(nx+4+0+0-2,1-4+3,indvars(2))

d2_rhs_u_dxdx_6_0_nxp4p0p01m4p3k = 1.5_wp*d2_rhs_u_dxdx_6_0_nxp4p0p01m4p3k_nxp4p0p0p01m4p3k-&
          2.0_wp*d2_rhs_u_dxdx_6_0_nxp4p0p01m4p3k_nxp4p0p0m11m4p3k+&
          0.5_wp*d2_rhs_u_dxdx_6_0_nxp4p0p01m4p3k_nxp4p0p0m21m4p3k

d2_rhs_u_dxdx_6_0_nxp4p0p01m4p3k = d2_rhs_u_dxdx_6_0_nxp4p0p01m4p3k*param_float(1)

d2_rhs_u_dxdx_6_0_nxp4p0m11m4p3k_nxp4p0m1m11m4p3k = q(nx+4+0-1-1,1-4+3,indvars(2))

d2_rhs_u_dxdx_6_0_nxp4p0m11m4p3k_nxp4p0m1p11m4p3k = q(nx+4+0-1+1,1-4+3,indvars(2))

d2_rhs_u_dxdx_6_0_nxp4p0m11m4p3k = -&
          0.5_wp*d2_rhs_u_dxdx_6_0_nxp4p0m11m4p3k_nxp4p0m1m11m4p3k+&
          0.5_wp*d2_rhs_u_dxdx_6_0_nxp4p0m11m4p3k_nxp4p0m1p11m4p3k

d2_rhs_u_dxdx_6_0_nxp4p0m11m4p3k = d2_rhs_u_dxdx_6_0_nxp4p0m11m4p3k*param_float(1)

d2_rhs_u_dxdx_6_0_nxp4p0m21m4p3k_nxp4p0m2m11m4p3k = q(nx+4+0-2-1,1-4+3,indvars(2))

d2_rhs_u_dxdx_6_0_nxp4p0m21m4p3k_nxp4p0m2p11m4p3k = q(nx+4+0-2+1,1-4+3,indvars(2))

d2_rhs_u_dxdx_6_0_nxp4p0m21m4p3k = -&
          0.5_wp*d2_rhs_u_dxdx_6_0_nxp4p0m21m4p3k_nxp4p0m2m11m4p3k+&
          0.5_wp*d2_rhs_u_dxdx_6_0_nxp4p0m21m4p3k_nxp4p0m2p11m4p3k

d2_rhs_u_dxdx_6_0_nxp4p0m21m4p3k = d2_rhs_u_dxdx_6_0_nxp4p0m21m4p3k*param_float(1)

d2_rhs_u_dxdy_6_0_nxp4p0p01m4p3k_nxp4p0p01m4p3m1k = q(nx+4+0+0,1-4+3-1,indvars(3))

d2_rhs_u_dxdy_6_0_nxp4p0p01m4p3k_nxp4p0p01m4p3p1k = q(nx+4+0+0,1-4+3+1,indvars(3))

d2_rhs_u_dxdy_6_0_nxp4p0p01m4p3k = -&
          0.5_wp*d2_rhs_u_dxdy_6_0_nxp4p0p01m4p3k_nxp4p0p01m4p3m1k+&
          0.5_wp*d2_rhs_u_dxdy_6_0_nxp4p0p01m4p3k_nxp4p0p01m4p3p1k

d2_rhs_u_dxdy_6_0_nxp4p0p01m4p3k = d2_rhs_u_dxdy_6_0_nxp4p0p01m4p3k*param_float(2)

d2_rhs_u_dxdy_6_0_nxp4p0m11m4p3k_nxp4p0m11m4p3m1k = q(nx+4+0-1,1-4+3-1,indvars(3))

d2_rhs_u_dxdy_6_0_nxp4p0m11m4p3k_nxp4p0m11m4p3p1k = q(nx+4+0-1,1-4+3+1,indvars(3))

d2_rhs_u_dxdy_6_0_nxp4p0m11m4p3k = -&
          0.5_wp*d2_rhs_u_dxdy_6_0_nxp4p0m11m4p3k_nxp4p0m11m4p3m1k+&
          0.5_wp*d2_rhs_u_dxdy_6_0_nxp4p0m11m4p3k_nxp4p0m11m4p3p1k

d2_rhs_u_dxdy_6_0_nxp4p0m11m4p3k = d2_rhs_u_dxdy_6_0_nxp4p0m11m4p3k*param_float(2)

d2_rhs_u_dxdy_6_0_nxp4p0m21m4p3k_nxp4p0m21m4p3m1k = q(nx+4+0-2,1-4+3-1,indvars(3))

d2_rhs_u_dxdy_6_0_nxp4p0m21m4p3k_nxp4p0m21m4p3p1k = q(nx+4+0-2,1-4+3+1,indvars(3))

d2_rhs_u_dxdy_6_0_nxp4p0m21m4p3k = -&
          0.5_wp*d2_rhs_u_dxdy_6_0_nxp4p0m21m4p3k_nxp4p0m21m4p3m1k+&
          0.5_wp*d2_rhs_u_dxdy_6_0_nxp4p0m21m4p3k_nxp4p0m21m4p3p1k

d2_rhs_u_dxdy_6_0_nxp4p0m21m4p3k = d2_rhs_u_dxdy_6_0_nxp4p0m21m4p3k*param_float(2)

d1_rhs_u_dx_6_nxp4p0p01m4p3k = -((1+&
                    param_float(21 + 5))/(((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+4+0+0,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0+0,1-4+3,indvars(1)))**3/((q(nx+4+0+0,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0+0,1-4+3,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+4+0+0,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0+0,1-4+3,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0+0,1-4+3,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0p01m4p3k)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0+0,1-4+3,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0p01m4p3k)+&
                    qst(nx+4+0+0,1-4+3,indvarsst(11))*(d2_rhs_u_dxdy_6_0_nxp4p0p01m4p3k)))

d1_rhs_u_dx_6_nxp4p0m11m4p3k = -((1+&
                    param_float(21 + 5))/(((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+4+0-1,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0-1,1-4+3,indvars(1)))**3/((q(nx+4+0-1,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0-1,1-4+3,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+4+0-1,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0-1,1-4+3,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0-1,1-4+3,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0m11m4p3k)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0-1,1-4+3,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0m11m4p3k)+&
                    qst(nx+4+0-1,1-4+3,indvarsst(11))*(d2_rhs_u_dxdy_6_0_nxp4p0m11m4p3k)))

d1_rhs_u_dx_6_nxp4p0m21m4p3k = -((1+&
                    param_float(21 + 5))/(((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+4+0-2,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0-2,1-4+3,indvars(1)))**3/((q(nx+4+0-2,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0-2,1-4+3,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+4+0-2,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0-2,1-4+3,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0-2,1-4+3,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0m21m4p3k)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0-2,1-4+3,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0m21m4p3k)+&
                    qst(nx+4+0-2,1-4+3,indvarsst(11))*(d2_rhs_u_dxdy_6_0_nxp4p0m21m4p3k)))

d1_rhs_u_dx_6_nxp4p01m4p3k = 1.5_wp*d1_rhs_u_dx_6_nxp4p0p01m4p3k-&
          2.0_wp*d1_rhs_u_dx_6_nxp4p0m11m4p3k+&
          0.5_wp*d1_rhs_u_dx_6_nxp4p0m21m4p3k

d1_rhs_u_dx_6_nxp4p01m4p3k = d1_rhs_u_dx_6_nxp4p01m4p3k*param_float(1)

d1_rhs_u_dy_0_nxp4p01m4p3m1k = q(nx+4+0,1-4+3-1,indvars(1))*q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(3))

d1_rhs_u_dy_0_nxp4p01m4p3p1k = q(nx+4+0,1-4+3+1,indvars(1))*q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(3))

d1_rhs_u_dy_0_nxp4p01m4p3k = -&
          0.5_wp*d1_rhs_u_dy_0_nxp4p01m4p3m1k+&
          0.5_wp*d1_rhs_u_dy_0_nxp4p01m4p3p1k

d1_rhs_u_dy_0_nxp4p01m4p3k = d1_rhs_u_dy_0_nxp4p01m4p3k*param_float(2)

d2_rhs_u_dydx_1_0_nxp4p01m4p3m3k_nxp4p0p01m4p3m3k = q(nx+4+0+0,1-4+3-3,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p01m4p3m3k_nxp4p0m11m4p3m3k = q(nx+4+0-1,1-4+3-3,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p01m4p3m3k_nxp4p0m21m4p3m3k = q(nx+4+0-2,1-4+3-3,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p01m4p3m3k = 1.5_wp*d2_rhs_u_dydx_1_0_nxp4p01m4p3m3k_nxp4p0p01m4p3m3k-&
          2.0_wp*d2_rhs_u_dydx_1_0_nxp4p01m4p3m3k_nxp4p0m11m4p3m3k+&
          0.5_wp*d2_rhs_u_dydx_1_0_nxp4p01m4p3m3k_nxp4p0m21m4p3m3k

d2_rhs_u_dydx_1_0_nxp4p01m4p3m3k = d2_rhs_u_dydx_1_0_nxp4p01m4p3m3k*param_float(1)

d2_rhs_u_dydx_1_0_nxp4p01m4p3m2k_nxp4p0p01m4p3m2k = q(nx+4+0+0,1-4+3-2,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p01m4p3m2k_nxp4p0m11m4p3m2k = q(nx+4+0-1,1-4+3-2,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p01m4p3m2k_nxp4p0m21m4p3m2k = q(nx+4+0-2,1-4+3-2,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p01m4p3m2k = 1.5_wp*d2_rhs_u_dydx_1_0_nxp4p01m4p3m2k_nxp4p0p01m4p3m2k-&
          2.0_wp*d2_rhs_u_dydx_1_0_nxp4p01m4p3m2k_nxp4p0m11m4p3m2k+&
          0.5_wp*d2_rhs_u_dydx_1_0_nxp4p01m4p3m2k_nxp4p0m21m4p3m2k

d2_rhs_u_dydx_1_0_nxp4p01m4p3m2k = d2_rhs_u_dydx_1_0_nxp4p01m4p3m2k*param_float(1)

d2_rhs_u_dydx_1_0_nxp4p01m4p3m1k_nxp4p0p01m4p3m1k = q(nx+4+0+0,1-4+3-1,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p01m4p3m1k_nxp4p0m11m4p3m1k = q(nx+4+0-1,1-4+3-1,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p01m4p3m1k_nxp4p0m21m4p3m1k = q(nx+4+0-2,1-4+3-1,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p01m4p3m1k = 1.5_wp*d2_rhs_u_dydx_1_0_nxp4p01m4p3m1k_nxp4p0p01m4p3m1k-&
          2.0_wp*d2_rhs_u_dydx_1_0_nxp4p01m4p3m1k_nxp4p0m11m4p3m1k+&
          0.5_wp*d2_rhs_u_dydx_1_0_nxp4p01m4p3m1k_nxp4p0m21m4p3m1k

d2_rhs_u_dydx_1_0_nxp4p01m4p3m1k = d2_rhs_u_dydx_1_0_nxp4p01m4p3m1k*param_float(1)

d2_rhs_u_dydx_1_0_nxp4p01m4p3p1k_nxp4p0p01m4p3p1k = q(nx+4+0+0,1-4+3+1,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p01m4p3p1k_nxp4p0m11m4p3p1k = q(nx+4+0-1,1-4+3+1,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p01m4p3p1k_nxp4p0m21m4p3p1k = q(nx+4+0-2,1-4+3+1,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p01m4p3p1k = 1.5_wp*d2_rhs_u_dydx_1_0_nxp4p01m4p3p1k_nxp4p0p01m4p3p1k-&
          2.0_wp*d2_rhs_u_dydx_1_0_nxp4p01m4p3p1k_nxp4p0m11m4p3p1k+&
          0.5_wp*d2_rhs_u_dydx_1_0_nxp4p01m4p3p1k_nxp4p0m21m4p3p1k

d2_rhs_u_dydx_1_0_nxp4p01m4p3p1k = d2_rhs_u_dydx_1_0_nxp4p01m4p3p1k*param_float(1)

d2_rhs_u_dydx_1_0_nxp4p01m4p3p2k_nxp4p0p01m4p3p2k = q(nx+4+0+0,1-4+3+2,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p01m4p3p2k_nxp4p0m11m4p3p2k = q(nx+4+0-1,1-4+3+2,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p01m4p3p2k_nxp4p0m21m4p3p2k = q(nx+4+0-2,1-4+3+2,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p01m4p3p2k = 1.5_wp*d2_rhs_u_dydx_1_0_nxp4p01m4p3p2k_nxp4p0p01m4p3p2k-&
          2.0_wp*d2_rhs_u_dydx_1_0_nxp4p01m4p3p2k_nxp4p0m11m4p3p2k+&
          0.5_wp*d2_rhs_u_dydx_1_0_nxp4p01m4p3p2k_nxp4p0m21m4p3p2k

d2_rhs_u_dydx_1_0_nxp4p01m4p3p2k = d2_rhs_u_dydx_1_0_nxp4p01m4p3p2k*param_float(1)

d2_rhs_u_dydx_1_0_nxp4p01m4p3p3k_nxp4p0p01m4p3p3k = q(nx+4+0+0,1-4+3+3,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p01m4p3p3k_nxp4p0m11m4p3p3k = q(nx+4+0-1,1-4+3+3,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p01m4p3p3k_nxp4p0m21m4p3p3k = q(nx+4+0-2,1-4+3+3,indvars(3))

d2_rhs_u_dydx_1_0_nxp4p01m4p3p3k = 1.5_wp*d2_rhs_u_dydx_1_0_nxp4p01m4p3p3k_nxp4p0p01m4p3p3k-&
          2.0_wp*d2_rhs_u_dydx_1_0_nxp4p01m4p3p3k_nxp4p0m11m4p3p3k+&
          0.5_wp*d2_rhs_u_dydx_1_0_nxp4p01m4p3p3k_nxp4p0m21m4p3p3k

d2_rhs_u_dydx_1_0_nxp4p01m4p3p3k = d2_rhs_u_dydx_1_0_nxp4p01m4p3p3k*param_float(1)

d2_rhs_u_dydy_1_0_nxp4p01m4p3m3k_nxp4p01m4p3m3p0k = q(nx+4+0,1-4+3-3+0,indvars(2))

d2_rhs_u_dydy_1_0_nxp4p01m4p3m3k_nxp4p01m4p3m3p1k = q(nx+4+0,1-4+3-3+1,indvars(2))

d2_rhs_u_dydy_1_0_nxp4p01m4p3m3k_nxp4p01m4p3m3p2k = q(nx+4+0,1-4+3-3+2,indvars(2))

d2_rhs_u_dydy_1_0_nxp4p01m4p3m3k = -&
          1.5_wp*d2_rhs_u_dydy_1_0_nxp4p01m4p3m3k_nxp4p01m4p3m3p0k+&
          2.0_wp*d2_rhs_u_dydy_1_0_nxp4p01m4p3m3k_nxp4p01m4p3m3p1k-&
          0.5_wp*d2_rhs_u_dydy_1_0_nxp4p01m4p3m3k_nxp4p01m4p3m3p2k

d2_rhs_u_dydy_1_0_nxp4p01m4p3m3k = d2_rhs_u_dydy_1_0_nxp4p01m4p3m3k*param_float(2)

d2_rhs_u_dydy_1_0_nxp4p01m4p3m2k_nxp4p01m4p3m2m1k = q(nx+4+0,1-4+3-2-1,indvars(2))

d2_rhs_u_dydy_1_0_nxp4p01m4p3m2k_nxp4p01m4p3m2p1k = q(nx+4+0,1-4+3-2+1,indvars(2))

d2_rhs_u_dydy_1_0_nxp4p01m4p3m2k = -&
          0.5_wp*d2_rhs_u_dydy_1_0_nxp4p01m4p3m2k_nxp4p01m4p3m2m1k+&
          0.5_wp*d2_rhs_u_dydy_1_0_nxp4p01m4p3m2k_nxp4p01m4p3m2p1k

d2_rhs_u_dydy_1_0_nxp4p01m4p3m2k = d2_rhs_u_dydy_1_0_nxp4p01m4p3m2k*param_float(2)

d2_rhs_u_dydy_1_0_nxp4p01m4p3m1k_nxp4p01m4p3m1m1k = q(nx+4+0,1-4+3-1-1,indvars(2))

d2_rhs_u_dydy_1_0_nxp4p01m4p3m1k_nxp4p01m4p3m1p1k = q(nx+4+0,1-4+3-1+1,indvars(2))

d2_rhs_u_dydy_1_0_nxp4p01m4p3m1k = -&
          0.5_wp*d2_rhs_u_dydy_1_0_nxp4p01m4p3m1k_nxp4p01m4p3m1m1k+&
          0.5_wp*d2_rhs_u_dydy_1_0_nxp4p01m4p3m1k_nxp4p01m4p3m1p1k

d2_rhs_u_dydy_1_0_nxp4p01m4p3m1k = d2_rhs_u_dydy_1_0_nxp4p01m4p3m1k*param_float(2)

d2_rhs_u_dydy_1_0_nxp4p01m4p3p1k_nxp4p01m4p3p1m1k = q(nx+4+0,1-4+3+1-1,indvars(2))

d2_rhs_u_dydy_1_0_nxp4p01m4p3p1k_nxp4p01m4p3p1p1k = q(nx+4+0,1-4+3+1+1,indvars(2))

d2_rhs_u_dydy_1_0_nxp4p01m4p3p1k = -&
          0.5_wp*d2_rhs_u_dydy_1_0_nxp4p01m4p3p1k_nxp4p01m4p3p1m1k+&
          0.5_wp*d2_rhs_u_dydy_1_0_nxp4p01m4p3p1k_nxp4p01m4p3p1p1k

d2_rhs_u_dydy_1_0_nxp4p01m4p3p1k = d2_rhs_u_dydy_1_0_nxp4p01m4p3p1k*param_float(2)

d2_rhs_u_dydy_1_0_nxp4p01m4p3p2k_nxp4p01m4p3p2m1k = q(nx+4+0,1-4+3+2-1,indvars(2))

d2_rhs_u_dydy_1_0_nxp4p01m4p3p2k_nxp4p01m4p3p2p1k = q(nx+4+0,1-4+3+2+1,indvars(2))

d2_rhs_u_dydy_1_0_nxp4p01m4p3p2k = -&
          0.5_wp*d2_rhs_u_dydy_1_0_nxp4p01m4p3p2k_nxp4p01m4p3p2m1k+&
          0.5_wp*d2_rhs_u_dydy_1_0_nxp4p01m4p3p2k_nxp4p01m4p3p2p1k

d2_rhs_u_dydy_1_0_nxp4p01m4p3p2k = d2_rhs_u_dydy_1_0_nxp4p01m4p3p2k*param_float(2)

d2_rhs_u_dydy_1_0_nxp4p01m4p3p3k_nxp4p01m4p3p3m1k = q(nx+4+0,1-4+3+3-1,indvars(2))

d2_rhs_u_dydy_1_0_nxp4p01m4p3p3k_nxp4p01m4p3p3p1k = q(nx+4+0,1-4+3+3+1,indvars(2))

d2_rhs_u_dydy_1_0_nxp4p01m4p3p3k = -&
          0.5_wp*d2_rhs_u_dydy_1_0_nxp4p01m4p3p3k_nxp4p01m4p3p3m1k+&
          0.5_wp*d2_rhs_u_dydy_1_0_nxp4p01m4p3p3k_nxp4p01m4p3p3p1k

d2_rhs_u_dydy_1_0_nxp4p01m4p3p3k = d2_rhs_u_dydy_1_0_nxp4p01m4p3p3k*param_float(2)

d1_rhs_u_dy_1_nxp4p01m4p3m1k = -((1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+4+0,1-4+3-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0,1-4+3-1,indvars(1)))**3/((q(nx+4+0,1-4+3-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0,1-4+3-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+4+0,1-4+3-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0,1-4+3-1,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0,1-4+3-1,indvarsst(11))*(d2_rhs_u_dydy_1_0_nxp4p01m4p3m1k)+&
                    qst(nx+4+0,1-4+3-1,indvarsst(10))*(d2_rhs_u_dydx_1_0_nxp4p01m4p3m1k))

d1_rhs_u_dy_1_nxp4p01m4p3p1k = -((1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+4+0,1-4+3+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0,1-4+3+1,indvars(1)))**3/((q(nx+4+0,1-4+3+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0,1-4+3+1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+4+0,1-4+3+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0,1-4+3+1,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0,1-4+3+1,indvarsst(11))*(d2_rhs_u_dydy_1_0_nxp4p01m4p3p1k)+&
                    qst(nx+4+0,1-4+3+1,indvarsst(10))*(d2_rhs_u_dydx_1_0_nxp4p01m4p3p1k))

d1_rhs_u_dy_1_nxp4p01m4p3k = -&
          0.5_wp*d1_rhs_u_dy_1_nxp4p01m4p3m1k+&
          0.5_wp*d1_rhs_u_dy_1_nxp4p01m4p3p1k

d1_rhs_u_dy_1_nxp4p01m4p3k = d1_rhs_u_dy_1_nxp4p01m4p3k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 3 None d(rho u)/dt ***********
!                                                           
!***********************************************************


rhs(nx+4+0,1-4+3,indvars(2)) =   -  ( q(nx+4+0,1-4+3,indvars(2))*(1.0_wp/(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1)))**2*((q(nx+4+0,1-4+3,indvars(2))*(-&
                    (param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1)))**2*d1_rhs_rho_dx_1_nxp4p01m4p3k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p01m4p3k))*qst(nx+4+0,1-4+3,indvarsst(10))+&
                    0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1)))+&
                    q(nx+4+0,1-4+3,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1)))*q(nx+4+0,1-4+3,indvars(1))*d1_rhs_rho_dx_3_nxp4p01m4p3k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p01m4p3k))*qst(nx+4+0,1-4+3,indvarsst(10))+&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1)))*(1-&
                    qface_i(1-4+3,2)*qface_i(1-4+3,2))/param_float(20 + 5)*((param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))-&
                    param_float(26 + 5)))))+&
                    q(nx+4+0,1-4+3,indvars(1))*(1.0_wp/(2*q(nx+4+0,1-4+3,indvars(1))*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1))))*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1)))+&
                    q(nx+4+0,1-4+3,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1)))*q(nx+4+0,1-4+3,indvars(1))*d1_rhs_rho_dx_3_nxp4p01m4p3k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p01m4p3k))*qst(nx+4+0,1-4+3,indvarsst(10))-&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1)))*(1-&
                    qface_i(1-4+3,2)*qface_i(1-4+3,2))/param_float(20 + 5)*((param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))-&
                    param_float(26 + 5))))+&
                    d1_rhs_u_dy_0_nxp4p01m4p3k+&
                    qst(nx+4+0,1-4+3,indvarsst(10))*(d1_rhs_u_dx_6_nxp4p01m4p3k)+&
                    qst(nx+4+0,1-4+3,indvarsst(11))*(d1_rhs_u_dy_1_nxp4p01m4p3k) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 3 None d(rho v)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! v*(1.0_wp/c**2*((u*(-c**2*[rho]_1x+1.0_wp*[p]_1x))*deltaxI+0.5_wp*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI+esse*c*(1-M_imax*M_imax)/L_ref*(p-P0))))+rho*((1.0_wp*u*[v]_1x)*deltaxI)+[rho*v*v]_1y+deltaxI*([-visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x))]_1x)+deltayI*([-visc_t*(2.0_wp*deltayI*({v}_1y)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y)))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_rhs_v_dx_4_nxp4p0p01m4p3k = q(nx+4+0+0,1-4+3,indvars(3))

d1_rhs_v_dx_4_nxp4p0m11m4p3k = q(nx+4+0-1,1-4+3,indvars(3))

d1_rhs_v_dx_4_nxp4p0m21m4p3k = q(nx+4+0-2,1-4+3,indvars(3))

d1_rhs_v_dx_4_nxp4p01m4p3k = 1.5_wp*d1_rhs_v_dx_4_nxp4p0p01m4p3k-&
          2.0_wp*d1_rhs_v_dx_4_nxp4p0m11m4p3k+&
          0.5_wp*d1_rhs_v_dx_4_nxp4p0m21m4p3k

d1_rhs_v_dx_4_nxp4p01m4p3k = d1_rhs_v_dx_4_nxp4p01m4p3k*param_float(1)

d2_rhs_v_dxdx_5_0_nxp4p0p01m4p3k_nxp4p0p0p01m4p3k = q(nx+4+0+0+0,1-4+3,indvars(3))

d2_rhs_v_dxdx_5_0_nxp4p0p01m4p3k_nxp4p0p0m11m4p3k = q(nx+4+0+0-1,1-4+3,indvars(3))

d2_rhs_v_dxdx_5_0_nxp4p0p01m4p3k_nxp4p0p0m21m4p3k = q(nx+4+0+0-2,1-4+3,indvars(3))

d2_rhs_v_dxdx_5_0_nxp4p0p01m4p3k = 1.5_wp*d2_rhs_v_dxdx_5_0_nxp4p0p01m4p3k_nxp4p0p0p01m4p3k-&
          2.0_wp*d2_rhs_v_dxdx_5_0_nxp4p0p01m4p3k_nxp4p0p0m11m4p3k+&
          0.5_wp*d2_rhs_v_dxdx_5_0_nxp4p0p01m4p3k_nxp4p0p0m21m4p3k

d2_rhs_v_dxdx_5_0_nxp4p0p01m4p3k = d2_rhs_v_dxdx_5_0_nxp4p0p01m4p3k*param_float(1)

d2_rhs_v_dxdx_5_0_nxp4p0m11m4p3k_nxp4p0m1m11m4p3k = q(nx+4+0-1-1,1-4+3,indvars(3))

d2_rhs_v_dxdx_5_0_nxp4p0m11m4p3k_nxp4p0m1p11m4p3k = q(nx+4+0-1+1,1-4+3,indvars(3))

d2_rhs_v_dxdx_5_0_nxp4p0m11m4p3k = -&
          0.5_wp*d2_rhs_v_dxdx_5_0_nxp4p0m11m4p3k_nxp4p0m1m11m4p3k+&
          0.5_wp*d2_rhs_v_dxdx_5_0_nxp4p0m11m4p3k_nxp4p0m1p11m4p3k

d2_rhs_v_dxdx_5_0_nxp4p0m11m4p3k = d2_rhs_v_dxdx_5_0_nxp4p0m11m4p3k*param_float(1)

d2_rhs_v_dxdx_5_0_nxp4p0m21m4p3k_nxp4p0m2m11m4p3k = q(nx+4+0-2-1,1-4+3,indvars(3))

d2_rhs_v_dxdx_5_0_nxp4p0m21m4p3k_nxp4p0m2p11m4p3k = q(nx+4+0-2+1,1-4+3,indvars(3))

d2_rhs_v_dxdx_5_0_nxp4p0m21m4p3k = -&
          0.5_wp*d2_rhs_v_dxdx_5_0_nxp4p0m21m4p3k_nxp4p0m2m11m4p3k+&
          0.5_wp*d2_rhs_v_dxdx_5_0_nxp4p0m21m4p3k_nxp4p0m2p11m4p3k

d2_rhs_v_dxdx_5_0_nxp4p0m21m4p3k = d2_rhs_v_dxdx_5_0_nxp4p0m21m4p3k*param_float(1)

d2_rhs_v_dxdy_5_0_nxp4p0p01m4p3k_nxp4p0p01m4p3m1k = q(nx+4+0+0,1-4+3-1,indvars(2))

d2_rhs_v_dxdy_5_0_nxp4p0p01m4p3k_nxp4p0p01m4p3p1k = q(nx+4+0+0,1-4+3+1,indvars(2))

d2_rhs_v_dxdy_5_0_nxp4p0p01m4p3k = -&
          0.5_wp*d2_rhs_v_dxdy_5_0_nxp4p0p01m4p3k_nxp4p0p01m4p3m1k+&
          0.5_wp*d2_rhs_v_dxdy_5_0_nxp4p0p01m4p3k_nxp4p0p01m4p3p1k

d2_rhs_v_dxdy_5_0_nxp4p0p01m4p3k = d2_rhs_v_dxdy_5_0_nxp4p0p01m4p3k*param_float(2)

d2_rhs_v_dxdy_5_0_nxp4p0m11m4p3k_nxp4p0m11m4p3m1k = q(nx+4+0-1,1-4+3-1,indvars(2))

d2_rhs_v_dxdy_5_0_nxp4p0m11m4p3k_nxp4p0m11m4p3p1k = q(nx+4+0-1,1-4+3+1,indvars(2))

d2_rhs_v_dxdy_5_0_nxp4p0m11m4p3k = -&
          0.5_wp*d2_rhs_v_dxdy_5_0_nxp4p0m11m4p3k_nxp4p0m11m4p3m1k+&
          0.5_wp*d2_rhs_v_dxdy_5_0_nxp4p0m11m4p3k_nxp4p0m11m4p3p1k

d2_rhs_v_dxdy_5_0_nxp4p0m11m4p3k = d2_rhs_v_dxdy_5_0_nxp4p0m11m4p3k*param_float(2)

d2_rhs_v_dxdy_5_0_nxp4p0m21m4p3k_nxp4p0m21m4p3m1k = q(nx+4+0-2,1-4+3-1,indvars(2))

d2_rhs_v_dxdy_5_0_nxp4p0m21m4p3k_nxp4p0m21m4p3p1k = q(nx+4+0-2,1-4+3+1,indvars(2))

d2_rhs_v_dxdy_5_0_nxp4p0m21m4p3k = -&
          0.5_wp*d2_rhs_v_dxdy_5_0_nxp4p0m21m4p3k_nxp4p0m21m4p3m1k+&
          0.5_wp*d2_rhs_v_dxdy_5_0_nxp4p0m21m4p3k_nxp4p0m21m4p3p1k

d2_rhs_v_dxdy_5_0_nxp4p0m21m4p3k = d2_rhs_v_dxdy_5_0_nxp4p0m21m4p3k*param_float(2)

d1_rhs_v_dx_5_nxp4p0p01m4p3k = -((1+&
                    param_float(21 + 5))/(((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+4+0+0,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0+0,1-4+3,indvars(1)))**3/((q(nx+4+0+0,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0+0,1-4+3,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+4+0+0,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0+0,1-4+3,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0+0,1-4+3,indvarsst(11))*(d2_rhs_v_dxdy_5_0_nxp4p0p01m4p3k)+&
                    qst(nx+4+0+0,1-4+3,indvarsst(10))*(d2_rhs_v_dxdx_5_0_nxp4p0p01m4p3k))

d1_rhs_v_dx_5_nxp4p0m11m4p3k = -((1+&
                    param_float(21 + 5))/(((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+4+0-1,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0-1,1-4+3,indvars(1)))**3/((q(nx+4+0-1,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0-1,1-4+3,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+4+0-1,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0-1,1-4+3,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0-1,1-4+3,indvarsst(11))*(d2_rhs_v_dxdy_5_0_nxp4p0m11m4p3k)+&
                    qst(nx+4+0-1,1-4+3,indvarsst(10))*(d2_rhs_v_dxdx_5_0_nxp4p0m11m4p3k))

d1_rhs_v_dx_5_nxp4p0m21m4p3k = -((1+&
                    param_float(21 + 5))/(((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+4+0-2,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0-2,1-4+3,indvars(1)))**3/((q(nx+4+0-2,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0-2,1-4+3,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+4+0-2,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0-2,1-4+3,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0-2,1-4+3,indvarsst(11))*(d2_rhs_v_dxdy_5_0_nxp4p0m21m4p3k)+&
                    qst(nx+4+0-2,1-4+3,indvarsst(10))*(d2_rhs_v_dxdx_5_0_nxp4p0m21m4p3k))

d1_rhs_v_dx_5_nxp4p01m4p3k = 1.5_wp*d1_rhs_v_dx_5_nxp4p0p01m4p3k-&
          2.0_wp*d1_rhs_v_dx_5_nxp4p0m11m4p3k+&
          0.5_wp*d1_rhs_v_dx_5_nxp4p0m21m4p3k

d1_rhs_v_dx_5_nxp4p01m4p3k = d1_rhs_v_dx_5_nxp4p01m4p3k*param_float(1)

d1_rhs_v_dy_0_nxp4p01m4p3m1k = q(nx+4+0,1-4+3-1,indvars(1))*q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3))

d1_rhs_v_dy_0_nxp4p01m4p3p1k = q(nx+4+0,1-4+3+1,indvars(1))*q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3))

d1_rhs_v_dy_0_nxp4p01m4p3k = -&
          0.5_wp*d1_rhs_v_dy_0_nxp4p01m4p3m1k+&
          0.5_wp*d1_rhs_v_dy_0_nxp4p01m4p3p1k

d1_rhs_v_dy_0_nxp4p01m4p3k = d1_rhs_v_dy_0_nxp4p01m4p3k*param_float(2)

d2_rhs_v_dydx_1_0_nxp4p01m4p3m3k_nxp4p0p01m4p3m3k = q(nx+4+0+0,1-4+3-3,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p01m4p3m3k_nxp4p0m11m4p3m3k = q(nx+4+0-1,1-4+3-3,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p01m4p3m3k_nxp4p0m21m4p3m3k = q(nx+4+0-2,1-4+3-3,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p01m4p3m3k = 1.5_wp*d2_rhs_v_dydx_1_0_nxp4p01m4p3m3k_nxp4p0p01m4p3m3k-&
          2.0_wp*d2_rhs_v_dydx_1_0_nxp4p01m4p3m3k_nxp4p0m11m4p3m3k+&
          0.5_wp*d2_rhs_v_dydx_1_0_nxp4p01m4p3m3k_nxp4p0m21m4p3m3k

d2_rhs_v_dydx_1_0_nxp4p01m4p3m3k = d2_rhs_v_dydx_1_0_nxp4p01m4p3m3k*param_float(1)

d2_rhs_v_dydx_1_0_nxp4p01m4p3m2k_nxp4p0p01m4p3m2k = q(nx+4+0+0,1-4+3-2,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p01m4p3m2k_nxp4p0m11m4p3m2k = q(nx+4+0-1,1-4+3-2,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p01m4p3m2k_nxp4p0m21m4p3m2k = q(nx+4+0-2,1-4+3-2,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p01m4p3m2k = 1.5_wp*d2_rhs_v_dydx_1_0_nxp4p01m4p3m2k_nxp4p0p01m4p3m2k-&
          2.0_wp*d2_rhs_v_dydx_1_0_nxp4p01m4p3m2k_nxp4p0m11m4p3m2k+&
          0.5_wp*d2_rhs_v_dydx_1_0_nxp4p01m4p3m2k_nxp4p0m21m4p3m2k

d2_rhs_v_dydx_1_0_nxp4p01m4p3m2k = d2_rhs_v_dydx_1_0_nxp4p01m4p3m2k*param_float(1)

d2_rhs_v_dydx_1_0_nxp4p01m4p3m1k_nxp4p0p01m4p3m1k = q(nx+4+0+0,1-4+3-1,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p01m4p3m1k_nxp4p0m11m4p3m1k = q(nx+4+0-1,1-4+3-1,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p01m4p3m1k_nxp4p0m21m4p3m1k = q(nx+4+0-2,1-4+3-1,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p01m4p3m1k = 1.5_wp*d2_rhs_v_dydx_1_0_nxp4p01m4p3m1k_nxp4p0p01m4p3m1k-&
          2.0_wp*d2_rhs_v_dydx_1_0_nxp4p01m4p3m1k_nxp4p0m11m4p3m1k+&
          0.5_wp*d2_rhs_v_dydx_1_0_nxp4p01m4p3m1k_nxp4p0m21m4p3m1k

d2_rhs_v_dydx_1_0_nxp4p01m4p3m1k = d2_rhs_v_dydx_1_0_nxp4p01m4p3m1k*param_float(1)

d2_rhs_v_dydx_1_0_nxp4p01m4p3p1k_nxp4p0p01m4p3p1k = q(nx+4+0+0,1-4+3+1,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p01m4p3p1k_nxp4p0m11m4p3p1k = q(nx+4+0-1,1-4+3+1,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p01m4p3p1k_nxp4p0m21m4p3p1k = q(nx+4+0-2,1-4+3+1,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p01m4p3p1k = 1.5_wp*d2_rhs_v_dydx_1_0_nxp4p01m4p3p1k_nxp4p0p01m4p3p1k-&
          2.0_wp*d2_rhs_v_dydx_1_0_nxp4p01m4p3p1k_nxp4p0m11m4p3p1k+&
          0.5_wp*d2_rhs_v_dydx_1_0_nxp4p01m4p3p1k_nxp4p0m21m4p3p1k

d2_rhs_v_dydx_1_0_nxp4p01m4p3p1k = d2_rhs_v_dydx_1_0_nxp4p01m4p3p1k*param_float(1)

d2_rhs_v_dydx_1_0_nxp4p01m4p3p2k_nxp4p0p01m4p3p2k = q(nx+4+0+0,1-4+3+2,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p01m4p3p2k_nxp4p0m11m4p3p2k = q(nx+4+0-1,1-4+3+2,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p01m4p3p2k_nxp4p0m21m4p3p2k = q(nx+4+0-2,1-4+3+2,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p01m4p3p2k = 1.5_wp*d2_rhs_v_dydx_1_0_nxp4p01m4p3p2k_nxp4p0p01m4p3p2k-&
          2.0_wp*d2_rhs_v_dydx_1_0_nxp4p01m4p3p2k_nxp4p0m11m4p3p2k+&
          0.5_wp*d2_rhs_v_dydx_1_0_nxp4p01m4p3p2k_nxp4p0m21m4p3p2k

d2_rhs_v_dydx_1_0_nxp4p01m4p3p2k = d2_rhs_v_dydx_1_0_nxp4p01m4p3p2k*param_float(1)

d2_rhs_v_dydx_1_0_nxp4p01m4p3p3k_nxp4p0p01m4p3p3k = q(nx+4+0+0,1-4+3+3,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p01m4p3p3k_nxp4p0m11m4p3p3k = q(nx+4+0-1,1-4+3+3,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p01m4p3p3k_nxp4p0m21m4p3p3k = q(nx+4+0-2,1-4+3+3,indvars(2))

d2_rhs_v_dydx_1_0_nxp4p01m4p3p3k = 1.5_wp*d2_rhs_v_dydx_1_0_nxp4p01m4p3p3k_nxp4p0p01m4p3p3k-&
          2.0_wp*d2_rhs_v_dydx_1_0_nxp4p01m4p3p3k_nxp4p0m11m4p3p3k+&
          0.5_wp*d2_rhs_v_dydx_1_0_nxp4p01m4p3p3k_nxp4p0m21m4p3p3k

d2_rhs_v_dydx_1_0_nxp4p01m4p3p3k = d2_rhs_v_dydx_1_0_nxp4p01m4p3p3k*param_float(1)

d2_rhs_v_dydy_1_0_nxp4p01m4p3m3k_nxp4p01m4p3m3p0k = q(nx+4+0,1-4+3-3+0,indvars(3))

d2_rhs_v_dydy_1_0_nxp4p01m4p3m3k_nxp4p01m4p3m3p1k = q(nx+4+0,1-4+3-3+1,indvars(3))

d2_rhs_v_dydy_1_0_nxp4p01m4p3m3k_nxp4p01m4p3m3p2k = q(nx+4+0,1-4+3-3+2,indvars(3))

d2_rhs_v_dydy_1_0_nxp4p01m4p3m3k = -&
          1.5_wp*d2_rhs_v_dydy_1_0_nxp4p01m4p3m3k_nxp4p01m4p3m3p0k+&
          2.0_wp*d2_rhs_v_dydy_1_0_nxp4p01m4p3m3k_nxp4p01m4p3m3p1k-&
          0.5_wp*d2_rhs_v_dydy_1_0_nxp4p01m4p3m3k_nxp4p01m4p3m3p2k

d2_rhs_v_dydy_1_0_nxp4p01m4p3m3k = d2_rhs_v_dydy_1_0_nxp4p01m4p3m3k*param_float(2)

d2_rhs_v_dydy_1_0_nxp4p01m4p3m2k_nxp4p01m4p3m2m1k = q(nx+4+0,1-4+3-2-1,indvars(3))

d2_rhs_v_dydy_1_0_nxp4p01m4p3m2k_nxp4p01m4p3m2p1k = q(nx+4+0,1-4+3-2+1,indvars(3))

d2_rhs_v_dydy_1_0_nxp4p01m4p3m2k = -&
          0.5_wp*d2_rhs_v_dydy_1_0_nxp4p01m4p3m2k_nxp4p01m4p3m2m1k+&
          0.5_wp*d2_rhs_v_dydy_1_0_nxp4p01m4p3m2k_nxp4p01m4p3m2p1k

d2_rhs_v_dydy_1_0_nxp4p01m4p3m2k = d2_rhs_v_dydy_1_0_nxp4p01m4p3m2k*param_float(2)

d2_rhs_v_dydy_1_0_nxp4p01m4p3m1k_nxp4p01m4p3m1m1k = q(nx+4+0,1-4+3-1-1,indvars(3))

d2_rhs_v_dydy_1_0_nxp4p01m4p3m1k_nxp4p01m4p3m1p1k = q(nx+4+0,1-4+3-1+1,indvars(3))

d2_rhs_v_dydy_1_0_nxp4p01m4p3m1k = -&
          0.5_wp*d2_rhs_v_dydy_1_0_nxp4p01m4p3m1k_nxp4p01m4p3m1m1k+&
          0.5_wp*d2_rhs_v_dydy_1_0_nxp4p01m4p3m1k_nxp4p01m4p3m1p1k

d2_rhs_v_dydy_1_0_nxp4p01m4p3m1k = d2_rhs_v_dydy_1_0_nxp4p01m4p3m1k*param_float(2)

d2_rhs_v_dydy_1_0_nxp4p01m4p3p1k_nxp4p01m4p3p1m1k = q(nx+4+0,1-4+3+1-1,indvars(3))

d2_rhs_v_dydy_1_0_nxp4p01m4p3p1k_nxp4p01m4p3p1p1k = q(nx+4+0,1-4+3+1+1,indvars(3))

d2_rhs_v_dydy_1_0_nxp4p01m4p3p1k = -&
          0.5_wp*d2_rhs_v_dydy_1_0_nxp4p01m4p3p1k_nxp4p01m4p3p1m1k+&
          0.5_wp*d2_rhs_v_dydy_1_0_nxp4p01m4p3p1k_nxp4p01m4p3p1p1k

d2_rhs_v_dydy_1_0_nxp4p01m4p3p1k = d2_rhs_v_dydy_1_0_nxp4p01m4p3p1k*param_float(2)

d2_rhs_v_dydy_1_0_nxp4p01m4p3p2k_nxp4p01m4p3p2m1k = q(nx+4+0,1-4+3+2-1,indvars(3))

d2_rhs_v_dydy_1_0_nxp4p01m4p3p2k_nxp4p01m4p3p2p1k = q(nx+4+0,1-4+3+2+1,indvars(3))

d2_rhs_v_dydy_1_0_nxp4p01m4p3p2k = -&
          0.5_wp*d2_rhs_v_dydy_1_0_nxp4p01m4p3p2k_nxp4p01m4p3p2m1k+&
          0.5_wp*d2_rhs_v_dydy_1_0_nxp4p01m4p3p2k_nxp4p01m4p3p2p1k

d2_rhs_v_dydy_1_0_nxp4p01m4p3p2k = d2_rhs_v_dydy_1_0_nxp4p01m4p3p2k*param_float(2)

d2_rhs_v_dydy_1_0_nxp4p01m4p3p3k_nxp4p01m4p3p3m1k = q(nx+4+0,1-4+3+3-1,indvars(3))

d2_rhs_v_dydy_1_0_nxp4p01m4p3p3k_nxp4p01m4p3p3p1k = q(nx+4+0,1-4+3+3+1,indvars(3))

d2_rhs_v_dydy_1_0_nxp4p01m4p3p3k = -&
          0.5_wp*d2_rhs_v_dydy_1_0_nxp4p01m4p3p3k_nxp4p01m4p3p3m1k+&
          0.5_wp*d2_rhs_v_dydy_1_0_nxp4p01m4p3p3k_nxp4p01m4p3p3p1k

d2_rhs_v_dydy_1_0_nxp4p01m4p3p3k = d2_rhs_v_dydy_1_0_nxp4p01m4p3p3k*param_float(2)

d1_rhs_v_dy_1_nxp4p01m4p3m1k = -((1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+4+0,1-4+3-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0,1-4+3-1,indvars(1)))**3/((q(nx+4+0,1-4+3-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0,1-4+3-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+4+0,1-4+3-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0,1-4+3-1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0,1-4+3-1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp4p01m4p3m1k)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0,1-4+3-1,indvarsst(10))*(d2_rhs_v_dydx_1_0_nxp4p01m4p3m1k)+&
                    qst(nx+4+0,1-4+3-1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp4p01m4p3m1k)))

d1_rhs_v_dy_1_nxp4p01m4p3p1k = -((1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+4+0,1-4+3+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0,1-4+3+1,indvars(1)))**3/((q(nx+4+0,1-4+3+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0,1-4+3+1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+4+0,1-4+3+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0,1-4+3+1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0,1-4+3+1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp4p01m4p3p1k)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0,1-4+3+1,indvarsst(10))*(d2_rhs_v_dydx_1_0_nxp4p01m4p3p1k)+&
                    qst(nx+4+0,1-4+3+1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp4p01m4p3p1k)))

d1_rhs_v_dy_1_nxp4p01m4p3k = -&
          0.5_wp*d1_rhs_v_dy_1_nxp4p01m4p3m1k+&
          0.5_wp*d1_rhs_v_dy_1_nxp4p01m4p3p1k

d1_rhs_v_dy_1_nxp4p01m4p3k = d1_rhs_v_dy_1_nxp4p01m4p3k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 3 None d(rho v)/dt ***********
!                                                           
!***********************************************************


rhs(nx+4+0,1-4+3,indvars(3)) =   -  ( q(nx+4+0,1-4+3,indvars(3))*(1.0_wp/(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1)))**2*((q(nx+4+0,1-4+3,indvars(2))*(-&
                    (param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1)))**2*d1_rhs_rho_dx_1_nxp4p01m4p3k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p01m4p3k))*qst(nx+4+0,1-4+3,indvarsst(10))+&
                    0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1)))+&
                    q(nx+4+0,1-4+3,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1)))*q(nx+4+0,1-4+3,indvars(1))*d1_rhs_rho_dx_3_nxp4p01m4p3k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p01m4p3k))*qst(nx+4+0,1-4+3,indvarsst(10))+&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1)))*(1-&
                    qface_i(1-4+3,2)*qface_i(1-4+3,2))/param_float(20 + 5)*((param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))-&
                    param_float(26 + 5)))))+&
                    q(nx+4+0,1-4+3,indvars(1))*((1.0_wp*q(nx+4+0,1-4+3,indvars(2))*d1_rhs_v_dx_4_nxp4p01m4p3k)*qst(nx+4+0,1-4+3,indvarsst(10)))+&
                    d1_rhs_v_dy_0_nxp4p01m4p3k+&
                    qst(nx+4+0,1-4+3,indvarsst(10))*(d1_rhs_v_dx_5_nxp4p01m4p3k)+&
                    qst(nx+4+0,1-4+3,indvarsst(11))*(d1_rhs_v_dy_1_nxp4p01m4p3k) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 3 None d(rho et)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! 0.5_wp*(u**2+v**2)*(1.0_wp/c**2*((u*(-c**2*[rho]_1x+1.0_wp*[p]_1x))*deltaxI+0.5_wp*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI+esse*c*(1-M_imax*M_imax)/L_ref*(p-P0))))+1.0_wp/gamma_m1*(0.5_wp*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI+esse*c*(1-M_imax*M_imax)/L_ref*(p-P0)))+rho*u*(1.0_wp/(2*rho*c)*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI-esse*c*(1-M_imax*M_imax)/L_ref*(p-P0)))+rho*v*((1.0_wp*u*[v]_1x)*deltaxI)+[(rho*et+p)*v]_1y+deltaxI*([-kappa*deltaxI*({T}_1x)-u*(visc_t*(2.0_wp*deltaxI*({u}_1x)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y))))-v*(visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x)))]_1x)+deltayI*([-kappa*deltayI*({T}_1y)-u*(visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x)))-v*(visc_t*(2.0_wp*deltayI*({v}_1y)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y))))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_rhs_et_dxdx_9_0_nxp4p0p01m4p3k_nxp4p0p0p01m4p3k = ((q(nx+4+0+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0+0,1-4+3,indvars(2))*q(nx+4+0+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0+0,1-4+3,indvars(3))*q(nx+4+0+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dxdx_9_0_nxp4p0p01m4p3k_nxp4p0p0m11m4p3k = ((q(nx+4+0+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0-1,1-4+3,indvars(2))*q(nx+4+0+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0+0-1,1-4+3,indvars(3))*q(nx+4+0+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dxdx_9_0_nxp4p0p01m4p3k_nxp4p0p0m21m4p3k = ((q(nx+4+0+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0-2,1-4+3,indvars(2))*q(nx+4+0+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0+0-2,1-4+3,indvars(3))*q(nx+4+0+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dxdx_9_0_nxp4p0p01m4p3k = 1.5_wp*d2_rhs_et_dxdx_9_0_nxp4p0p01m4p3k_nxp4p0p0p01m4p3k-&
          2.0_wp*d2_rhs_et_dxdx_9_0_nxp4p0p01m4p3k_nxp4p0p0m11m4p3k+&
          0.5_wp*d2_rhs_et_dxdx_9_0_nxp4p0p01m4p3k_nxp4p0p0m21m4p3k

d2_rhs_et_dxdx_9_0_nxp4p0p01m4p3k = d2_rhs_et_dxdx_9_0_nxp4p0p01m4p3k*param_float(1)

d2_rhs_et_dxdx_9_0_nxp4p0m11m4p3k_nxp4p0m1m11m4p3k = ((q(nx+4+0-1-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1-1,1-4+3,indvars(2))*q(nx+4+0-1-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1-1,1-4+3,indvars(3))*q(nx+4+0-1-1,1-4+3,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dxdx_9_0_nxp4p0m11m4p3k_nxp4p0m1p11m4p3k = ((q(nx+4+0-1+1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1+1,1-4+3,indvars(2))*q(nx+4+0-1+1,1-4+3,indvars(2))+&
                    q(nx+4+0-1+1,1-4+3,indvars(3))*q(nx+4+0-1+1,1-4+3,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dxdx_9_0_nxp4p0m11m4p3k = -&
          0.5_wp*d2_rhs_et_dxdx_9_0_nxp4p0m11m4p3k_nxp4p0m1m11m4p3k+&
          0.5_wp*d2_rhs_et_dxdx_9_0_nxp4p0m11m4p3k_nxp4p0m1p11m4p3k

d2_rhs_et_dxdx_9_0_nxp4p0m11m4p3k = d2_rhs_et_dxdx_9_0_nxp4p0m11m4p3k*param_float(1)

d2_rhs_et_dxdx_9_0_nxp4p0m21m4p3k_nxp4p0m2m11m4p3k = ((q(nx+4+0-2-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2-1,1-4+3,indvars(2))*q(nx+4+0-2-1,1-4+3,indvars(2))+&
                    q(nx+4+0-2-1,1-4+3,indvars(3))*q(nx+4+0-2-1,1-4+3,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dxdx_9_0_nxp4p0m21m4p3k_nxp4p0m2p11m4p3k = ((q(nx+4+0-2+1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2+1,1-4+3,indvars(2))*q(nx+4+0-2+1,1-4+3,indvars(2))+&
                    q(nx+4+0-2+1,1-4+3,indvars(3))*q(nx+4+0-2+1,1-4+3,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dxdx_9_0_nxp4p0m21m4p3k = -&
          0.5_wp*d2_rhs_et_dxdx_9_0_nxp4p0m21m4p3k_nxp4p0m2m11m4p3k+&
          0.5_wp*d2_rhs_et_dxdx_9_0_nxp4p0m21m4p3k_nxp4p0m2p11m4p3k

d2_rhs_et_dxdx_9_0_nxp4p0m21m4p3k = d2_rhs_et_dxdx_9_0_nxp4p0m21m4p3k*param_float(1)

d1_rhs_et_dx_9_nxp4p0p01m4p3k = -param_float(2 + 5)*qst(nx+4+0+0,1-4+3,indvarsst(10))*(d2_rhs_et_dxdx_9_0_nxp4p0p01m4p3k)-&
                    q(nx+4+0+0,1-4+3,indvars(2))*(((1+&
                    param_float(21 + 5))/(((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+4+0+0,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0+0,1-4+3,indvars(1)))**3/((q(nx+4+0+0,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0+0,1-4+3,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+4+0+0,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0+0,1-4+3,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0+0,1-4+3,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0p01m4p3k)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0+0,1-4+3,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0p01m4p3k)+&
                    qst(nx+4+0+0,1-4+3,indvarsst(11))*(d2_rhs_u_dxdy_6_0_nxp4p0p01m4p3k))))-&
                    q(nx+4+0+0,1-4+3,indvars(3))*(((1+&
                    param_float(21 + 5))/(((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+4+0+0,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0+0,1-4+3,indvars(1)))**3/((q(nx+4+0+0,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0+0,1-4+3,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+4+0+0,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0+0,1-4+3,indvars(2))*q(nx+4+0+0,1-4+3,indvars(2))+&
                    q(nx+4+0+0,1-4+3,indvars(3))*q(nx+4+0+0,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0+0,1-4+3,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0+0,1-4+3,indvarsst(11))*(d2_rhs_v_dxdy_5_0_nxp4p0p01m4p3k)+&
                    qst(nx+4+0+0,1-4+3,indvarsst(10))*(d2_rhs_v_dxdx_5_0_nxp4p0p01m4p3k)))

d1_rhs_et_dx_9_nxp4p0m11m4p3k = -param_float(2 + 5)*qst(nx+4+0-1,1-4+3,indvarsst(10))*(d2_rhs_et_dxdx_9_0_nxp4p0m11m4p3k)-&
                    q(nx+4+0-1,1-4+3,indvars(2))*(((1+&
                    param_float(21 + 5))/(((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+4+0-1,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0-1,1-4+3,indvars(1)))**3/((q(nx+4+0-1,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0-1,1-4+3,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+4+0-1,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0-1,1-4+3,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0-1,1-4+3,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0m11m4p3k)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0-1,1-4+3,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0m11m4p3k)+&
                    qst(nx+4+0-1,1-4+3,indvarsst(11))*(d2_rhs_u_dxdy_6_0_nxp4p0m11m4p3k))))-&
                    q(nx+4+0-1,1-4+3,indvars(3))*(((1+&
                    param_float(21 + 5))/(((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+4+0-1,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0-1,1-4+3,indvars(1)))**3/((q(nx+4+0-1,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0-1,1-4+3,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+4+0-1,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-1,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-1,1-4+3,indvars(2))*q(nx+4+0-1,1-4+3,indvars(2))+&
                    q(nx+4+0-1,1-4+3,indvars(3))*q(nx+4+0-1,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0-1,1-4+3,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0-1,1-4+3,indvarsst(11))*(d2_rhs_v_dxdy_5_0_nxp4p0m11m4p3k)+&
                    qst(nx+4+0-1,1-4+3,indvarsst(10))*(d2_rhs_v_dxdx_5_0_nxp4p0m11m4p3k)))

d1_rhs_et_dx_9_nxp4p0m21m4p3k = -param_float(2 + 5)*qst(nx+4+0-2,1-4+3,indvarsst(10))*(d2_rhs_et_dxdx_9_0_nxp4p0m21m4p3k)-&
                    q(nx+4+0-2,1-4+3,indvars(2))*(((1+&
                    param_float(21 + 5))/(((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+4+0-2,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0-2,1-4+3,indvars(1)))**3/((q(nx+4+0-2,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0-2,1-4+3,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+4+0-2,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0-2,1-4+3,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0-2,1-4+3,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0m21m4p3k)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0-2,1-4+3,indvarsst(10))*(d2_rhs_u_dxdx_6_0_nxp4p0m21m4p3k)+&
                    qst(nx+4+0-2,1-4+3,indvarsst(11))*(d2_rhs_u_dxdy_6_0_nxp4p0m21m4p3k))))-&
                    q(nx+4+0-2,1-4+3,indvars(3))*(((1+&
                    param_float(21 + 5))/(((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+4+0-2,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0-2,1-4+3,indvars(1)))**3/((q(nx+4+0-2,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0-2,1-4+3,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+4+0-2,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0-2,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0-2,1-4+3,indvars(2))*q(nx+4+0-2,1-4+3,indvars(2))+&
                    q(nx+4+0-2,1-4+3,indvars(3))*q(nx+4+0-2,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0-2,1-4+3,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0-2,1-4+3,indvarsst(11))*(d2_rhs_v_dxdy_5_0_nxp4p0m21m4p3k)+&
                    qst(nx+4+0-2,1-4+3,indvarsst(10))*(d2_rhs_v_dxdx_5_0_nxp4p0m21m4p3k)))

d1_rhs_et_dx_9_nxp4p01m4p3k = 1.5_wp*d1_rhs_et_dx_9_nxp4p0p01m4p3k-&
          2.0_wp*d1_rhs_et_dx_9_nxp4p0m11m4p3k+&
          0.5_wp*d1_rhs_et_dx_9_nxp4p0m21m4p3k

d1_rhs_et_dx_9_nxp4p01m4p3k = d1_rhs_et_dx_9_nxp4p01m4p3k*param_float(1)

d1_rhs_et_dy_0_nxp4p01m4p3m1k = (q(nx+4+0,1-4+3-1,indvars(1))*q(nx+4+0,1-4+3-1,indvars(4))+&
                    (param_float(3 + 5))*q(nx+4+0,1-4+3-1,indvars(1))*((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3))))))*q(nx+4+0,1-4+3-1,indvars(3))

d1_rhs_et_dy_0_nxp4p01m4p3p1k = (q(nx+4+0,1-4+3+1,indvars(1))*q(nx+4+0,1-4+3+1,indvars(4))+&
                    (param_float(3 + 5))*q(nx+4+0,1-4+3+1,indvars(1))*((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3))))))*q(nx+4+0,1-4+3+1,indvars(3))

d1_rhs_et_dy_0_nxp4p01m4p3k = -&
          0.5_wp*d1_rhs_et_dy_0_nxp4p01m4p3m1k+&
          0.5_wp*d1_rhs_et_dy_0_nxp4p01m4p3p1k

d1_rhs_et_dy_0_nxp4p01m4p3k = d1_rhs_et_dy_0_nxp4p01m4p3k*param_float(2)

d2_rhs_et_dydy_1_0_nxp4p01m4p3m3k_nxp4p01m4p3m3p0k = ((q(nx+4+0,1-4+3-3+0,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-3+0,indvars(2))*q(nx+4+0,1-4+3-3+0,indvars(2))+&
                    q(nx+4+0,1-4+3-3+0,indvars(3))*q(nx+4+0,1-4+3-3+0,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_1_0_nxp4p01m4p3m3k_nxp4p01m4p3m3p1k = ((q(nx+4+0,1-4+3-3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-3+1,indvars(2))*q(nx+4+0,1-4+3-3+1,indvars(2))+&
                    q(nx+4+0,1-4+3-3+1,indvars(3))*q(nx+4+0,1-4+3-3+1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_1_0_nxp4p01m4p3m3k_nxp4p01m4p3m3p2k = ((q(nx+4+0,1-4+3-3+2,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-3+2,indvars(2))*q(nx+4+0,1-4+3-3+2,indvars(2))+&
                    q(nx+4+0,1-4+3-3+2,indvars(3))*q(nx+4+0,1-4+3-3+2,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_1_0_nxp4p01m4p3m3k = -&
          1.5_wp*d2_rhs_et_dydy_1_0_nxp4p01m4p3m3k_nxp4p01m4p3m3p0k+&
          2.0_wp*d2_rhs_et_dydy_1_0_nxp4p01m4p3m3k_nxp4p01m4p3m3p1k-&
          0.5_wp*d2_rhs_et_dydy_1_0_nxp4p01m4p3m3k_nxp4p01m4p3m3p2k

d2_rhs_et_dydy_1_0_nxp4p01m4p3m3k = d2_rhs_et_dydy_1_0_nxp4p01m4p3m3k*param_float(2)

d2_rhs_et_dydy_1_0_nxp4p01m4p3m2k_nxp4p01m4p3m2m1k = ((q(nx+4+0,1-4+3-2-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-2-1,indvars(2))*q(nx+4+0,1-4+3-2-1,indvars(2))+&
                    q(nx+4+0,1-4+3-2-1,indvars(3))*q(nx+4+0,1-4+3-2-1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_1_0_nxp4p01m4p3m2k_nxp4p01m4p3m2p1k = ((q(nx+4+0,1-4+3-2+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-2+1,indvars(2))*q(nx+4+0,1-4+3-2+1,indvars(2))+&
                    q(nx+4+0,1-4+3-2+1,indvars(3))*q(nx+4+0,1-4+3-2+1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_1_0_nxp4p01m4p3m2k = -&
          0.5_wp*d2_rhs_et_dydy_1_0_nxp4p01m4p3m2k_nxp4p01m4p3m2m1k+&
          0.5_wp*d2_rhs_et_dydy_1_0_nxp4p01m4p3m2k_nxp4p01m4p3m2p1k

d2_rhs_et_dydy_1_0_nxp4p01m4p3m2k = d2_rhs_et_dydy_1_0_nxp4p01m4p3m2k*param_float(2)

d2_rhs_et_dydy_1_0_nxp4p01m4p3m1k_nxp4p01m4p3m1m1k = ((q(nx+4+0,1-4+3-1-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1-1,indvars(2))*q(nx+4+0,1-4+3-1-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1-1,indvars(3))*q(nx+4+0,1-4+3-1-1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_1_0_nxp4p01m4p3m1k_nxp4p01m4p3m1p1k = ((q(nx+4+0,1-4+3-1+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1+1,indvars(2))*q(nx+4+0,1-4+3-1+1,indvars(2))+&
                    q(nx+4+0,1-4+3-1+1,indvars(3))*q(nx+4+0,1-4+3-1+1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_1_0_nxp4p01m4p3m1k = -&
          0.5_wp*d2_rhs_et_dydy_1_0_nxp4p01m4p3m1k_nxp4p01m4p3m1m1k+&
          0.5_wp*d2_rhs_et_dydy_1_0_nxp4p01m4p3m1k_nxp4p01m4p3m1p1k

d2_rhs_et_dydy_1_0_nxp4p01m4p3m1k = d2_rhs_et_dydy_1_0_nxp4p01m4p3m1k*param_float(2)

d2_rhs_et_dydy_1_0_nxp4p01m4p3p1k_nxp4p01m4p3p1m1k = ((q(nx+4+0,1-4+3+1-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1-1,indvars(2))*q(nx+4+0,1-4+3+1-1,indvars(2))+&
                    q(nx+4+0,1-4+3+1-1,indvars(3))*q(nx+4+0,1-4+3+1-1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_1_0_nxp4p01m4p3p1k_nxp4p01m4p3p1p1k = ((q(nx+4+0,1-4+3+1+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1+1,indvars(2))*q(nx+4+0,1-4+3+1+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1+1,indvars(3))*q(nx+4+0,1-4+3+1+1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_1_0_nxp4p01m4p3p1k = -&
          0.5_wp*d2_rhs_et_dydy_1_0_nxp4p01m4p3p1k_nxp4p01m4p3p1m1k+&
          0.5_wp*d2_rhs_et_dydy_1_0_nxp4p01m4p3p1k_nxp4p01m4p3p1p1k

d2_rhs_et_dydy_1_0_nxp4p01m4p3p1k = d2_rhs_et_dydy_1_0_nxp4p01m4p3p1k*param_float(2)

d2_rhs_et_dydy_1_0_nxp4p01m4p3p2k_nxp4p01m4p3p2m1k = ((q(nx+4+0,1-4+3+2-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+2-1,indvars(2))*q(nx+4+0,1-4+3+2-1,indvars(2))+&
                    q(nx+4+0,1-4+3+2-1,indvars(3))*q(nx+4+0,1-4+3+2-1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_1_0_nxp4p01m4p3p2k_nxp4p01m4p3p2p1k = ((q(nx+4+0,1-4+3+2+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+2+1,indvars(2))*q(nx+4+0,1-4+3+2+1,indvars(2))+&
                    q(nx+4+0,1-4+3+2+1,indvars(3))*q(nx+4+0,1-4+3+2+1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_1_0_nxp4p01m4p3p2k = -&
          0.5_wp*d2_rhs_et_dydy_1_0_nxp4p01m4p3p2k_nxp4p01m4p3p2m1k+&
          0.5_wp*d2_rhs_et_dydy_1_0_nxp4p01m4p3p2k_nxp4p01m4p3p2p1k

d2_rhs_et_dydy_1_0_nxp4p01m4p3p2k = d2_rhs_et_dydy_1_0_nxp4p01m4p3p2k*param_float(2)

d2_rhs_et_dydy_1_0_nxp4p01m4p3p3k_nxp4p01m4p3p3m1k = ((q(nx+4+0,1-4+3+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+3-1,indvars(2))*q(nx+4+0,1-4+3+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3+3-1,indvars(3))*q(nx+4+0,1-4+3+3-1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_1_0_nxp4p01m4p3p3k_nxp4p01m4p3p3p1k = ((q(nx+4+0,1-4+3+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+3+1,indvars(2))*q(nx+4+0,1-4+3+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+3+1,indvars(3))*q(nx+4+0,1-4+3+3+1,indvars(3)))))/param_float(4 + 5)

d2_rhs_et_dydy_1_0_nxp4p01m4p3p3k = -&
          0.5_wp*d2_rhs_et_dydy_1_0_nxp4p01m4p3p3k_nxp4p01m4p3p3m1k+&
          0.5_wp*d2_rhs_et_dydy_1_0_nxp4p01m4p3p3k_nxp4p01m4p3p3p1k

d2_rhs_et_dydy_1_0_nxp4p01m4p3p3k = d2_rhs_et_dydy_1_0_nxp4p01m4p3p3k*param_float(2)

d1_rhs_et_dy_1_nxp4p01m4p3m1k = -param_float(2 + 5)*qst(nx+4+0,1-4+3-1,indvarsst(11))*(d2_rhs_et_dydy_1_0_nxp4p01m4p3m1k)-&
                    q(nx+4+0,1-4+3-1,indvars(2))*(((1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+4+0,1-4+3-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0,1-4+3-1,indvars(1)))**3/((q(nx+4+0,1-4+3-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0,1-4+3-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+4+0,1-4+3-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0,1-4+3-1,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0,1-4+3-1,indvarsst(11))*(d2_rhs_u_dydy_1_0_nxp4p01m4p3m1k)+&
                    qst(nx+4+0,1-4+3-1,indvarsst(10))*(d2_rhs_u_dydx_1_0_nxp4p01m4p3m1k)))-&
                    q(nx+4+0,1-4+3-1,indvars(3))*(((1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+4+0,1-4+3-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0,1-4+3-1,indvars(1)))**3/((q(nx+4+0,1-4+3-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0,1-4+3-1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+4+0,1-4+3-1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3-1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3-1,indvars(2))*q(nx+4+0,1-4+3-1,indvars(2))+&
                    q(nx+4+0,1-4+3-1,indvars(3))*q(nx+4+0,1-4+3-1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0,1-4+3-1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0,1-4+3-1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp4p01m4p3m1k)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0,1-4+3-1,indvarsst(10))*(d2_rhs_v_dydx_1_0_nxp4p01m4p3m1k)+&
                    qst(nx+4+0,1-4+3-1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp4p01m4p3m1k))))

d1_rhs_et_dy_1_nxp4p01m4p3p1k = -param_float(2 + 5)*qst(nx+4+0,1-4+3+1,indvarsst(11))*(d2_rhs_et_dydy_1_0_nxp4p01m4p3p1k)-&
                    q(nx+4+0,1-4+3+1,indvars(2))*(((1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+4+0,1-4+3+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0,1-4+3+1,indvars(1)))**3/((q(nx+4+0,1-4+3+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0,1-4+3+1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+4+0,1-4+3+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0,1-4+3+1,indvars(1))))*param_float(1 + 5)*(qst(nx+4+0,1-4+3+1,indvarsst(11))*(d2_rhs_u_dydy_1_0_nxp4p01m4p3p1k)+&
                    qst(nx+4+0,1-4+3+1,indvarsst(10))*(d2_rhs_u_dydx_1_0_nxp4p01m4p3p1k)))-&
                    q(nx+4+0,1-4+3+1,indvars(3))*(((1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+4+0,1-4+3+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0,1-4+3+1,indvars(1)))**3/((q(nx+4+0,1-4+3+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0,1-4+3+1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+4+0,1-4+3+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,1-4+3+1,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3+1,indvars(2))*q(nx+4+0,1-4+3+1,indvars(2))+&
                    q(nx+4+0,1-4+3+1,indvars(3))*q(nx+4+0,1-4+3+1,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0,1-4+3+1,indvars(1))))*param_float(1 + 5)*(2.0_wp*qst(nx+4+0,1-4+3+1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp4p01m4p3p1k)-&
                    2.0_wp/3.0_wp*(qst(nx+4+0,1-4+3+1,indvarsst(10))*(d2_rhs_v_dydx_1_0_nxp4p01m4p3p1k)+&
                    qst(nx+4+0,1-4+3+1,indvarsst(11))*(d2_rhs_v_dydy_1_0_nxp4p01m4p3p1k))))

d1_rhs_et_dy_1_nxp4p01m4p3k = -&
          0.5_wp*d1_rhs_et_dy_1_nxp4p01m4p3m1k+&
          0.5_wp*d1_rhs_et_dy_1_nxp4p01m4p3p1k

d1_rhs_et_dy_1_nxp4p01m4p3k = d1_rhs_et_dy_1_nxp4p01m4p3k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 3 None d(rho et)/dt **********
!                                                           
!***********************************************************


rhs(nx+4+0,1-4+3,indvars(4)) =   -  ( 0.5_wp*(q(nx+4+0,1-4+3,indvars(2))**2+&
                    q(nx+4+0,1-4+3,indvars(3))**2)*(1.0_wp/(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1)))**2*((q(nx+4+0,1-4+3,indvars(2))*(-&
                    (param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1)))**2*d1_rhs_rho_dx_1_nxp4p01m4p3k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p01m4p3k))*qst(nx+4+0,1-4+3,indvarsst(10))+&
                    0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1)))+&
                    q(nx+4+0,1-4+3,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1)))*q(nx+4+0,1-4+3,indvars(1))*d1_rhs_rho_dx_3_nxp4p01m4p3k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p01m4p3k))*qst(nx+4+0,1-4+3,indvarsst(10))+&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1)))*(1-&
                    qface_i(1-4+3,2)*qface_i(1-4+3,2))/param_float(20 + 5)*((param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))-&
                    param_float(26 + 5)))))+&
                    1.0_wp/param_float(3 + 5)*(0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1)))+&
                    q(nx+4+0,1-4+3,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1)))*q(nx+4+0,1-4+3,indvars(1))*d1_rhs_rho_dx_3_nxp4p01m4p3k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p01m4p3k))*qst(nx+4+0,1-4+3,indvarsst(10))+&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1)))*(1-&
                    qface_i(1-4+3,2)*qface_i(1-4+3,2))/param_float(20 + 5)*((param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))-&
                    param_float(26 + 5))))+&
                    q(nx+4+0,1-4+3,indvars(1))*q(nx+4+0,1-4+3,indvars(2))*(1.0_wp/(2*q(nx+4+0,1-4+3,indvars(1))*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1))))*((((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1)))+&
                    q(nx+4+0,1-4+3,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1)))*q(nx+4+0,1-4+3,indvars(1))*d1_rhs_rho_dx_3_nxp4p01m4p3k+&
                    1.0_wp*d1_rhs_rho_dx_2_nxp4p01m4p3k))*qst(nx+4+0,1-4+3,indvarsst(10))-&
                    param_float(19 + 5)*(param_float(23 + 5)*(param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))/q(nx+4+0,1-4+3,indvars(1)))*(1-&
                    qface_i(1-4+3,2)*qface_i(1-4+3,2))/param_float(20 + 5)*((param_float(3 + 5))*q(nx+4+0,1-4+3,indvars(1))*((q(nx+4+0,1-4+3,indvars(4))-&
                    0.5_wp*(q(nx+4+0,1-4+3,indvars(2))*q(nx+4+0,1-4+3,indvars(2))+&
                    q(nx+4+0,1-4+3,indvars(3))*q(nx+4+0,1-4+3,indvars(3)))))-&
                    param_float(26 + 5))))+&
                    q(nx+4+0,1-4+3,indvars(1))*q(nx+4+0,1-4+3,indvars(3))*((1.0_wp*q(nx+4+0,1-4+3,indvars(2))*d1_rhs_v_dx_4_nxp4p01m4p3k)*qst(nx+4+0,1-4+3,indvarsst(10)))+&
                    d1_rhs_et_dy_0_nxp4p01m4p3k+&
                    qst(nx+4+0,1-4+3,indvarsst(10))*(d1_rhs_et_dx_9_nxp4p01m4p3k)+&
                    qst(nx+4+0,1-4+3,indvarsst(11))*(d1_rhs_et_dy_1_nxp4p01m4p3k) ) 

