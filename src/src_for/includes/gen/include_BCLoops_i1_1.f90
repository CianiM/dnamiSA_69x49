

!***********************************************************
!                                                           
! Start building layers for BC : i1 None None **************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 1 None None ************************************
!                                                           
!***********************************************************


     do j=idloop(3),idloop(4) 


!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None d(rho)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([rho*u]_1x)+deltayI*([rho*v]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_conv_rho_dx_0_1m2p1m1jk = q(1-2+1-1,j,indvars(1))*q(1-2+1-1,j,indvars(2))

d1_conv_rho_dx_0_1m2p1p1jk = q(1-2+1+1,j,indvars(1))*q(1-2+1+1,j,indvars(2))

d1_conv_rho_dx_0_1m2p1jk = -&
          0.5_wp*d1_conv_rho_dx_0_1m2p1m1jk+&
          0.5_wp*d1_conv_rho_dx_0_1m2p1p1jk

d1_conv_rho_dx_0_1m2p1jk = d1_conv_rho_dx_0_1m2p1jk*param_float(1)

d1_conv_rho_dy_0_1m2p1jm1k = q(1-2+1,j-1,indvars(1))*q(1-2+1,j-1,indvars(3))

d1_conv_rho_dy_0_1m2p1jp1k = q(1-2+1,j+1,indvars(1))*q(1-2+1,j+1,indvars(3))

d1_conv_rho_dy_0_1m2p1jk = -&
          0.5_wp*d1_conv_rho_dy_0_1m2p1jm1k+&
          0.5_wp*d1_conv_rho_dy_0_1m2p1jp1k

d1_conv_rho_dy_0_1m2p1jk = d1_conv_rho_dy_0_1m2p1jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None d(rho)/dt **********
!                                                           
!***********************************************************


rhs(1-2+1,j,indvars(1)) =   -  ( qst(1-2+1,j,indvarsst(10))*(d1_conv_rho_dx_0_1m2p1jk)+&
                    qst(1-2+1,j,indvarsst(11))*(d1_conv_rho_dy_0_1m2p1jk) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None d(rho u)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([rho*u*u+p]_1x)+deltayI*([rho*v*u]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_conv_rhou_dx_0_1m2p1m1jk = q(1-2+1-1,j,indvars(1))*q(1-2+1-1,j,indvars(2))*q(1-2+1-1,j,indvars(2))+(param_float(3 + 5))*q(1-2+1-1,j,indvars(1))*((q(1-2+1-1,j,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,j,indvars(2))*q(1-2+1-1,j,indvars(2))+&
                    q(1-2+1-1,j,indvars(3))*q(1-2+1-1,j,indvars(3)))))

d1_conv_rhou_dx_0_1m2p1p1jk = q(1-2+1+1,j,indvars(1))*q(1-2+1+1,j,indvars(2))*q(1-2+1+1,j,indvars(2))+(param_float(3 + 5))*q(1-2+1+1,j,indvars(1))*((q(1-2+1+1,j,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,j,indvars(2))*q(1-2+1+1,j,indvars(2))+&
                    q(1-2+1+1,j,indvars(3))*q(1-2+1+1,j,indvars(3)))))

d1_conv_rhou_dx_0_1m2p1jk = -&
          0.5_wp*d1_conv_rhou_dx_0_1m2p1m1jk+&
          0.5_wp*d1_conv_rhou_dx_0_1m2p1p1jk

d1_conv_rhou_dx_0_1m2p1jk = d1_conv_rhou_dx_0_1m2p1jk*param_float(1)

d1_conv_rhou_dy_0_1m2p1jm1k = q(1-2+1,j-1,indvars(1))*q(1-2+1,j-1,indvars(3))*q(1-2+1,j-1,indvars(2))

d1_conv_rhou_dy_0_1m2p1jp1k = q(1-2+1,j+1,indvars(1))*q(1-2+1,j+1,indvars(3))*q(1-2+1,j+1,indvars(2))

d1_conv_rhou_dy_0_1m2p1jk = -&
          0.5_wp*d1_conv_rhou_dy_0_1m2p1jm1k+&
          0.5_wp*d1_conv_rhou_dy_0_1m2p1jp1k

d1_conv_rhou_dy_0_1m2p1jk = d1_conv_rhou_dy_0_1m2p1jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None d(rho u)/dt ********
!                                                           
!***********************************************************


rhs(1-2+1,j,indvars(2)) =   -  ( qst(1-2+1,j,indvarsst(10))*(d1_conv_rhou_dx_0_1m2p1jk)+&
                    qst(1-2+1,j,indvarsst(11))*(d1_conv_rhou_dy_0_1m2p1jk) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None d(rho v)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([rho*u*v]_1x)+deltayI*([rho*v*v+p]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_conv_rhov_dx_0_1m2p1m1jk = q(1-2+1-1,j,indvars(1))*q(1-2+1-1,j,indvars(2))*q(1-2+1-1,j,indvars(3))

d1_conv_rhov_dx_0_1m2p1p1jk = q(1-2+1+1,j,indvars(1))*q(1-2+1+1,j,indvars(2))*q(1-2+1+1,j,indvars(3))

d1_conv_rhov_dx_0_1m2p1jk = -&
          0.5_wp*d1_conv_rhov_dx_0_1m2p1m1jk+&
          0.5_wp*d1_conv_rhov_dx_0_1m2p1p1jk

d1_conv_rhov_dx_0_1m2p1jk = d1_conv_rhov_dx_0_1m2p1jk*param_float(1)

d1_conv_rhov_dy_0_1m2p1jm1k = q(1-2+1,j-1,indvars(1))*q(1-2+1,j-1,indvars(3))*q(1-2+1,j-1,indvars(3))+(param_float(3 + 5))*q(1-2+1,j-1,indvars(1))*((q(1-2+1,j-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,j-1,indvars(2))*q(1-2+1,j-1,indvars(2))+&
                    q(1-2+1,j-1,indvars(3))*q(1-2+1,j-1,indvars(3)))))

d1_conv_rhov_dy_0_1m2p1jp1k = q(1-2+1,j+1,indvars(1))*q(1-2+1,j+1,indvars(3))*q(1-2+1,j+1,indvars(3))+(param_float(3 + 5))*q(1-2+1,j+1,indvars(1))*((q(1-2+1,j+1,indvars(4))-&
                    0.5_wp*(q(1-2+1,j+1,indvars(2))*q(1-2+1,j+1,indvars(2))+&
                    q(1-2+1,j+1,indvars(3))*q(1-2+1,j+1,indvars(3)))))

d1_conv_rhov_dy_0_1m2p1jk = -&
          0.5_wp*d1_conv_rhov_dy_0_1m2p1jm1k+&
          0.5_wp*d1_conv_rhov_dy_0_1m2p1jp1k

d1_conv_rhov_dy_0_1m2p1jk = d1_conv_rhov_dy_0_1m2p1jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None d(rho v)/dt ********
!                                                           
!***********************************************************


rhs(1-2+1,j,indvars(3)) =   -  ( qst(1-2+1,j,indvarsst(10))*(d1_conv_rhov_dx_0_1m2p1jk)+&
                    qst(1-2+1,j,indvarsst(11))*(d1_conv_rhov_dy_0_1m2p1jk) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None d(rho et)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([(rho*et+p)*u]_1x)+deltayI*([(rho*et+p)*v]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_conv_rhoet_dx_0_1m2p1m1jk = (q(1-2+1-1,j,indvars(1))*q(1-2+1-1,j,indvars(4))+&
                    (param_float(3 + 5))*q(1-2+1-1,j,indvars(1))*((q(1-2+1-1,j,indvars(4))-&
                    0.5_wp*(q(1-2+1-1,j,indvars(2))*q(1-2+1-1,j,indvars(2))+&
                    q(1-2+1-1,j,indvars(3))*q(1-2+1-1,j,indvars(3))))))*q(1-2+1-1,j,indvars(2))

d1_conv_rhoet_dx_0_1m2p1p1jk = (q(1-2+1+1,j,indvars(1))*q(1-2+1+1,j,indvars(4))+&
                    (param_float(3 + 5))*q(1-2+1+1,j,indvars(1))*((q(1-2+1+1,j,indvars(4))-&
                    0.5_wp*(q(1-2+1+1,j,indvars(2))*q(1-2+1+1,j,indvars(2))+&
                    q(1-2+1+1,j,indvars(3))*q(1-2+1+1,j,indvars(3))))))*q(1-2+1+1,j,indvars(2))

d1_conv_rhoet_dx_0_1m2p1jk = -&
          0.5_wp*d1_conv_rhoet_dx_0_1m2p1m1jk+&
          0.5_wp*d1_conv_rhoet_dx_0_1m2p1p1jk

d1_conv_rhoet_dx_0_1m2p1jk = d1_conv_rhoet_dx_0_1m2p1jk*param_float(1)

d1_conv_rhoet_dy_0_1m2p1jm1k = (q(1-2+1,j-1,indvars(1))*q(1-2+1,j-1,indvars(4))+&
                    (param_float(3 + 5))*q(1-2+1,j-1,indvars(1))*((q(1-2+1,j-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,j-1,indvars(2))*q(1-2+1,j-1,indvars(2))+&
                    q(1-2+1,j-1,indvars(3))*q(1-2+1,j-1,indvars(3))))))*q(1-2+1,j-1,indvars(3))

d1_conv_rhoet_dy_0_1m2p1jp1k = (q(1-2+1,j+1,indvars(1))*q(1-2+1,j+1,indvars(4))+&
                    (param_float(3 + 5))*q(1-2+1,j+1,indvars(1))*((q(1-2+1,j+1,indvars(4))-&
                    0.5_wp*(q(1-2+1,j+1,indvars(2))*q(1-2+1,j+1,indvars(2))+&
                    q(1-2+1,j+1,indvars(3))*q(1-2+1,j+1,indvars(3))))))*q(1-2+1,j+1,indvars(3))

d1_conv_rhoet_dy_0_1m2p1jk = -&
          0.5_wp*d1_conv_rhoet_dy_0_1m2p1jm1k+&
          0.5_wp*d1_conv_rhoet_dy_0_1m2p1jp1k

d1_conv_rhoet_dy_0_1m2p1jk = d1_conv_rhoet_dy_0_1m2p1jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None d(rho et)/dt *******
!                                                           
!***********************************************************


rhs(1-2+1,j,indvars(4)) =   -  ( qst(1-2+1,j,indvarsst(10))*(d1_conv_rhoet_dx_0_1m2p1jk)+&
                    qst(1-2+1,j,indvarsst(11))*(d1_conv_rhoet_dy_0_1m2p1jk) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None d(rho nut)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([(rho*nut)*u]_1x)+deltayI*([(rho*nut)*v]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_conv_rhonut_dx_0_1m2p1m1jk = (q(1-2+1-1,j,indvars(1))*q(1-2+1-1,j,indvars(5)))*q(1-2+1-1,j,indvars(2))

d1_conv_rhonut_dx_0_1m2p1p1jk = (q(1-2+1+1,j,indvars(1))*q(1-2+1+1,j,indvars(5)))*q(1-2+1+1,j,indvars(2))

d1_conv_rhonut_dx_0_1m2p1jk = -&
          0.5_wp*d1_conv_rhonut_dx_0_1m2p1m1jk+&
          0.5_wp*d1_conv_rhonut_dx_0_1m2p1p1jk

d1_conv_rhonut_dx_0_1m2p1jk = d1_conv_rhonut_dx_0_1m2p1jk*param_float(1)

d1_conv_rhonut_dy_0_1m2p1jm1k = (q(1-2+1,j-1,indvars(1))*q(1-2+1,j-1,indvars(5)))*q(1-2+1,j-1,indvars(3))

d1_conv_rhonut_dy_0_1m2p1jp1k = (q(1-2+1,j+1,indvars(1))*q(1-2+1,j+1,indvars(5)))*q(1-2+1,j+1,indvars(3))

d1_conv_rhonut_dy_0_1m2p1jk = -&
          0.5_wp*d1_conv_rhonut_dy_0_1m2p1jm1k+&
          0.5_wp*d1_conv_rhonut_dy_0_1m2p1jp1k

d1_conv_rhonut_dy_0_1m2p1jk = d1_conv_rhonut_dy_0_1m2p1jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None d(rho nut)/dt ******
!                                                           
!***********************************************************


rhs(1-2+1,j,indvars(5)) =   -  ( qst(1-2+1,j,indvarsst(10))*(d1_conv_rhonut_dx_0_1m2p1jk)+&
                    qst(1-2+1,j,indvarsst(11))*(d1_conv_rhonut_dy_0_1m2p1jk) ) 

     enddo


!***********************************************************
!                                                           
! Start building layers for BC : i1 None None **************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 1 None None ************************************
!                                                           
!***********************************************************


     do j=idloop(3),idloop(4) 


!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None d(rho u)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([-visc_t*(2.0_wp*deltaxI*({u}_1x)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y)))]_1x)+deltayI*([-visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_dif_rhou_dxdx_0_0_1m2p1m1jk_1m2p1m1p0jk = q(1-2+1-1+0,j,indvars(2))

d2_dif_rhou_dxdx_0_0_1m2p1m1jk_1m2p1m1p1jk = q(1-2+1-1+1,j,indvars(2))

d2_dif_rhou_dxdx_0_0_1m2p1m1jk_1m2p1m1p2jk = q(1-2+1-1+2,j,indvars(2))

d2_dif_rhou_dxdx_0_0_1m2p1m1jk = -&
          1.5_wp*d2_dif_rhou_dxdx_0_0_1m2p1m1jk_1m2p1m1p0jk+&
          2.0_wp*d2_dif_rhou_dxdx_0_0_1m2p1m1jk_1m2p1m1p1jk-&
          0.5_wp*d2_dif_rhou_dxdx_0_0_1m2p1m1jk_1m2p1m1p2jk

d2_dif_rhou_dxdx_0_0_1m2p1m1jk = d2_dif_rhou_dxdx_0_0_1m2p1m1jk*param_float(1)

d2_dif_rhou_dxdx_0_0_1m2p1p1jk_1m2p1p1m1jk = q(1-2+1+1-1,j,indvars(2))

d2_dif_rhou_dxdx_0_0_1m2p1p1jk_1m2p1p1p1jk = q(1-2+1+1+1,j,indvars(2))

d2_dif_rhou_dxdx_0_0_1m2p1p1jk = -&
          0.5_wp*d2_dif_rhou_dxdx_0_0_1m2p1p1jk_1m2p1p1m1jk+&
          0.5_wp*d2_dif_rhou_dxdx_0_0_1m2p1p1jk_1m2p1p1p1jk

d2_dif_rhou_dxdx_0_0_1m2p1p1jk = d2_dif_rhou_dxdx_0_0_1m2p1p1jk*param_float(1)

d2_dif_rhou_dxdy_0_0_1m2p1m1jk_1m2p1m1jm1k = q(1-2+1-1,j-1,indvars(3))

d2_dif_rhou_dxdy_0_0_1m2p1m1jk_1m2p1m1jp1k = q(1-2+1-1,j+1,indvars(3))

d2_dif_rhou_dxdy_0_0_1m2p1m1jk = -&
          0.5_wp*d2_dif_rhou_dxdy_0_0_1m2p1m1jk_1m2p1m1jm1k+&
          0.5_wp*d2_dif_rhou_dxdy_0_0_1m2p1m1jk_1m2p1m1jp1k

d2_dif_rhou_dxdy_0_0_1m2p1m1jk = d2_dif_rhou_dxdy_0_0_1m2p1m1jk*param_float(2)

d2_dif_rhou_dxdy_0_0_1m2p1p1jk_1m2p1p1jm1k = q(1-2+1+1,j-1,indvars(3))

d2_dif_rhou_dxdy_0_0_1m2p1p1jk_1m2p1p1jp1k = q(1-2+1+1,j+1,indvars(3))

d2_dif_rhou_dxdy_0_0_1m2p1p1jk = -&
          0.5_wp*d2_dif_rhou_dxdy_0_0_1m2p1p1jk_1m2p1p1jm1k+&
          0.5_wp*d2_dif_rhou_dxdy_0_0_1m2p1p1jk_1m2p1p1jp1k

d2_dif_rhou_dxdy_0_0_1m2p1p1jk = d2_dif_rhou_dxdy_0_0_1m2p1p1jk*param_float(2)

d1_dif_rhou_dx_0_1m2p1m1jk = -(1.0_wp)*(1.0_wp+&
                    ((q(1-2+1-1,j,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1-1,j,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1-1,j,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(1-2+1-1,j,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_1m2p1m1jk)-&
                    2.0_wp/3.0_wp*(qst(1-2+1-1,j,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_1m2p1m1jk)+&
                    qst(1-2+1-1,j,indvarsst(11))*(d2_dif_rhou_dxdy_0_0_1m2p1m1jk)))

d1_dif_rhou_dx_0_1m2p1p1jk = -(1.0_wp)*(1.0_wp+&
                    ((q(1-2+1+1,j,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1+1,j,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1+1,j,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(1-2+1+1,j,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_1m2p1p1jk)-&
                    2.0_wp/3.0_wp*(qst(1-2+1+1,j,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_1m2p1p1jk)+&
                    qst(1-2+1+1,j,indvarsst(11))*(d2_dif_rhou_dxdy_0_0_1m2p1p1jk)))

d1_dif_rhou_dx_0_1m2p1jk = -&
          0.5_wp*d1_dif_rhou_dx_0_1m2p1m1jk+&
          0.5_wp*d1_dif_rhou_dx_0_1m2p1p1jk

d1_dif_rhou_dx_0_1m2p1jk = d1_dif_rhou_dx_0_1m2p1jk*param_float(1)

d2_dif_rhou_dydx_0_0_1m2p1jm1k_1m2p1m1jm1k = q(1-2+1-1,j-1,indvars(3))

d2_dif_rhou_dydx_0_0_1m2p1jm1k_1m2p1p1jm1k = q(1-2+1+1,j-1,indvars(3))

d2_dif_rhou_dydx_0_0_1m2p1jm1k = -&
          0.5_wp*d2_dif_rhou_dydx_0_0_1m2p1jm1k_1m2p1m1jm1k+&
          0.5_wp*d2_dif_rhou_dydx_0_0_1m2p1jm1k_1m2p1p1jm1k

d2_dif_rhou_dydx_0_0_1m2p1jm1k = d2_dif_rhou_dydx_0_0_1m2p1jm1k*param_float(1)

d2_dif_rhou_dydx_0_0_1m2p1jp1k_1m2p1m1jp1k = q(1-2+1-1,j+1,indvars(3))

d2_dif_rhou_dydx_0_0_1m2p1jp1k_1m2p1p1jp1k = q(1-2+1+1,j+1,indvars(3))

d2_dif_rhou_dydx_0_0_1m2p1jp1k = -&
          0.5_wp*d2_dif_rhou_dydx_0_0_1m2p1jp1k_1m2p1m1jp1k+&
          0.5_wp*d2_dif_rhou_dydx_0_0_1m2p1jp1k_1m2p1p1jp1k

d2_dif_rhou_dydx_0_0_1m2p1jp1k = d2_dif_rhou_dydx_0_0_1m2p1jp1k*param_float(1)

d2_dif_rhou_dydy_0_0_1m2p1jm1k_1m2p1jm1m1k = q(1-2+1,j-1-1,indvars(2))

d2_dif_rhou_dydy_0_0_1m2p1jm1k_1m2p1jm1p1k = q(1-2+1,j-1+1,indvars(2))

d2_dif_rhou_dydy_0_0_1m2p1jm1k = -&
          0.5_wp*d2_dif_rhou_dydy_0_0_1m2p1jm1k_1m2p1jm1m1k+&
          0.5_wp*d2_dif_rhou_dydy_0_0_1m2p1jm1k_1m2p1jm1p1k

d2_dif_rhou_dydy_0_0_1m2p1jm1k = d2_dif_rhou_dydy_0_0_1m2p1jm1k*param_float(2)

d2_dif_rhou_dydy_0_0_1m2p1jp1k_1m2p1jp1m1k = q(1-2+1,j+1-1,indvars(2))

d2_dif_rhou_dydy_0_0_1m2p1jp1k_1m2p1jp1p1k = q(1-2+1,j+1+1,indvars(2))

d2_dif_rhou_dydy_0_0_1m2p1jp1k = -&
          0.5_wp*d2_dif_rhou_dydy_0_0_1m2p1jp1k_1m2p1jp1m1k+&
          0.5_wp*d2_dif_rhou_dydy_0_0_1m2p1jp1k_1m2p1jp1p1k

d2_dif_rhou_dydy_0_0_1m2p1jp1k = d2_dif_rhou_dydy_0_0_1m2p1jp1k*param_float(2)

d1_dif_rhou_dy_0_1m2p1jm1k = -(1.0_wp)*(1.0_wp+&
                    ((q(1-2+1,j-1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1,j-1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1,j-1,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(1-2+1,j-1,indvarsst(11))*(d2_dif_rhou_dydy_0_0_1m2p1jm1k)+&
                    qst(1-2+1,j-1,indvarsst(10))*(d2_dif_rhou_dydx_0_0_1m2p1jm1k))

d1_dif_rhou_dy_0_1m2p1jp1k = -(1.0_wp)*(1.0_wp+&
                    ((q(1-2+1,j+1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1,j+1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1,j+1,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(1-2+1,j+1,indvarsst(11))*(d2_dif_rhou_dydy_0_0_1m2p1jp1k)+&
                    qst(1-2+1,j+1,indvarsst(10))*(d2_dif_rhou_dydx_0_0_1m2p1jp1k))

d1_dif_rhou_dy_0_1m2p1jk = -&
          0.5_wp*d1_dif_rhou_dy_0_1m2p1jm1k+&
          0.5_wp*d1_dif_rhou_dy_0_1m2p1jp1k

d1_dif_rhou_dy_0_1m2p1jk = d1_dif_rhou_dy_0_1m2p1jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None d(rho u)/dt ********
!                                                           
!***********************************************************


rhs(1-2+1,j,indvars(2)) = rhs(1-2+1,j,indvars(2))  -  ( qst(1-2+1,j,indvarsst(10))*(d1_dif_rhou_dx_0_1m2p1jk)+&
                    qst(1-2+1,j,indvarsst(11))*(d1_dif_rhou_dy_0_1m2p1jk) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None d(rho v)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([-visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x))]_1x)+deltayI*([-visc_t*(2.0_wp*deltayI*({v}_1y)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y)))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_dif_rhov_dxdx_0_0_1m2p1m1jk_1m2p1m1p0jk = q(1-2+1-1+0,j,indvars(3))

d2_dif_rhov_dxdx_0_0_1m2p1m1jk_1m2p1m1p1jk = q(1-2+1-1+1,j,indvars(3))

d2_dif_rhov_dxdx_0_0_1m2p1m1jk_1m2p1m1p2jk = q(1-2+1-1+2,j,indvars(3))

d2_dif_rhov_dxdx_0_0_1m2p1m1jk = -&
          1.5_wp*d2_dif_rhov_dxdx_0_0_1m2p1m1jk_1m2p1m1p0jk+&
          2.0_wp*d2_dif_rhov_dxdx_0_0_1m2p1m1jk_1m2p1m1p1jk-&
          0.5_wp*d2_dif_rhov_dxdx_0_0_1m2p1m1jk_1m2p1m1p2jk

d2_dif_rhov_dxdx_0_0_1m2p1m1jk = d2_dif_rhov_dxdx_0_0_1m2p1m1jk*param_float(1)

d2_dif_rhov_dxdx_0_0_1m2p1p1jk_1m2p1p1m1jk = q(1-2+1+1-1,j,indvars(3))

d2_dif_rhov_dxdx_0_0_1m2p1p1jk_1m2p1p1p1jk = q(1-2+1+1+1,j,indvars(3))

d2_dif_rhov_dxdx_0_0_1m2p1p1jk = -&
          0.5_wp*d2_dif_rhov_dxdx_0_0_1m2p1p1jk_1m2p1p1m1jk+&
          0.5_wp*d2_dif_rhov_dxdx_0_0_1m2p1p1jk_1m2p1p1p1jk

d2_dif_rhov_dxdx_0_0_1m2p1p1jk = d2_dif_rhov_dxdx_0_0_1m2p1p1jk*param_float(1)

d2_dif_rhov_dxdy_0_0_1m2p1m1jk_1m2p1m1jm1k = q(1-2+1-1,j-1,indvars(2))

d2_dif_rhov_dxdy_0_0_1m2p1m1jk_1m2p1m1jp1k = q(1-2+1-1,j+1,indvars(2))

d2_dif_rhov_dxdy_0_0_1m2p1m1jk = -&
          0.5_wp*d2_dif_rhov_dxdy_0_0_1m2p1m1jk_1m2p1m1jm1k+&
          0.5_wp*d2_dif_rhov_dxdy_0_0_1m2p1m1jk_1m2p1m1jp1k

d2_dif_rhov_dxdy_0_0_1m2p1m1jk = d2_dif_rhov_dxdy_0_0_1m2p1m1jk*param_float(2)

d2_dif_rhov_dxdy_0_0_1m2p1p1jk_1m2p1p1jm1k = q(1-2+1+1,j-1,indvars(2))

d2_dif_rhov_dxdy_0_0_1m2p1p1jk_1m2p1p1jp1k = q(1-2+1+1,j+1,indvars(2))

d2_dif_rhov_dxdy_0_0_1m2p1p1jk = -&
          0.5_wp*d2_dif_rhov_dxdy_0_0_1m2p1p1jk_1m2p1p1jm1k+&
          0.5_wp*d2_dif_rhov_dxdy_0_0_1m2p1p1jk_1m2p1p1jp1k

d2_dif_rhov_dxdy_0_0_1m2p1p1jk = d2_dif_rhov_dxdy_0_0_1m2p1p1jk*param_float(2)

d1_dif_rhov_dx_0_1m2p1m1jk = -(1.0_wp)*(1.0_wp+&
                    ((q(1-2+1-1,j,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1-1,j,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1-1,j,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(1-2+1-1,j,indvarsst(11))*(d2_dif_rhov_dxdy_0_0_1m2p1m1jk)+&
                    qst(1-2+1-1,j,indvarsst(10))*(d2_dif_rhov_dxdx_0_0_1m2p1m1jk))

d1_dif_rhov_dx_0_1m2p1p1jk = -(1.0_wp)*(1.0_wp+&
                    ((q(1-2+1+1,j,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1+1,j,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1+1,j,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(1-2+1+1,j,indvarsst(11))*(d2_dif_rhov_dxdy_0_0_1m2p1p1jk)+&
                    qst(1-2+1+1,j,indvarsst(10))*(d2_dif_rhov_dxdx_0_0_1m2p1p1jk))

d1_dif_rhov_dx_0_1m2p1jk = -&
          0.5_wp*d1_dif_rhov_dx_0_1m2p1m1jk+&
          0.5_wp*d1_dif_rhov_dx_0_1m2p1p1jk

d1_dif_rhov_dx_0_1m2p1jk = d1_dif_rhov_dx_0_1m2p1jk*param_float(1)

d2_dif_rhov_dydx_0_0_1m2p1jm1k_1m2p1m1jm1k = q(1-2+1-1,j-1,indvars(2))

d2_dif_rhov_dydx_0_0_1m2p1jm1k_1m2p1p1jm1k = q(1-2+1+1,j-1,indvars(2))

d2_dif_rhov_dydx_0_0_1m2p1jm1k = -&
          0.5_wp*d2_dif_rhov_dydx_0_0_1m2p1jm1k_1m2p1m1jm1k+&
          0.5_wp*d2_dif_rhov_dydx_0_0_1m2p1jm1k_1m2p1p1jm1k

d2_dif_rhov_dydx_0_0_1m2p1jm1k = d2_dif_rhov_dydx_0_0_1m2p1jm1k*param_float(1)

d2_dif_rhov_dydx_0_0_1m2p1jp1k_1m2p1m1jp1k = q(1-2+1-1,j+1,indvars(2))

d2_dif_rhov_dydx_0_0_1m2p1jp1k_1m2p1p1jp1k = q(1-2+1+1,j+1,indvars(2))

d2_dif_rhov_dydx_0_0_1m2p1jp1k = -&
          0.5_wp*d2_dif_rhov_dydx_0_0_1m2p1jp1k_1m2p1m1jp1k+&
          0.5_wp*d2_dif_rhov_dydx_0_0_1m2p1jp1k_1m2p1p1jp1k

d2_dif_rhov_dydx_0_0_1m2p1jp1k = d2_dif_rhov_dydx_0_0_1m2p1jp1k*param_float(1)

d2_dif_rhov_dydy_0_0_1m2p1jm1k_1m2p1jm1m1k = q(1-2+1,j-1-1,indvars(3))

d2_dif_rhov_dydy_0_0_1m2p1jm1k_1m2p1jm1p1k = q(1-2+1,j-1+1,indvars(3))

d2_dif_rhov_dydy_0_0_1m2p1jm1k = -&
          0.5_wp*d2_dif_rhov_dydy_0_0_1m2p1jm1k_1m2p1jm1m1k+&
          0.5_wp*d2_dif_rhov_dydy_0_0_1m2p1jm1k_1m2p1jm1p1k

d2_dif_rhov_dydy_0_0_1m2p1jm1k = d2_dif_rhov_dydy_0_0_1m2p1jm1k*param_float(2)

d2_dif_rhov_dydy_0_0_1m2p1jp1k_1m2p1jp1m1k = q(1-2+1,j+1-1,indvars(3))

d2_dif_rhov_dydy_0_0_1m2p1jp1k_1m2p1jp1p1k = q(1-2+1,j+1+1,indvars(3))

d2_dif_rhov_dydy_0_0_1m2p1jp1k = -&
          0.5_wp*d2_dif_rhov_dydy_0_0_1m2p1jp1k_1m2p1jp1m1k+&
          0.5_wp*d2_dif_rhov_dydy_0_0_1m2p1jp1k_1m2p1jp1p1k

d2_dif_rhov_dydy_0_0_1m2p1jp1k = d2_dif_rhov_dydy_0_0_1m2p1jp1k*param_float(2)

d1_dif_rhov_dy_0_1m2p1jm1k = -(1.0_wp)*(1.0_wp+&
                    ((q(1-2+1,j-1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1,j-1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1,j-1,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(1-2+1,j-1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_1m2p1jm1k)-&
                    2.0_wp/3.0_wp*(qst(1-2+1,j-1,indvarsst(10))*(d2_dif_rhov_dydx_0_0_1m2p1jm1k)+&
                    qst(1-2+1,j-1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_1m2p1jm1k)))

d1_dif_rhov_dy_0_1m2p1jp1k = -(1.0_wp)*(1.0_wp+&
                    ((q(1-2+1,j+1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1,j+1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1,j+1,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(1-2+1,j+1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_1m2p1jp1k)-&
                    2.0_wp/3.0_wp*(qst(1-2+1,j+1,indvarsst(10))*(d2_dif_rhov_dydx_0_0_1m2p1jp1k)+&
                    qst(1-2+1,j+1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_1m2p1jp1k)))

d1_dif_rhov_dy_0_1m2p1jk = -&
          0.5_wp*d1_dif_rhov_dy_0_1m2p1jm1k+&
          0.5_wp*d1_dif_rhov_dy_0_1m2p1jp1k

d1_dif_rhov_dy_0_1m2p1jk = d1_dif_rhov_dy_0_1m2p1jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None d(rho v)/dt ********
!                                                           
!***********************************************************


rhs(1-2+1,j,indvars(3)) = rhs(1-2+1,j,indvars(3))  -  ( qst(1-2+1,j,indvarsst(10))*(d1_dif_rhov_dx_0_1m2p1jk)+&
                    qst(1-2+1,j,indvarsst(11))*(d1_dif_rhov_dy_0_1m2p1jk) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None d(rho et)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([-kappa*deltaxI*({T}_1x)-u*(visc_t*(2.0_wp*deltaxI*({u}_1x)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y))))-v*(visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x)))]_1x)+deltayI*([-kappa*deltayI*({T}_1y)-u*(visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x)))-v*(visc_t*(2.0_wp*deltayI*({v}_1y)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y))))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_dif_rhoet_dxdx_0_0_1m2p1m1jk_1m2p1m1p0jk = ((q(1-2+1-1+0,j,indvars(4))-&
                    0.5_wp*(q(1-2+1-1+0,j,indvars(2))*q(1-2+1-1+0,j,indvars(2))+&
                    q(1-2+1-1+0,j,indvars(3))*q(1-2+1-1+0,j,indvars(3)))))/(param_float(4 + 5))

d2_dif_rhoet_dxdx_0_0_1m2p1m1jk_1m2p1m1p1jk = ((q(1-2+1-1+1,j,indvars(4))-&
                    0.5_wp*(q(1-2+1-1+1,j,indvars(2))*q(1-2+1-1+1,j,indvars(2))+&
                    q(1-2+1-1+1,j,indvars(3))*q(1-2+1-1+1,j,indvars(3)))))/(param_float(4 + 5))

d2_dif_rhoet_dxdx_0_0_1m2p1m1jk_1m2p1m1p2jk = ((q(1-2+1-1+2,j,indvars(4))-&
                    0.5_wp*(q(1-2+1-1+2,j,indvars(2))*q(1-2+1-1+2,j,indvars(2))+&
                    q(1-2+1-1+2,j,indvars(3))*q(1-2+1-1+2,j,indvars(3)))))/(param_float(4 + 5))

d2_dif_rhoet_dxdx_0_0_1m2p1m1jk = -&
          1.5_wp*d2_dif_rhoet_dxdx_0_0_1m2p1m1jk_1m2p1m1p0jk+&
          2.0_wp*d2_dif_rhoet_dxdx_0_0_1m2p1m1jk_1m2p1m1p1jk-&
          0.5_wp*d2_dif_rhoet_dxdx_0_0_1m2p1m1jk_1m2p1m1p2jk

d2_dif_rhoet_dxdx_0_0_1m2p1m1jk = d2_dif_rhoet_dxdx_0_0_1m2p1m1jk*param_float(1)

d2_dif_rhoet_dxdx_0_0_1m2p1p1jk_1m2p1p1m1jk = ((q(1-2+1+1-1,j,indvars(4))-&
                    0.5_wp*(q(1-2+1+1-1,j,indvars(2))*q(1-2+1+1-1,j,indvars(2))+&
                    q(1-2+1+1-1,j,indvars(3))*q(1-2+1+1-1,j,indvars(3)))))/(param_float(4 + 5))

d2_dif_rhoet_dxdx_0_0_1m2p1p1jk_1m2p1p1p1jk = ((q(1-2+1+1+1,j,indvars(4))-&
                    0.5_wp*(q(1-2+1+1+1,j,indvars(2))*q(1-2+1+1+1,j,indvars(2))+&
                    q(1-2+1+1+1,j,indvars(3))*q(1-2+1+1+1,j,indvars(3)))))/(param_float(4 + 5))

d2_dif_rhoet_dxdx_0_0_1m2p1p1jk = -&
          0.5_wp*d2_dif_rhoet_dxdx_0_0_1m2p1p1jk_1m2p1p1m1jk+&
          0.5_wp*d2_dif_rhoet_dxdx_0_0_1m2p1p1jk_1m2p1p1p1jk

d2_dif_rhoet_dxdx_0_0_1m2p1p1jk = d2_dif_rhoet_dxdx_0_0_1m2p1p1jk*param_float(1)

d1_dif_rhoet_dx_0_1m2p1m1jk = -param_float(2 + 5)*qst(1-2+1-1,j,indvarsst(10))*(d2_dif_rhoet_dxdx_0_0_1m2p1m1jk)-&
                    q(1-2+1-1,j,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(1-2+1-1,j,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1-1,j,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1-1,j,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(1-2+1-1,j,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_1m2p1m1jk)-&
                    2.0_wp/3.0_wp*(qst(1-2+1-1,j,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_1m2p1m1jk)+&
                    qst(1-2+1-1,j,indvarsst(11))*(d2_dif_rhou_dxdy_0_0_1m2p1m1jk))))-&
                    q(1-2+1-1,j,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(1-2+1-1,j,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1-1,j,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1-1,j,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(1-2+1-1,j,indvarsst(11))*(d2_dif_rhov_dxdy_0_0_1m2p1m1jk)+&
                    qst(1-2+1-1,j,indvarsst(10))*(d2_dif_rhov_dxdx_0_0_1m2p1m1jk)))

d1_dif_rhoet_dx_0_1m2p1p1jk = -param_float(2 + 5)*qst(1-2+1+1,j,indvarsst(10))*(d2_dif_rhoet_dxdx_0_0_1m2p1p1jk)-&
                    q(1-2+1+1,j,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(1-2+1+1,j,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1+1,j,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1+1,j,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(1-2+1+1,j,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_1m2p1p1jk)-&
                    2.0_wp/3.0_wp*(qst(1-2+1+1,j,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_1m2p1p1jk)+&
                    qst(1-2+1+1,j,indvarsst(11))*(d2_dif_rhou_dxdy_0_0_1m2p1p1jk))))-&
                    q(1-2+1+1,j,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(1-2+1+1,j,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1+1,j,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1+1,j,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(1-2+1+1,j,indvarsst(11))*(d2_dif_rhov_dxdy_0_0_1m2p1p1jk)+&
                    qst(1-2+1+1,j,indvarsst(10))*(d2_dif_rhov_dxdx_0_0_1m2p1p1jk)))

d1_dif_rhoet_dx_0_1m2p1jk = -&
          0.5_wp*d1_dif_rhoet_dx_0_1m2p1m1jk+&
          0.5_wp*d1_dif_rhoet_dx_0_1m2p1p1jk

d1_dif_rhoet_dx_0_1m2p1jk = d1_dif_rhoet_dx_0_1m2p1jk*param_float(1)

d2_dif_rhoet_dydy_0_0_1m2p1jm1k_1m2p1jm1m1k = ((q(1-2+1,j-1-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,j-1-1,indvars(2))*q(1-2+1,j-1-1,indvars(2))+&
                    q(1-2+1,j-1-1,indvars(3))*q(1-2+1,j-1-1,indvars(3)))))/(param_float(4 + 5))

d2_dif_rhoet_dydy_0_0_1m2p1jm1k_1m2p1jm1p1k = ((q(1-2+1,j-1+1,indvars(4))-&
                    0.5_wp*(q(1-2+1,j-1+1,indvars(2))*q(1-2+1,j-1+1,indvars(2))+&
                    q(1-2+1,j-1+1,indvars(3))*q(1-2+1,j-1+1,indvars(3)))))/(param_float(4 + 5))

d2_dif_rhoet_dydy_0_0_1m2p1jm1k = -&
          0.5_wp*d2_dif_rhoet_dydy_0_0_1m2p1jm1k_1m2p1jm1m1k+&
          0.5_wp*d2_dif_rhoet_dydy_0_0_1m2p1jm1k_1m2p1jm1p1k

d2_dif_rhoet_dydy_0_0_1m2p1jm1k = d2_dif_rhoet_dydy_0_0_1m2p1jm1k*param_float(2)

d2_dif_rhoet_dydy_0_0_1m2p1jp1k_1m2p1jp1m1k = ((q(1-2+1,j+1-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,j+1-1,indvars(2))*q(1-2+1,j+1-1,indvars(2))+&
                    q(1-2+1,j+1-1,indvars(3))*q(1-2+1,j+1-1,indvars(3)))))/(param_float(4 + 5))

d2_dif_rhoet_dydy_0_0_1m2p1jp1k_1m2p1jp1p1k = ((q(1-2+1,j+1+1,indvars(4))-&
                    0.5_wp*(q(1-2+1,j+1+1,indvars(2))*q(1-2+1,j+1+1,indvars(2))+&
                    q(1-2+1,j+1+1,indvars(3))*q(1-2+1,j+1+1,indvars(3)))))/(param_float(4 + 5))

d2_dif_rhoet_dydy_0_0_1m2p1jp1k = -&
          0.5_wp*d2_dif_rhoet_dydy_0_0_1m2p1jp1k_1m2p1jp1m1k+&
          0.5_wp*d2_dif_rhoet_dydy_0_0_1m2p1jp1k_1m2p1jp1p1k

d2_dif_rhoet_dydy_0_0_1m2p1jp1k = d2_dif_rhoet_dydy_0_0_1m2p1jp1k*param_float(2)

d1_dif_rhoet_dy_0_1m2p1jm1k = -param_float(2 + 5)*qst(1-2+1,j-1,indvarsst(11))*(d2_dif_rhoet_dydy_0_0_1m2p1jm1k)-&
                    q(1-2+1,j-1,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(1-2+1,j-1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1,j-1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1,j-1,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(1-2+1,j-1,indvarsst(11))*(d2_dif_rhou_dydy_0_0_1m2p1jm1k)+&
                    qst(1-2+1,j-1,indvarsst(10))*(d2_dif_rhou_dydx_0_0_1m2p1jm1k)))-&
                    q(1-2+1,j-1,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(1-2+1,j-1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1,j-1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1,j-1,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(1-2+1,j-1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_1m2p1jm1k)-&
                    2.0_wp/3.0_wp*(qst(1-2+1,j-1,indvarsst(10))*(d2_dif_rhov_dydx_0_0_1m2p1jm1k)+&
                    qst(1-2+1,j-1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_1m2p1jm1k))))

d1_dif_rhoet_dy_0_1m2p1jp1k = -param_float(2 + 5)*qst(1-2+1,j+1,indvarsst(11))*(d2_dif_rhoet_dydy_0_0_1m2p1jp1k)-&
                    q(1-2+1,j+1,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(1-2+1,j+1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1,j+1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1,j+1,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(1-2+1,j+1,indvarsst(11))*(d2_dif_rhou_dydy_0_0_1m2p1jp1k)+&
                    qst(1-2+1,j+1,indvarsst(10))*(d2_dif_rhou_dydx_0_0_1m2p1jp1k)))-&
                    q(1-2+1,j+1,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(1-2+1,j+1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1,j+1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1,j+1,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(1-2+1,j+1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_1m2p1jp1k)-&
                    2.0_wp/3.0_wp*(qst(1-2+1,j+1,indvarsst(10))*(d2_dif_rhov_dydx_0_0_1m2p1jp1k)+&
                    qst(1-2+1,j+1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_1m2p1jp1k))))

d1_dif_rhoet_dy_0_1m2p1jk = -&
          0.5_wp*d1_dif_rhoet_dy_0_1m2p1jm1k+&
          0.5_wp*d1_dif_rhoet_dy_0_1m2p1jp1k

d1_dif_rhoet_dy_0_1m2p1jk = d1_dif_rhoet_dy_0_1m2p1jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None d(rho et)/dt *******
!                                                           
!***********************************************************


rhs(1-2+1,j,indvars(4)) = rhs(1-2+1,j,indvars(4))  -  ( qst(1-2+1,j,indvarsst(10))*(d1_dif_rhoet_dx_0_1m2p1jk)+&
                    qst(1-2+1,j,indvarsst(11))*(d1_dif_rhoet_dy_0_1m2p1jk) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None d(rho nut)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([-(ReI*(1.0_wp+chi)*sigmaI*deltaxI*({nut}_1x))]_1x)+deltayI*([-(ReI*(1.0_wp+chi)*sigmaI*deltayI*({nut}_1y))]_1y)-ReI*Cb2*sigmaI*((deltaxI)**2*([rho*nut]_1x)*([nut]_1x)+(deltayI)**2*([rho*nut]_1y)*([nut]_1y))-Cb1*(1.0_wp-ft2)*SS*rho*nut+ReI*(Cw1*fw-Cb1/k**2.0_wp*ft2)*rho*(nut/eta)**2.0_wp
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_dif_rhonut_dxdx_0_0_1m2p1m1jk_1m2p1m1p0jk = q(1-2+1-1+0,j,indvars(5))

d2_dif_rhonut_dxdx_0_0_1m2p1m1jk_1m2p1m1p1jk = q(1-2+1-1+1,j,indvars(5))

d2_dif_rhonut_dxdx_0_0_1m2p1m1jk_1m2p1m1p2jk = q(1-2+1-1+2,j,indvars(5))

d2_dif_rhonut_dxdx_0_0_1m2p1m1jk = -&
          1.5_wp*d2_dif_rhonut_dxdx_0_0_1m2p1m1jk_1m2p1m1p0jk+&
          2.0_wp*d2_dif_rhonut_dxdx_0_0_1m2p1m1jk_1m2p1m1p1jk-&
          0.5_wp*d2_dif_rhonut_dxdx_0_0_1m2p1m1jk_1m2p1m1p2jk

d2_dif_rhonut_dxdx_0_0_1m2p1m1jk = d2_dif_rhonut_dxdx_0_0_1m2p1m1jk*param_float(1)

d2_dif_rhonut_dxdx_0_0_1m2p1p1jk_1m2p1p1m1jk = q(1-2+1+1-1,j,indvars(5))

d2_dif_rhonut_dxdx_0_0_1m2p1p1jk_1m2p1p1p1jk = q(1-2+1+1+1,j,indvars(5))

d2_dif_rhonut_dxdx_0_0_1m2p1p1jk = -&
          0.5_wp*d2_dif_rhonut_dxdx_0_0_1m2p1p1jk_1m2p1p1m1jk+&
          0.5_wp*d2_dif_rhonut_dxdx_0_0_1m2p1p1jk_1m2p1p1p1jk

d2_dif_rhonut_dxdx_0_0_1m2p1p1jk = d2_dif_rhonut_dxdx_0_0_1m2p1p1jk*param_float(1)

d1_dif_rhonut_dx_1_1m2p1m1jk = q(1-2+1-1,j,indvars(1))*q(1-2+1-1,j,indvars(5))

d1_dif_rhonut_dx_1_1m2p1p1jk = q(1-2+1+1,j,indvars(1))*q(1-2+1+1,j,indvars(5))

d1_dif_rhonut_dx_1_1m2p1jk = -&
          0.5_wp*d1_dif_rhonut_dx_1_1m2p1m1jk+&
          0.5_wp*d1_dif_rhonut_dx_1_1m2p1p1jk

d1_dif_rhonut_dx_1_1m2p1jk = d1_dif_rhonut_dx_1_1m2p1jk*param_float(1)

d1_dif_rhonut_dx_2_1m2p1m1jk = q(1-2+1-1,j,indvars(5))

d1_dif_rhonut_dx_2_1m2p1p1jk = q(1-2+1+1,j,indvars(5))

d1_dif_rhonut_dx_2_1m2p1jk = -&
          0.5_wp*d1_dif_rhonut_dx_2_1m2p1m1jk+&
          0.5_wp*d1_dif_rhonut_dx_2_1m2p1p1jk

d1_dif_rhonut_dx_2_1m2p1jk = d1_dif_rhonut_dx_2_1m2p1jk*param_float(1)

d1_dif_rhonut_dx_0_1m2p1m1jk = -(param_float(1 + 5)*(1.0_wp+&
                    (q(1-2+1-1,j,indvars(5))/1.0_wp))*param_float(18 + 5)*qst(1-2+1-1,j,indvarsst(10))*(d2_dif_rhonut_dxdx_0_0_1m2p1m1jk))

d1_dif_rhonut_dx_0_1m2p1p1jk = -(param_float(1 + 5)*(1.0_wp+&
                    (q(1-2+1+1,j,indvars(5))/1.0_wp))*param_float(18 + 5)*qst(1-2+1+1,j,indvarsst(10))*(d2_dif_rhonut_dxdx_0_0_1m2p1p1jk))

d1_dif_rhonut_dx_0_1m2p1jk = -&
          0.5_wp*d1_dif_rhonut_dx_0_1m2p1m1jk+&
          0.5_wp*d1_dif_rhonut_dx_0_1m2p1p1jk

d1_dif_rhonut_dx_0_1m2p1jk = d1_dif_rhonut_dx_0_1m2p1jk*param_float(1)

d2_dif_rhonut_dydy_0_0_1m2p1jm1k_1m2p1jm1m1k = q(1-2+1,j-1-1,indvars(5))

d2_dif_rhonut_dydy_0_0_1m2p1jm1k_1m2p1jm1p1k = q(1-2+1,j-1+1,indvars(5))

d2_dif_rhonut_dydy_0_0_1m2p1jm1k = -&
          0.5_wp*d2_dif_rhonut_dydy_0_0_1m2p1jm1k_1m2p1jm1m1k+&
          0.5_wp*d2_dif_rhonut_dydy_0_0_1m2p1jm1k_1m2p1jm1p1k

d2_dif_rhonut_dydy_0_0_1m2p1jm1k = d2_dif_rhonut_dydy_0_0_1m2p1jm1k*param_float(2)

d2_dif_rhonut_dydy_0_0_1m2p1jp1k_1m2p1jp1m1k = q(1-2+1,j+1-1,indvars(5))

d2_dif_rhonut_dydy_0_0_1m2p1jp1k_1m2p1jp1p1k = q(1-2+1,j+1+1,indvars(5))

d2_dif_rhonut_dydy_0_0_1m2p1jp1k = -&
          0.5_wp*d2_dif_rhonut_dydy_0_0_1m2p1jp1k_1m2p1jp1m1k+&
          0.5_wp*d2_dif_rhonut_dydy_0_0_1m2p1jp1k_1m2p1jp1p1k

d2_dif_rhonut_dydy_0_0_1m2p1jp1k = d2_dif_rhonut_dydy_0_0_1m2p1jp1k*param_float(2)

d1_dif_rhonut_dy_1_1m2p1jm1k = q(1-2+1,j-1,indvars(1))*q(1-2+1,j-1,indvars(5))

d1_dif_rhonut_dy_1_1m2p1jp1k = q(1-2+1,j+1,indvars(1))*q(1-2+1,j+1,indvars(5))

d1_dif_rhonut_dy_1_1m2p1jk = -&
          0.5_wp*d1_dif_rhonut_dy_1_1m2p1jm1k+&
          0.5_wp*d1_dif_rhonut_dy_1_1m2p1jp1k

d1_dif_rhonut_dy_1_1m2p1jk = d1_dif_rhonut_dy_1_1m2p1jk*param_float(2)

d1_dif_rhonut_dy_2_1m2p1jm1k = q(1-2+1,j-1,indvars(5))

d1_dif_rhonut_dy_2_1m2p1jp1k = q(1-2+1,j+1,indvars(5))

d1_dif_rhonut_dy_2_1m2p1jk = -&
          0.5_wp*d1_dif_rhonut_dy_2_1m2p1jm1k+&
          0.5_wp*d1_dif_rhonut_dy_2_1m2p1jp1k

d1_dif_rhonut_dy_2_1m2p1jk = d1_dif_rhonut_dy_2_1m2p1jk*param_float(2)

d1_dif_rhonut_dy_0_1m2p1jm1k = -(param_float(1 + 5)*(1.0_wp+&
                    (q(1-2+1,j-1,indvars(5))/1.0_wp))*param_float(18 + 5)*qst(1-2+1,j-1,indvarsst(11))*(d2_dif_rhonut_dydy_0_0_1m2p1jm1k))

d1_dif_rhonut_dy_0_1m2p1jp1k = -(param_float(1 + 5)*(1.0_wp+&
                    (q(1-2+1,j+1,indvars(5))/1.0_wp))*param_float(18 + 5)*qst(1-2+1,j+1,indvarsst(11))*(d2_dif_rhonut_dydy_0_0_1m2p1jp1k))

d1_dif_rhonut_dy_0_1m2p1jk = -&
          0.5_wp*d1_dif_rhonut_dy_0_1m2p1jm1k+&
          0.5_wp*d1_dif_rhonut_dy_0_1m2p1jp1k

d1_dif_rhonut_dy_0_1m2p1jk = d1_dif_rhonut_dy_0_1m2p1jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None d(rho nut)/dt ******
!                                                           
!***********************************************************


rhs(1-2+1,j,indvars(5)) = rhs(1-2+1,j,indvars(5))  -  ( qst(1-2+1,j,indvarsst(10))*(d1_dif_rhonut_dx_0_1m2p1jk)+&
                    qst(1-2+1,j,indvarsst(11))*(d1_dif_rhonut_dy_0_1m2p1jk)-&
                    param_float(1 + 5)*param_float(7 + 5)*param_float(18 + 5)*((qst(1-2+1,j,indvarsst(10)))**2*(d1_dif_rhonut_dx_1_1m2p1jk)*(d1_dif_rhonut_dx_2_1m2p1jk)+&
                    (qst(1-2+1,j,indvarsst(11)))**2*(d1_dif_rhonut_dy_1_1m2p1jk)*(d1_dif_rhonut_dy_2_1m2p1jk))-&
                    param_float(6 + 5)*(1.0_wp-&
                    (param_float(16 + 5)*exp(-&
                    param_float(17 + 5)*(q(1-2+1,j,indvars(5))/1.0_wp)**2.0_wp)))*qst(1-2+1,j,indvarsst(13))*q(1-2+1,j,indvars(1))*q(1-2+1,j,indvars(5))+&
                    param_float(1 + 5)*(param_float(10 + 5)*qst(1-2+1,j,indvarsst(16))-&
                    param_float(6 + 5)/param_float(9 + 5)**2.0_wp*(param_float(16 + 5)*exp(-&
                    param_float(17 + 5)*(q(1-2+1,j,indvars(5))/1.0_wp)**2.0_wp)))*q(1-2+1,j,indvars(1))*(q(1-2+1,j,indvars(5))/qst(1-2+1,j,indvarsst(2)))**2.0_wp ) 

     enddo
