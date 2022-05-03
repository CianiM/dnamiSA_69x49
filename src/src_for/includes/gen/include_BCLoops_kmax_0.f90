

!***********************************************************
!                                                           
! Start building layers for BC : None None kmax ************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: None None 0 ************************************
!                                                           
!***********************************************************


     do j=idloop(3),idloop(4) 
 
      do i=idloop(1),idloop(2) 


!***********************************************************
!                                                           
! building source terms in RHS for layer None None 0 d(rho)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([rho*u]_1x)+deltayI*([rho*v]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_conv_rho_dx_0_im1jnzp2 = q(i-1,j,indvars(1))*q(i-1,j,indvars(2))

d1_conv_rho_dx_0_ip1jnzp2 = q(i+1,j,indvars(1))*q(i+1,j,indvars(2))

d1_conv_rho_dx_0_ijnzp2 = -&
          0.5_wp*d1_conv_rho_dx_0_im1jnzp2+&
          0.5_wp*d1_conv_rho_dx_0_ip1jnzp2

d1_conv_rho_dx_0_ijnzp2 = d1_conv_rho_dx_0_ijnzp2*param_float(1)

d1_conv_rho_dy_0_ijm1nzp2 = q(i,j-1,indvars(1))*q(i,j-1,indvars(3))

d1_conv_rho_dy_0_ijp1nzp2 = q(i,j+1,indvars(1))*q(i,j+1,indvars(3))

d1_conv_rho_dy_0_ijnzp2 = -&
          0.5_wp*d1_conv_rho_dy_0_ijm1nzp2+&
          0.5_wp*d1_conv_rho_dy_0_ijp1nzp2

d1_conv_rho_dy_0_ijnzp2 = d1_conv_rho_dy_0_ijnzp2*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None None 0 d(rho)/dt **********
!                                                           
!***********************************************************


rhs(i,j,indvars(1)) =   -  ( qst(i,j,indvarsst(10))*(d1_conv_rho_dx_0_ijnzp2)+&
                    qst(i,j,indvarsst(11))*(d1_conv_rho_dy_0_ijnzp2) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer None None 0 d(rho u)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([rho*u*u+p]_1x)+deltayI*([rho*v*u]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_conv_rhou_dx_0_im1jnzp2 = q(i-1,j,indvars(1))*q(i-1,j,indvars(2))*q(i-1,j,indvars(2))+(param_float(3 + 5))*q(i-1,j,indvars(1))*((q(i-1,j,indvars(4))-&
                    0.5_wp*(q(i-1,j,indvars(2))*q(i-1,j,indvars(2))+&
                    q(i-1,j,indvars(3))*q(i-1,j,indvars(3)))))

d1_conv_rhou_dx_0_ip1jnzp2 = q(i+1,j,indvars(1))*q(i+1,j,indvars(2))*q(i+1,j,indvars(2))+(param_float(3 + 5))*q(i+1,j,indvars(1))*((q(i+1,j,indvars(4))-&
                    0.5_wp*(q(i+1,j,indvars(2))*q(i+1,j,indvars(2))+&
                    q(i+1,j,indvars(3))*q(i+1,j,indvars(3)))))

d1_conv_rhou_dx_0_ijnzp2 = -&
          0.5_wp*d1_conv_rhou_dx_0_im1jnzp2+&
          0.5_wp*d1_conv_rhou_dx_0_ip1jnzp2

d1_conv_rhou_dx_0_ijnzp2 = d1_conv_rhou_dx_0_ijnzp2*param_float(1)

d1_conv_rhou_dy_0_ijm1nzp2 = q(i,j-1,indvars(1))*q(i,j-1,indvars(3))*q(i,j-1,indvars(2))

d1_conv_rhou_dy_0_ijp1nzp2 = q(i,j+1,indvars(1))*q(i,j+1,indvars(3))*q(i,j+1,indvars(2))

d1_conv_rhou_dy_0_ijnzp2 = -&
          0.5_wp*d1_conv_rhou_dy_0_ijm1nzp2+&
          0.5_wp*d1_conv_rhou_dy_0_ijp1nzp2

d1_conv_rhou_dy_0_ijnzp2 = d1_conv_rhou_dy_0_ijnzp2*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None None 0 d(rho u)/dt ********
!                                                           
!***********************************************************


rhs(i,j,indvars(2)) =   -  ( qst(i,j,indvarsst(10))*(d1_conv_rhou_dx_0_ijnzp2)+&
                    qst(i,j,indvarsst(11))*(d1_conv_rhou_dy_0_ijnzp2) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer None None 0 d(rho v)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([rho*u*v]_1x)+deltayI*([rho*v*v+p]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_conv_rhov_dx_0_im1jnzp2 = q(i-1,j,indvars(1))*q(i-1,j,indvars(2))*q(i-1,j,indvars(3))

d1_conv_rhov_dx_0_ip1jnzp2 = q(i+1,j,indvars(1))*q(i+1,j,indvars(2))*q(i+1,j,indvars(3))

d1_conv_rhov_dx_0_ijnzp2 = -&
          0.5_wp*d1_conv_rhov_dx_0_im1jnzp2+&
          0.5_wp*d1_conv_rhov_dx_0_ip1jnzp2

d1_conv_rhov_dx_0_ijnzp2 = d1_conv_rhov_dx_0_ijnzp2*param_float(1)

d1_conv_rhov_dy_0_ijm1nzp2 = q(i,j-1,indvars(1))*q(i,j-1,indvars(3))*q(i,j-1,indvars(3))+(param_float(3 + 5))*q(i,j-1,indvars(1))*((q(i,j-1,indvars(4))-&
                    0.5_wp*(q(i,j-1,indvars(2))*q(i,j-1,indvars(2))+&
                    q(i,j-1,indvars(3))*q(i,j-1,indvars(3)))))

d1_conv_rhov_dy_0_ijp1nzp2 = q(i,j+1,indvars(1))*q(i,j+1,indvars(3))*q(i,j+1,indvars(3))+(param_float(3 + 5))*q(i,j+1,indvars(1))*((q(i,j+1,indvars(4))-&
                    0.5_wp*(q(i,j+1,indvars(2))*q(i,j+1,indvars(2))+&
                    q(i,j+1,indvars(3))*q(i,j+1,indvars(3)))))

d1_conv_rhov_dy_0_ijnzp2 = -&
          0.5_wp*d1_conv_rhov_dy_0_ijm1nzp2+&
          0.5_wp*d1_conv_rhov_dy_0_ijp1nzp2

d1_conv_rhov_dy_0_ijnzp2 = d1_conv_rhov_dy_0_ijnzp2*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None None 0 d(rho v)/dt ********
!                                                           
!***********************************************************


rhs(i,j,indvars(3)) =   -  ( qst(i,j,indvarsst(10))*(d1_conv_rhov_dx_0_ijnzp2)+&
                    qst(i,j,indvarsst(11))*(d1_conv_rhov_dy_0_ijnzp2) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer None None 0 d(rho et)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([(rho*et+p)*u]_1x)+deltayI*([(rho*et+p)*v]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_conv_rhoet_dx_0_im1jnzp2 = (q(i-1,j,indvars(1))*q(i-1,j,indvars(4))+&
                    (param_float(3 + 5))*q(i-1,j,indvars(1))*((q(i-1,j,indvars(4))-&
                    0.5_wp*(q(i-1,j,indvars(2))*q(i-1,j,indvars(2))+&
                    q(i-1,j,indvars(3))*q(i-1,j,indvars(3))))))*q(i-1,j,indvars(2))

d1_conv_rhoet_dx_0_ip1jnzp2 = (q(i+1,j,indvars(1))*q(i+1,j,indvars(4))+&
                    (param_float(3 + 5))*q(i+1,j,indvars(1))*((q(i+1,j,indvars(4))-&
                    0.5_wp*(q(i+1,j,indvars(2))*q(i+1,j,indvars(2))+&
                    q(i+1,j,indvars(3))*q(i+1,j,indvars(3))))))*q(i+1,j,indvars(2))

d1_conv_rhoet_dx_0_ijnzp2 = -&
          0.5_wp*d1_conv_rhoet_dx_0_im1jnzp2+&
          0.5_wp*d1_conv_rhoet_dx_0_ip1jnzp2

d1_conv_rhoet_dx_0_ijnzp2 = d1_conv_rhoet_dx_0_ijnzp2*param_float(1)

d1_conv_rhoet_dy_0_ijm1nzp2 = (q(i,j-1,indvars(1))*q(i,j-1,indvars(4))+&
                    (param_float(3 + 5))*q(i,j-1,indvars(1))*((q(i,j-1,indvars(4))-&
                    0.5_wp*(q(i,j-1,indvars(2))*q(i,j-1,indvars(2))+&
                    q(i,j-1,indvars(3))*q(i,j-1,indvars(3))))))*q(i,j-1,indvars(3))

d1_conv_rhoet_dy_0_ijp1nzp2 = (q(i,j+1,indvars(1))*q(i,j+1,indvars(4))+&
                    (param_float(3 + 5))*q(i,j+1,indvars(1))*((q(i,j+1,indvars(4))-&
                    0.5_wp*(q(i,j+1,indvars(2))*q(i,j+1,indvars(2))+&
                    q(i,j+1,indvars(3))*q(i,j+1,indvars(3))))))*q(i,j+1,indvars(3))

d1_conv_rhoet_dy_0_ijnzp2 = -&
          0.5_wp*d1_conv_rhoet_dy_0_ijm1nzp2+&
          0.5_wp*d1_conv_rhoet_dy_0_ijp1nzp2

d1_conv_rhoet_dy_0_ijnzp2 = d1_conv_rhoet_dy_0_ijnzp2*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None None 0 d(rho et)/dt *******
!                                                           
!***********************************************************


rhs(i,j,indvars(4)) =   -  ( qst(i,j,indvarsst(10))*(d1_conv_rhoet_dx_0_ijnzp2)+&
                    qst(i,j,indvarsst(11))*(d1_conv_rhoet_dy_0_ijnzp2) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer None None 0 d(rho nut)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([(rho*nut)*u]_1x)+deltayI*([(rho*nut)*v]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_conv_rhonut_dx_0_im1jnzp2 = (q(i-1,j,indvars(1))*q(i-1,j,indvars(5)))*q(i-1,j,indvars(2))

d1_conv_rhonut_dx_0_ip1jnzp2 = (q(i+1,j,indvars(1))*q(i+1,j,indvars(5)))*q(i+1,j,indvars(2))

d1_conv_rhonut_dx_0_ijnzp2 = -&
          0.5_wp*d1_conv_rhonut_dx_0_im1jnzp2+&
          0.5_wp*d1_conv_rhonut_dx_0_ip1jnzp2

d1_conv_rhonut_dx_0_ijnzp2 = d1_conv_rhonut_dx_0_ijnzp2*param_float(1)

d1_conv_rhonut_dy_0_ijm1nzp2 = (q(i,j-1,indvars(1))*q(i,j-1,indvars(5)))*q(i,j-1,indvars(3))

d1_conv_rhonut_dy_0_ijp1nzp2 = (q(i,j+1,indvars(1))*q(i,j+1,indvars(5)))*q(i,j+1,indvars(3))

d1_conv_rhonut_dy_0_ijnzp2 = -&
          0.5_wp*d1_conv_rhonut_dy_0_ijm1nzp2+&
          0.5_wp*d1_conv_rhonut_dy_0_ijp1nzp2

d1_conv_rhonut_dy_0_ijnzp2 = d1_conv_rhonut_dy_0_ijnzp2*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None None 0 d(rho nut)/dt ******
!                                                           
!***********************************************************


rhs(i,j,indvars(5)) =   -  ( qst(i,j,indvarsst(10))*(d1_conv_rhonut_dx_0_ijnzp2)+&
                    qst(i,j,indvarsst(11))*(d1_conv_rhonut_dy_0_ijnzp2) ) 

     enddo
   enddo


!***********************************************************
!                                                           
! Start building layers for BC : None None kmax ************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: None None 0 ************************************
!                                                           
!***********************************************************


     do j=idloop(3),idloop(4) 
 
      do i=idloop(1),idloop(2) 


!***********************************************************
!                                                           
! building source terms in RHS for layer None None 0 d(rho u)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([-visc_t*(2.0_wp*deltaxI*({u}_1x)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y)))]_1x)+deltayI*([-visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_dif_rhou_dxdx_0_0_im1jnzp2_im1m1jnzp2 = q(i-1-1,j,indvars(2))

d2_dif_rhou_dxdx_0_0_im1jnzp2_im1p1jnzp2 = q(i-1+1,j,indvars(2))

d2_dif_rhou_dxdx_0_0_im1jnzp2 = -&
          0.5_wp*d2_dif_rhou_dxdx_0_0_im1jnzp2_im1m1jnzp2+&
          0.5_wp*d2_dif_rhou_dxdx_0_0_im1jnzp2_im1p1jnzp2

d2_dif_rhou_dxdx_0_0_im1jnzp2 = d2_dif_rhou_dxdx_0_0_im1jnzp2*param_float(1)

d2_dif_rhou_dxdx_0_0_ip1jnzp2_ip1m1jnzp2 = q(i+1-1,j,indvars(2))

d2_dif_rhou_dxdx_0_0_ip1jnzp2_ip1p1jnzp2 = q(i+1+1,j,indvars(2))

d2_dif_rhou_dxdx_0_0_ip1jnzp2 = -&
          0.5_wp*d2_dif_rhou_dxdx_0_0_ip1jnzp2_ip1m1jnzp2+&
          0.5_wp*d2_dif_rhou_dxdx_0_0_ip1jnzp2_ip1p1jnzp2

d2_dif_rhou_dxdx_0_0_ip1jnzp2 = d2_dif_rhou_dxdx_0_0_ip1jnzp2*param_float(1)

d2_dif_rhou_dxdy_0_0_im1jnzp2_im1jm1nzp2 = q(i-1,j-1,indvars(3))

d2_dif_rhou_dxdy_0_0_im1jnzp2_im1jp1nzp2 = q(i-1,j+1,indvars(3))

d2_dif_rhou_dxdy_0_0_im1jnzp2 = -&
          0.5_wp*d2_dif_rhou_dxdy_0_0_im1jnzp2_im1jm1nzp2+&
          0.5_wp*d2_dif_rhou_dxdy_0_0_im1jnzp2_im1jp1nzp2

d2_dif_rhou_dxdy_0_0_im1jnzp2 = d2_dif_rhou_dxdy_0_0_im1jnzp2*param_float(2)

d2_dif_rhou_dxdy_0_0_ip1jnzp2_ip1jm1nzp2 = q(i+1,j-1,indvars(3))

d2_dif_rhou_dxdy_0_0_ip1jnzp2_ip1jp1nzp2 = q(i+1,j+1,indvars(3))

d2_dif_rhou_dxdy_0_0_ip1jnzp2 = -&
          0.5_wp*d2_dif_rhou_dxdy_0_0_ip1jnzp2_ip1jm1nzp2+&
          0.5_wp*d2_dif_rhou_dxdy_0_0_ip1jnzp2_ip1jp1nzp2

d2_dif_rhou_dxdy_0_0_ip1jnzp2 = d2_dif_rhou_dxdy_0_0_ip1jnzp2*param_float(2)

d1_dif_rhou_dx_0_im1jnzp2 = -(1.0_wp)*(1.0_wp+&
                    ((q(i-1,j,indvars(5))/1.0_wp)**3.0_wp/((q(i-1,j,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i-1,j,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(i-1,j,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_im1jnzp2)-&
                    2.0_wp/3.0_wp*(qst(i-1,j,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_im1jnzp2)+&
                    qst(i-1,j,indvarsst(11))*(d2_dif_rhou_dxdy_0_0_im1jnzp2)))

d1_dif_rhou_dx_0_ip1jnzp2 = -(1.0_wp)*(1.0_wp+&
                    ((q(i+1,j,indvars(5))/1.0_wp)**3.0_wp/((q(i+1,j,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i+1,j,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(i+1,j,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_ip1jnzp2)-&
                    2.0_wp/3.0_wp*(qst(i+1,j,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_ip1jnzp2)+&
                    qst(i+1,j,indvarsst(11))*(d2_dif_rhou_dxdy_0_0_ip1jnzp2)))

d1_dif_rhou_dx_0_ijnzp2 = -&
          0.5_wp*d1_dif_rhou_dx_0_im1jnzp2+&
          0.5_wp*d1_dif_rhou_dx_0_ip1jnzp2

d1_dif_rhou_dx_0_ijnzp2 = d1_dif_rhou_dx_0_ijnzp2*param_float(1)

d2_dif_rhou_dydx_0_0_ijm1nzp2_im1jm1nzp2 = q(i-1,j-1,indvars(3))

d2_dif_rhou_dydx_0_0_ijm1nzp2_ip1jm1nzp2 = q(i+1,j-1,indvars(3))

d2_dif_rhou_dydx_0_0_ijm1nzp2 = -&
          0.5_wp*d2_dif_rhou_dydx_0_0_ijm1nzp2_im1jm1nzp2+&
          0.5_wp*d2_dif_rhou_dydx_0_0_ijm1nzp2_ip1jm1nzp2

d2_dif_rhou_dydx_0_0_ijm1nzp2 = d2_dif_rhou_dydx_0_0_ijm1nzp2*param_float(1)

d2_dif_rhou_dydx_0_0_ijp1nzp2_im1jp1nzp2 = q(i-1,j+1,indvars(3))

d2_dif_rhou_dydx_0_0_ijp1nzp2_ip1jp1nzp2 = q(i+1,j+1,indvars(3))

d2_dif_rhou_dydx_0_0_ijp1nzp2 = -&
          0.5_wp*d2_dif_rhou_dydx_0_0_ijp1nzp2_im1jp1nzp2+&
          0.5_wp*d2_dif_rhou_dydx_0_0_ijp1nzp2_ip1jp1nzp2

d2_dif_rhou_dydx_0_0_ijp1nzp2 = d2_dif_rhou_dydx_0_0_ijp1nzp2*param_float(1)

d2_dif_rhou_dydy_0_0_ijm1nzp2_ijm1m1nzp2 = q(i,j-1-1,indvars(2))

d2_dif_rhou_dydy_0_0_ijm1nzp2_ijm1p1nzp2 = q(i,j-1+1,indvars(2))

d2_dif_rhou_dydy_0_0_ijm1nzp2 = -&
          0.5_wp*d2_dif_rhou_dydy_0_0_ijm1nzp2_ijm1m1nzp2+&
          0.5_wp*d2_dif_rhou_dydy_0_0_ijm1nzp2_ijm1p1nzp2

d2_dif_rhou_dydy_0_0_ijm1nzp2 = d2_dif_rhou_dydy_0_0_ijm1nzp2*param_float(2)

d2_dif_rhou_dydy_0_0_ijp1nzp2_ijp1m1nzp2 = q(i,j+1-1,indvars(2))

d2_dif_rhou_dydy_0_0_ijp1nzp2_ijp1p1nzp2 = q(i,j+1+1,indvars(2))

d2_dif_rhou_dydy_0_0_ijp1nzp2 = -&
          0.5_wp*d2_dif_rhou_dydy_0_0_ijp1nzp2_ijp1m1nzp2+&
          0.5_wp*d2_dif_rhou_dydy_0_0_ijp1nzp2_ijp1p1nzp2

d2_dif_rhou_dydy_0_0_ijp1nzp2 = d2_dif_rhou_dydy_0_0_ijp1nzp2*param_float(2)

d1_dif_rhou_dy_0_ijm1nzp2 = -(1.0_wp)*(1.0_wp+&
                    ((q(i,j-1,indvars(5))/1.0_wp)**3.0_wp/((q(i,j-1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i,j-1,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(i,j-1,indvarsst(11))*(d2_dif_rhou_dydy_0_0_ijm1nzp2)+&
                    qst(i,j-1,indvarsst(10))*(d2_dif_rhou_dydx_0_0_ijm1nzp2))

d1_dif_rhou_dy_0_ijp1nzp2 = -(1.0_wp)*(1.0_wp+&
                    ((q(i,j+1,indvars(5))/1.0_wp)**3.0_wp/((q(i,j+1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i,j+1,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(i,j+1,indvarsst(11))*(d2_dif_rhou_dydy_0_0_ijp1nzp2)+&
                    qst(i,j+1,indvarsst(10))*(d2_dif_rhou_dydx_0_0_ijp1nzp2))

d1_dif_rhou_dy_0_ijnzp2 = -&
          0.5_wp*d1_dif_rhou_dy_0_ijm1nzp2+&
          0.5_wp*d1_dif_rhou_dy_0_ijp1nzp2

d1_dif_rhou_dy_0_ijnzp2 = d1_dif_rhou_dy_0_ijnzp2*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None None 0 d(rho u)/dt ********
!                                                           
!***********************************************************


rhs(i,j,indvars(2)) = rhs(i,j,indvars(2))  -  ( qst(i,j,indvarsst(10))*(d1_dif_rhou_dx_0_ijnzp2)+&
                    qst(i,j,indvarsst(11))*(d1_dif_rhou_dy_0_ijnzp2) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer None None 0 d(rho v)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([-visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x))]_1x)+deltayI*([-visc_t*(2.0_wp*deltayI*({v}_1y)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y)))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_dif_rhov_dxdx_0_0_im1jnzp2_im1m1jnzp2 = q(i-1-1,j,indvars(3))

d2_dif_rhov_dxdx_0_0_im1jnzp2_im1p1jnzp2 = q(i-1+1,j,indvars(3))

d2_dif_rhov_dxdx_0_0_im1jnzp2 = -&
          0.5_wp*d2_dif_rhov_dxdx_0_0_im1jnzp2_im1m1jnzp2+&
          0.5_wp*d2_dif_rhov_dxdx_0_0_im1jnzp2_im1p1jnzp2

d2_dif_rhov_dxdx_0_0_im1jnzp2 = d2_dif_rhov_dxdx_0_0_im1jnzp2*param_float(1)

d2_dif_rhov_dxdx_0_0_ip1jnzp2_ip1m1jnzp2 = q(i+1-1,j,indvars(3))

d2_dif_rhov_dxdx_0_0_ip1jnzp2_ip1p1jnzp2 = q(i+1+1,j,indvars(3))

d2_dif_rhov_dxdx_0_0_ip1jnzp2 = -&
          0.5_wp*d2_dif_rhov_dxdx_0_0_ip1jnzp2_ip1m1jnzp2+&
          0.5_wp*d2_dif_rhov_dxdx_0_0_ip1jnzp2_ip1p1jnzp2

d2_dif_rhov_dxdx_0_0_ip1jnzp2 = d2_dif_rhov_dxdx_0_0_ip1jnzp2*param_float(1)

d2_dif_rhov_dxdy_0_0_im1jnzp2_im1jm1nzp2 = q(i-1,j-1,indvars(2))

d2_dif_rhov_dxdy_0_0_im1jnzp2_im1jp1nzp2 = q(i-1,j+1,indvars(2))

d2_dif_rhov_dxdy_0_0_im1jnzp2 = -&
          0.5_wp*d2_dif_rhov_dxdy_0_0_im1jnzp2_im1jm1nzp2+&
          0.5_wp*d2_dif_rhov_dxdy_0_0_im1jnzp2_im1jp1nzp2

d2_dif_rhov_dxdy_0_0_im1jnzp2 = d2_dif_rhov_dxdy_0_0_im1jnzp2*param_float(2)

d2_dif_rhov_dxdy_0_0_ip1jnzp2_ip1jm1nzp2 = q(i+1,j-1,indvars(2))

d2_dif_rhov_dxdy_0_0_ip1jnzp2_ip1jp1nzp2 = q(i+1,j+1,indvars(2))

d2_dif_rhov_dxdy_0_0_ip1jnzp2 = -&
          0.5_wp*d2_dif_rhov_dxdy_0_0_ip1jnzp2_ip1jm1nzp2+&
          0.5_wp*d2_dif_rhov_dxdy_0_0_ip1jnzp2_ip1jp1nzp2

d2_dif_rhov_dxdy_0_0_ip1jnzp2 = d2_dif_rhov_dxdy_0_0_ip1jnzp2*param_float(2)

d1_dif_rhov_dx_0_im1jnzp2 = -(1.0_wp)*(1.0_wp+&
                    ((q(i-1,j,indvars(5))/1.0_wp)**3.0_wp/((q(i-1,j,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i-1,j,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(i-1,j,indvarsst(11))*(d2_dif_rhov_dxdy_0_0_im1jnzp2)+&
                    qst(i-1,j,indvarsst(10))*(d2_dif_rhov_dxdx_0_0_im1jnzp2))

d1_dif_rhov_dx_0_ip1jnzp2 = -(1.0_wp)*(1.0_wp+&
                    ((q(i+1,j,indvars(5))/1.0_wp)**3.0_wp/((q(i+1,j,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i+1,j,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(i+1,j,indvarsst(11))*(d2_dif_rhov_dxdy_0_0_ip1jnzp2)+&
                    qst(i+1,j,indvarsst(10))*(d2_dif_rhov_dxdx_0_0_ip1jnzp2))

d1_dif_rhov_dx_0_ijnzp2 = -&
          0.5_wp*d1_dif_rhov_dx_0_im1jnzp2+&
          0.5_wp*d1_dif_rhov_dx_0_ip1jnzp2

d1_dif_rhov_dx_0_ijnzp2 = d1_dif_rhov_dx_0_ijnzp2*param_float(1)

d2_dif_rhov_dydx_0_0_ijm1nzp2_im1jm1nzp2 = q(i-1,j-1,indvars(2))

d2_dif_rhov_dydx_0_0_ijm1nzp2_ip1jm1nzp2 = q(i+1,j-1,indvars(2))

d2_dif_rhov_dydx_0_0_ijm1nzp2 = -&
          0.5_wp*d2_dif_rhov_dydx_0_0_ijm1nzp2_im1jm1nzp2+&
          0.5_wp*d2_dif_rhov_dydx_0_0_ijm1nzp2_ip1jm1nzp2

d2_dif_rhov_dydx_0_0_ijm1nzp2 = d2_dif_rhov_dydx_0_0_ijm1nzp2*param_float(1)

d2_dif_rhov_dydx_0_0_ijp1nzp2_im1jp1nzp2 = q(i-1,j+1,indvars(2))

d2_dif_rhov_dydx_0_0_ijp1nzp2_ip1jp1nzp2 = q(i+1,j+1,indvars(2))

d2_dif_rhov_dydx_0_0_ijp1nzp2 = -&
          0.5_wp*d2_dif_rhov_dydx_0_0_ijp1nzp2_im1jp1nzp2+&
          0.5_wp*d2_dif_rhov_dydx_0_0_ijp1nzp2_ip1jp1nzp2

d2_dif_rhov_dydx_0_0_ijp1nzp2 = d2_dif_rhov_dydx_0_0_ijp1nzp2*param_float(1)

d2_dif_rhov_dydy_0_0_ijm1nzp2_ijm1m1nzp2 = q(i,j-1-1,indvars(3))

d2_dif_rhov_dydy_0_0_ijm1nzp2_ijm1p1nzp2 = q(i,j-1+1,indvars(3))

d2_dif_rhov_dydy_0_0_ijm1nzp2 = -&
          0.5_wp*d2_dif_rhov_dydy_0_0_ijm1nzp2_ijm1m1nzp2+&
          0.5_wp*d2_dif_rhov_dydy_0_0_ijm1nzp2_ijm1p1nzp2

d2_dif_rhov_dydy_0_0_ijm1nzp2 = d2_dif_rhov_dydy_0_0_ijm1nzp2*param_float(2)

d2_dif_rhov_dydy_0_0_ijp1nzp2_ijp1m1nzp2 = q(i,j+1-1,indvars(3))

d2_dif_rhov_dydy_0_0_ijp1nzp2_ijp1p1nzp2 = q(i,j+1+1,indvars(3))

d2_dif_rhov_dydy_0_0_ijp1nzp2 = -&
          0.5_wp*d2_dif_rhov_dydy_0_0_ijp1nzp2_ijp1m1nzp2+&
          0.5_wp*d2_dif_rhov_dydy_0_0_ijp1nzp2_ijp1p1nzp2

d2_dif_rhov_dydy_0_0_ijp1nzp2 = d2_dif_rhov_dydy_0_0_ijp1nzp2*param_float(2)

d1_dif_rhov_dy_0_ijm1nzp2 = -(1.0_wp)*(1.0_wp+&
                    ((q(i,j-1,indvars(5))/1.0_wp)**3.0_wp/((q(i,j-1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i,j-1,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(i,j-1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_ijm1nzp2)-&
                    2.0_wp/3.0_wp*(qst(i,j-1,indvarsst(10))*(d2_dif_rhov_dydx_0_0_ijm1nzp2)+&
                    qst(i,j-1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_ijm1nzp2)))

d1_dif_rhov_dy_0_ijp1nzp2 = -(1.0_wp)*(1.0_wp+&
                    ((q(i,j+1,indvars(5))/1.0_wp)**3.0_wp/((q(i,j+1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i,j+1,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(i,j+1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_ijp1nzp2)-&
                    2.0_wp/3.0_wp*(qst(i,j+1,indvarsst(10))*(d2_dif_rhov_dydx_0_0_ijp1nzp2)+&
                    qst(i,j+1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_ijp1nzp2)))

d1_dif_rhov_dy_0_ijnzp2 = -&
          0.5_wp*d1_dif_rhov_dy_0_ijm1nzp2+&
          0.5_wp*d1_dif_rhov_dy_0_ijp1nzp2

d1_dif_rhov_dy_0_ijnzp2 = d1_dif_rhov_dy_0_ijnzp2*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None None 0 d(rho v)/dt ********
!                                                           
!***********************************************************


rhs(i,j,indvars(3)) = rhs(i,j,indvars(3))  -  ( qst(i,j,indvarsst(10))*(d1_dif_rhov_dx_0_ijnzp2)+&
                    qst(i,j,indvarsst(11))*(d1_dif_rhov_dy_0_ijnzp2) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer None None 0 d(rho et)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([-kappa*deltaxI*({T}_1x)-u*(visc_t*(2.0_wp*deltaxI*({u}_1x)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y))))-v*(visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x)))]_1x)+deltayI*([-kappa*deltayI*({T}_1y)-u*(visc_t*(deltayI*({u}_1y)+deltaxI*({v}_1x)))-v*(visc_t*(2.0_wp*deltayI*({v}_1y)-2.0_wp/3.0_wp*(deltaxI*({u}_1x)+deltayI*({v}_1y))))]_1y)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_dif_rhoet_dxdx_0_0_im1jnzp2_im1m1jnzp2 = ((q(i-1-1,j,indvars(4))-&
                    0.5_wp*(q(i-1-1,j,indvars(2))*q(i-1-1,j,indvars(2))+&
                    q(i-1-1,j,indvars(3))*q(i-1-1,j,indvars(3)))))/(param_float(4 + 5))

d2_dif_rhoet_dxdx_0_0_im1jnzp2_im1p1jnzp2 = ((q(i-1+1,j,indvars(4))-&
                    0.5_wp*(q(i-1+1,j,indvars(2))*q(i-1+1,j,indvars(2))+&
                    q(i-1+1,j,indvars(3))*q(i-1+1,j,indvars(3)))))/(param_float(4 + 5))

d2_dif_rhoet_dxdx_0_0_im1jnzp2 = -&
          0.5_wp*d2_dif_rhoet_dxdx_0_0_im1jnzp2_im1m1jnzp2+&
          0.5_wp*d2_dif_rhoet_dxdx_0_0_im1jnzp2_im1p1jnzp2

d2_dif_rhoet_dxdx_0_0_im1jnzp2 = d2_dif_rhoet_dxdx_0_0_im1jnzp2*param_float(1)

d2_dif_rhoet_dxdx_0_0_ip1jnzp2_ip1m1jnzp2 = ((q(i+1-1,j,indvars(4))-&
                    0.5_wp*(q(i+1-1,j,indvars(2))*q(i+1-1,j,indvars(2))+&
                    q(i+1-1,j,indvars(3))*q(i+1-1,j,indvars(3)))))/(param_float(4 + 5))

d2_dif_rhoet_dxdx_0_0_ip1jnzp2_ip1p1jnzp2 = ((q(i+1+1,j,indvars(4))-&
                    0.5_wp*(q(i+1+1,j,indvars(2))*q(i+1+1,j,indvars(2))+&
                    q(i+1+1,j,indvars(3))*q(i+1+1,j,indvars(3)))))/(param_float(4 + 5))

d2_dif_rhoet_dxdx_0_0_ip1jnzp2 = -&
          0.5_wp*d2_dif_rhoet_dxdx_0_0_ip1jnzp2_ip1m1jnzp2+&
          0.5_wp*d2_dif_rhoet_dxdx_0_0_ip1jnzp2_ip1p1jnzp2

d2_dif_rhoet_dxdx_0_0_ip1jnzp2 = d2_dif_rhoet_dxdx_0_0_ip1jnzp2*param_float(1)

d1_dif_rhoet_dx_0_im1jnzp2 = -param_float(2 + 5)*qst(i-1,j,indvarsst(10))*(d2_dif_rhoet_dxdx_0_0_im1jnzp2)-&
                    q(i-1,j,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(i-1,j,indvars(5))/1.0_wp)**3.0_wp/((q(i-1,j,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i-1,j,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(i-1,j,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_im1jnzp2)-&
                    2.0_wp/3.0_wp*(qst(i-1,j,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_im1jnzp2)+&
                    qst(i-1,j,indvarsst(11))*(d2_dif_rhou_dxdy_0_0_im1jnzp2))))-&
                    q(i-1,j,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(i-1,j,indvars(5))/1.0_wp)**3.0_wp/((q(i-1,j,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i-1,j,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(i-1,j,indvarsst(11))*(d2_dif_rhov_dxdy_0_0_im1jnzp2)+&
                    qst(i-1,j,indvarsst(10))*(d2_dif_rhov_dxdx_0_0_im1jnzp2)))

d1_dif_rhoet_dx_0_ip1jnzp2 = -param_float(2 + 5)*qst(i+1,j,indvarsst(10))*(d2_dif_rhoet_dxdx_0_0_ip1jnzp2)-&
                    q(i+1,j,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(i+1,j,indvars(5))/1.0_wp)**3.0_wp/((q(i+1,j,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i+1,j,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(i+1,j,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_ip1jnzp2)-&
                    2.0_wp/3.0_wp*(qst(i+1,j,indvarsst(10))*(d2_dif_rhou_dxdx_0_0_ip1jnzp2)+&
                    qst(i+1,j,indvarsst(11))*(d2_dif_rhou_dxdy_0_0_ip1jnzp2))))-&
                    q(i+1,j,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(i+1,j,indvars(5))/1.0_wp)**3.0_wp/((q(i+1,j,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i+1,j,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(i+1,j,indvarsst(11))*(d2_dif_rhov_dxdy_0_0_ip1jnzp2)+&
                    qst(i+1,j,indvarsst(10))*(d2_dif_rhov_dxdx_0_0_ip1jnzp2)))

d1_dif_rhoet_dx_0_ijnzp2 = -&
          0.5_wp*d1_dif_rhoet_dx_0_im1jnzp2+&
          0.5_wp*d1_dif_rhoet_dx_0_ip1jnzp2

d1_dif_rhoet_dx_0_ijnzp2 = d1_dif_rhoet_dx_0_ijnzp2*param_float(1)

d2_dif_rhoet_dydy_0_0_ijm1nzp2_ijm1m1nzp2 = ((q(i,j-1-1,indvars(4))-&
                    0.5_wp*(q(i,j-1-1,indvars(2))*q(i,j-1-1,indvars(2))+&
                    q(i,j-1-1,indvars(3))*q(i,j-1-1,indvars(3)))))/(param_float(4 + 5))

d2_dif_rhoet_dydy_0_0_ijm1nzp2_ijm1p1nzp2 = ((q(i,j-1+1,indvars(4))-&
                    0.5_wp*(q(i,j-1+1,indvars(2))*q(i,j-1+1,indvars(2))+&
                    q(i,j-1+1,indvars(3))*q(i,j-1+1,indvars(3)))))/(param_float(4 + 5))

d2_dif_rhoet_dydy_0_0_ijm1nzp2 = -&
          0.5_wp*d2_dif_rhoet_dydy_0_0_ijm1nzp2_ijm1m1nzp2+&
          0.5_wp*d2_dif_rhoet_dydy_0_0_ijm1nzp2_ijm1p1nzp2

d2_dif_rhoet_dydy_0_0_ijm1nzp2 = d2_dif_rhoet_dydy_0_0_ijm1nzp2*param_float(2)

d2_dif_rhoet_dydy_0_0_ijp1nzp2_ijp1m1nzp2 = ((q(i,j+1-1,indvars(4))-&
                    0.5_wp*(q(i,j+1-1,indvars(2))*q(i,j+1-1,indvars(2))+&
                    q(i,j+1-1,indvars(3))*q(i,j+1-1,indvars(3)))))/(param_float(4 + 5))

d2_dif_rhoet_dydy_0_0_ijp1nzp2_ijp1p1nzp2 = ((q(i,j+1+1,indvars(4))-&
                    0.5_wp*(q(i,j+1+1,indvars(2))*q(i,j+1+1,indvars(2))+&
                    q(i,j+1+1,indvars(3))*q(i,j+1+1,indvars(3)))))/(param_float(4 + 5))

d2_dif_rhoet_dydy_0_0_ijp1nzp2 = -&
          0.5_wp*d2_dif_rhoet_dydy_0_0_ijp1nzp2_ijp1m1nzp2+&
          0.5_wp*d2_dif_rhoet_dydy_0_0_ijp1nzp2_ijp1p1nzp2

d2_dif_rhoet_dydy_0_0_ijp1nzp2 = d2_dif_rhoet_dydy_0_0_ijp1nzp2*param_float(2)

d1_dif_rhoet_dy_0_ijm1nzp2 = -param_float(2 + 5)*qst(i,j-1,indvarsst(11))*(d2_dif_rhoet_dydy_0_0_ijm1nzp2)-&
                    q(i,j-1,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(i,j-1,indvars(5))/1.0_wp)**3.0_wp/((q(i,j-1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i,j-1,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(i,j-1,indvarsst(11))*(d2_dif_rhou_dydy_0_0_ijm1nzp2)+&
                    qst(i,j-1,indvarsst(10))*(d2_dif_rhou_dydx_0_0_ijm1nzp2)))-&
                    q(i,j-1,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(i,j-1,indvars(5))/1.0_wp)**3.0_wp/((q(i,j-1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i,j-1,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(i,j-1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_ijm1nzp2)-&
                    2.0_wp/3.0_wp*(qst(i,j-1,indvarsst(10))*(d2_dif_rhov_dydx_0_0_ijm1nzp2)+&
                    qst(i,j-1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_ijm1nzp2))))

d1_dif_rhoet_dy_0_ijp1nzp2 = -param_float(2 + 5)*qst(i,j+1,indvarsst(11))*(d2_dif_rhoet_dydy_0_0_ijp1nzp2)-&
                    q(i,j+1,indvars(2))*((1.0_wp)*(1.0_wp+&
                    ((q(i,j+1,indvars(5))/1.0_wp)**3.0_wp/((q(i,j+1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i,j+1,indvars(5))/1.0_wp))*param_float(1 + 5)*(qst(i,j+1,indvarsst(11))*(d2_dif_rhou_dydy_0_0_ijp1nzp2)+&
                    qst(i,j+1,indvarsst(10))*(d2_dif_rhou_dydx_0_0_ijp1nzp2)))-&
                    q(i,j+1,indvars(3))*((1.0_wp)*(1.0_wp+&
                    ((q(i,j+1,indvars(5))/1.0_wp)**3.0_wp/((q(i,j+1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i,j+1,indvars(5))/1.0_wp))*param_float(1 + 5)*(2.0_wp*qst(i,j+1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_ijp1nzp2)-&
                    2.0_wp/3.0_wp*(qst(i,j+1,indvarsst(10))*(d2_dif_rhov_dydx_0_0_ijp1nzp2)+&
                    qst(i,j+1,indvarsst(11))*(d2_dif_rhov_dydy_0_0_ijp1nzp2))))

d1_dif_rhoet_dy_0_ijnzp2 = -&
          0.5_wp*d1_dif_rhoet_dy_0_ijm1nzp2+&
          0.5_wp*d1_dif_rhoet_dy_0_ijp1nzp2

d1_dif_rhoet_dy_0_ijnzp2 = d1_dif_rhoet_dy_0_ijnzp2*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None None 0 d(rho et)/dt *******
!                                                           
!***********************************************************


rhs(i,j,indvars(4)) = rhs(i,j,indvars(4))  -  ( qst(i,j,indvarsst(10))*(d1_dif_rhoet_dx_0_ijnzp2)+&
                    qst(i,j,indvarsst(11))*(d1_dif_rhoet_dy_0_ijnzp2) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer None None 0 d(rho nut)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! deltaxI*([-(ReI*(1.0_wp+chi)*sigmaI*deltaxI*({nut}_1x))]_1x)+deltayI*([-(ReI*(1.0_wp+chi)*sigmaI*deltayI*({nut}_1y))]_1y)-ReI*Cb2*sigmaI*((deltaxI)**2*([rho*nut]_1x)*([nut]_1x)+(deltayI)**2*([rho*nut]_1y)*([nut]_1y))-Cb1*(1.0_wp-ft2)*SS*rho*nut+ReI*(Cw1*fw-Cb1/k**2.0_wp*ft2)*rho*(nut/eta)**2.0_wp
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_dif_rhonut_dxdx_0_0_im1jnzp2_im1m1jnzp2 = q(i-1-1,j,indvars(5))

d2_dif_rhonut_dxdx_0_0_im1jnzp2_im1p1jnzp2 = q(i-1+1,j,indvars(5))

d2_dif_rhonut_dxdx_0_0_im1jnzp2 = -&
          0.5_wp*d2_dif_rhonut_dxdx_0_0_im1jnzp2_im1m1jnzp2+&
          0.5_wp*d2_dif_rhonut_dxdx_0_0_im1jnzp2_im1p1jnzp2

d2_dif_rhonut_dxdx_0_0_im1jnzp2 = d2_dif_rhonut_dxdx_0_0_im1jnzp2*param_float(1)

d2_dif_rhonut_dxdx_0_0_ip1jnzp2_ip1m1jnzp2 = q(i+1-1,j,indvars(5))

d2_dif_rhonut_dxdx_0_0_ip1jnzp2_ip1p1jnzp2 = q(i+1+1,j,indvars(5))

d2_dif_rhonut_dxdx_0_0_ip1jnzp2 = -&
          0.5_wp*d2_dif_rhonut_dxdx_0_0_ip1jnzp2_ip1m1jnzp2+&
          0.5_wp*d2_dif_rhonut_dxdx_0_0_ip1jnzp2_ip1p1jnzp2

d2_dif_rhonut_dxdx_0_0_ip1jnzp2 = d2_dif_rhonut_dxdx_0_0_ip1jnzp2*param_float(1)

d1_dif_rhonut_dx_1_im1jnzp2 = q(i-1,j,indvars(1))*q(i-1,j,indvars(5))

d1_dif_rhonut_dx_1_ip1jnzp2 = q(i+1,j,indvars(1))*q(i+1,j,indvars(5))

d1_dif_rhonut_dx_1_ijnzp2 = -&
          0.5_wp*d1_dif_rhonut_dx_1_im1jnzp2+&
          0.5_wp*d1_dif_rhonut_dx_1_ip1jnzp2

d1_dif_rhonut_dx_1_ijnzp2 = d1_dif_rhonut_dx_1_ijnzp2*param_float(1)

d1_dif_rhonut_dx_2_im1jnzp2 = q(i-1,j,indvars(5))

d1_dif_rhonut_dx_2_ip1jnzp2 = q(i+1,j,indvars(5))

d1_dif_rhonut_dx_2_ijnzp2 = -&
          0.5_wp*d1_dif_rhonut_dx_2_im1jnzp2+&
          0.5_wp*d1_dif_rhonut_dx_2_ip1jnzp2

d1_dif_rhonut_dx_2_ijnzp2 = d1_dif_rhonut_dx_2_ijnzp2*param_float(1)

d1_dif_rhonut_dx_0_im1jnzp2 = -(param_float(1 + 5)*(1.0_wp+&
                    (q(i-1,j,indvars(5))/1.0_wp))*param_float(18 + 5)*qst(i-1,j,indvarsst(10))*(d2_dif_rhonut_dxdx_0_0_im1jnzp2))

d1_dif_rhonut_dx_0_ip1jnzp2 = -(param_float(1 + 5)*(1.0_wp+&
                    (q(i+1,j,indvars(5))/1.0_wp))*param_float(18 + 5)*qst(i+1,j,indvarsst(10))*(d2_dif_rhonut_dxdx_0_0_ip1jnzp2))

d1_dif_rhonut_dx_0_ijnzp2 = -&
          0.5_wp*d1_dif_rhonut_dx_0_im1jnzp2+&
          0.5_wp*d1_dif_rhonut_dx_0_ip1jnzp2

d1_dif_rhonut_dx_0_ijnzp2 = d1_dif_rhonut_dx_0_ijnzp2*param_float(1)

d2_dif_rhonut_dydy_0_0_ijm1nzp2_ijm1m1nzp2 = q(i,j-1-1,indvars(5))

d2_dif_rhonut_dydy_0_0_ijm1nzp2_ijm1p1nzp2 = q(i,j-1+1,indvars(5))

d2_dif_rhonut_dydy_0_0_ijm1nzp2 = -&
          0.5_wp*d2_dif_rhonut_dydy_0_0_ijm1nzp2_ijm1m1nzp2+&
          0.5_wp*d2_dif_rhonut_dydy_0_0_ijm1nzp2_ijm1p1nzp2

d2_dif_rhonut_dydy_0_0_ijm1nzp2 = d2_dif_rhonut_dydy_0_0_ijm1nzp2*param_float(2)

d2_dif_rhonut_dydy_0_0_ijp1nzp2_ijp1m1nzp2 = q(i,j+1-1,indvars(5))

d2_dif_rhonut_dydy_0_0_ijp1nzp2_ijp1p1nzp2 = q(i,j+1+1,indvars(5))

d2_dif_rhonut_dydy_0_0_ijp1nzp2 = -&
          0.5_wp*d2_dif_rhonut_dydy_0_0_ijp1nzp2_ijp1m1nzp2+&
          0.5_wp*d2_dif_rhonut_dydy_0_0_ijp1nzp2_ijp1p1nzp2

d2_dif_rhonut_dydy_0_0_ijp1nzp2 = d2_dif_rhonut_dydy_0_0_ijp1nzp2*param_float(2)

d1_dif_rhonut_dy_1_ijm1nzp2 = q(i,j-1,indvars(1))*q(i,j-1,indvars(5))

d1_dif_rhonut_dy_1_ijp1nzp2 = q(i,j+1,indvars(1))*q(i,j+1,indvars(5))

d1_dif_rhonut_dy_1_ijnzp2 = -&
          0.5_wp*d1_dif_rhonut_dy_1_ijm1nzp2+&
          0.5_wp*d1_dif_rhonut_dy_1_ijp1nzp2

d1_dif_rhonut_dy_1_ijnzp2 = d1_dif_rhonut_dy_1_ijnzp2*param_float(2)

d1_dif_rhonut_dy_2_ijm1nzp2 = q(i,j-1,indvars(5))

d1_dif_rhonut_dy_2_ijp1nzp2 = q(i,j+1,indvars(5))

d1_dif_rhonut_dy_2_ijnzp2 = -&
          0.5_wp*d1_dif_rhonut_dy_2_ijm1nzp2+&
          0.5_wp*d1_dif_rhonut_dy_2_ijp1nzp2

d1_dif_rhonut_dy_2_ijnzp2 = d1_dif_rhonut_dy_2_ijnzp2*param_float(2)

d1_dif_rhonut_dy_0_ijm1nzp2 = -(param_float(1 + 5)*(1.0_wp+&
                    (q(i,j-1,indvars(5))/1.0_wp))*param_float(18 + 5)*qst(i,j-1,indvarsst(11))*(d2_dif_rhonut_dydy_0_0_ijm1nzp2))

d1_dif_rhonut_dy_0_ijp1nzp2 = -(param_float(1 + 5)*(1.0_wp+&
                    (q(i,j+1,indvars(5))/1.0_wp))*param_float(18 + 5)*qst(i,j+1,indvarsst(11))*(d2_dif_rhonut_dydy_0_0_ijp1nzp2))

d1_dif_rhonut_dy_0_ijnzp2 = -&
          0.5_wp*d1_dif_rhonut_dy_0_ijm1nzp2+&
          0.5_wp*d1_dif_rhonut_dy_0_ijp1nzp2

d1_dif_rhonut_dy_0_ijnzp2 = d1_dif_rhonut_dy_0_ijnzp2*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None None 0 d(rho nut)/dt ******
!                                                           
!***********************************************************


rhs(i,j,indvars(5)) = rhs(i,j,indvars(5))  -  ( qst(i,j,indvarsst(10))*(d1_dif_rhonut_dx_0_ijnzp2)+&
                    qst(i,j,indvarsst(11))*(d1_dif_rhonut_dy_0_ijnzp2)-&
                    param_float(1 + 5)*param_float(7 + 5)*param_float(18 + 5)*((qst(i,j,indvarsst(10)))**2*(d1_dif_rhonut_dx_1_ijnzp2)*(d1_dif_rhonut_dx_2_ijnzp2)+&
                    (qst(i,j,indvarsst(11)))**2*(d1_dif_rhonut_dy_1_ijnzp2)*(d1_dif_rhonut_dy_2_ijnzp2))-&
                    param_float(6 + 5)*(1.0_wp-&
                    (param_float(16 + 5)*exp(-&
                    param_float(17 + 5)*(q(i,j,indvars(5))/1.0_wp)**2.0_wp)))*qst(i,j,indvarsst(13))*q(i,j,indvars(1))*q(i,j,indvars(5))+&
                    param_float(1 + 5)*(param_float(10 + 5)*qst(i,j,indvarsst(16))-&
                    param_float(6 + 5)/param_float(9 + 5)**2.0_wp*(param_float(16 + 5)*exp(-&
                    param_float(17 + 5)*(q(i,j,indvars(5))/1.0_wp)**2.0_wp)))*q(i,j,indvars(1))*(q(i,j,indvars(5))/qst(i,j,indvarsst(2)))**2.0_wp ) 

     enddo
   enddo
