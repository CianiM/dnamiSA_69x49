

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
! building source terms in RHS for layer 0 None None stemp *
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (((0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))**2+(0.5_wp*(deltaxI*([v]_1x)-deltayI*([u]_1y)))**2)*2)**0.5
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_nxp4p0p0jk = q(nx+4+0+0,j,indvars(3))

d1_stemp_dx_0_nxp4p0m1jk = q(nx+4+0-1,j,indvars(3))

d1_stemp_dx_0_nxp4p0m2jk = q(nx+4+0-2,j,indvars(3))

d1_stemp_dx_0_nxp4p0jk = 1.5_wp*d1_stemp_dx_0_nxp4p0p0jk-&
          2.0_wp*d1_stemp_dx_0_nxp4p0m1jk+&
          0.5_wp*d1_stemp_dx_0_nxp4p0m2jk

d1_stemp_dx_0_nxp4p0jk = d1_stemp_dx_0_nxp4p0jk*param_float(1)

d1_stemp_dy_0_nxp4p0jm1k = q(nx+4+0,j-1,indvars(2))

d1_stemp_dy_0_nxp4p0jp1k = q(nx+4+0,j+1,indvars(2))

d1_stemp_dy_0_nxp4p0jk = -&
          0.5_wp*d1_stemp_dy_0_nxp4p0jm1k+&
          0.5_wp*d1_stemp_dy_0_nxp4p0jp1k

d1_stemp_dy_0_nxp4p0jk = d1_stemp_dy_0_nxp4p0jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None stemp **************
!                                                           
!***********************************************************


qst(nx+4+0,j,indvarsst(4)) =  (((0.5_wp*(qst(nx+4+0,j,indvarsst(11))*(d1_stemp_dy_0_nxp4p0jk)-&
                    qst(nx+4+0,j,indvarsst(10))*(d1_stemp_dx_0_nxp4p0jk)))**2+&
                    (0.5_wp*(qst(nx+4+0,j,indvarsst(10))*(d1_stemp_dx_0_nxp4p0jk)-&
                    qst(nx+4+0,j,indvarsst(11))*(d1_stemp_dy_0_nxp4p0jk)))**2)*2)**0.5



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None SS ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (stemp+ReI*nut/(k**2*eta**2))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None SS *****************
!                                                           
!***********************************************************


qst(nx+4+0,j,indvarsst(12)) =  (qst(nx+4+0,j,indvarsst(4))+&
                    param_float(1 + 5)*q(nx+4+0,j,indvars(5))/(param_float(9 + 5)**2*qst(nx+4+0,j,indvarsst(2))**2))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None tau_wall 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (visc_t*([u]_1y)*deltayI)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None tau_wall ***********
!                                                           
!***********************************************************


qst(nx+4+0,j,indvarsst(13)) =  (((1+&
                    param_float(21 + 5))/(((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+4+0,j,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0,j,indvars(1)))**3/((q(nx+4+0,j,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0,j,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+4+0,j,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4+0,j,indvars(4))-&
                    0.5_wp*(q(nx+4+0,j,indvars(2))*q(nx+4+0,j,indvars(2))+&
                    q(nx+4+0,j,indvars(3))*q(nx+4+0,j,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4+0,j,indvars(1))))*param_float(1 + 5)*(d1_stemp_dy_0_nxp4p0jk)*qst(nx+4+0,j,indvarsst(11)))

     enddo


!***********************************************************
!                                                           
! Start building layers for BC : imax None None ************
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
! building source terms in RHS for layer 1 None None stemp *
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (((0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))**2+(0.5_wp*(deltaxI*([v]_1x)-deltayI*([u]_1y)))**2)*2)**0.5
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_nxp4m1m1jk = q(nx+4-1-1,j,indvars(3))

d1_stemp_dx_0_nxp4m1p1jk = q(nx+4-1+1,j,indvars(3))

d1_stemp_dx_0_nxp4m1jk = -&
          0.5_wp*d1_stemp_dx_0_nxp4m1m1jk+&
          0.5_wp*d1_stemp_dx_0_nxp4m1p1jk

d1_stemp_dx_0_nxp4m1jk = d1_stemp_dx_0_nxp4m1jk*param_float(1)

d1_stemp_dy_0_nxp4m1jm1k = q(nx+4-1,j-1,indvars(2))

d1_stemp_dy_0_nxp4m1jp1k = q(nx+4-1,j+1,indvars(2))

d1_stemp_dy_0_nxp4m1jk = -&
          0.5_wp*d1_stemp_dy_0_nxp4m1jm1k+&
          0.5_wp*d1_stemp_dy_0_nxp4m1jp1k

d1_stemp_dy_0_nxp4m1jk = d1_stemp_dy_0_nxp4m1jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None stemp **************
!                                                           
!***********************************************************


qst(nx+4-1,j,indvarsst(4)) =  (((0.5_wp*(qst(nx+4-1,j,indvarsst(11))*(d1_stemp_dy_0_nxp4m1jk)-&
                    qst(nx+4-1,j,indvarsst(10))*(d1_stemp_dx_0_nxp4m1jk)))**2+&
                    (0.5_wp*(qst(nx+4-1,j,indvarsst(10))*(d1_stemp_dx_0_nxp4m1jk)-&
                    qst(nx+4-1,j,indvarsst(11))*(d1_stemp_dy_0_nxp4m1jk)))**2)*2)**0.5



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None SS ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (stemp+ReI*nut/(k**2*eta**2))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None SS *****************
!                                                           
!***********************************************************


qst(nx+4-1,j,indvarsst(12)) =  (qst(nx+4-1,j,indvarsst(4))+&
                    param_float(1 + 5)*q(nx+4-1,j,indvars(5))/(param_float(9 + 5)**2*qst(nx+4-1,j,indvarsst(2))**2))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None tau_wall 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (visc_t*([u]_1y)*deltayI)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None tau_wall ***********
!                                                           
!***********************************************************


qst(nx+4-1,j,indvarsst(13)) =  (((1+&
                    param_float(21 + 5))/(((q(nx+4-1,j,indvars(4))-&
                    0.5_wp*(q(nx+4-1,j,indvars(2))*q(nx+4-1,j,indvars(2))+&
                    q(nx+4-1,j,indvars(3))*q(nx+4-1,j,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4-1,j,indvars(4))-&
                    0.5_wp*(q(nx+4-1,j,indvars(2))*q(nx+4-1,j,indvars(2))+&
                    q(nx+4-1,j,indvars(3))*q(nx+4-1,j,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+4-1,j,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4-1,j,indvars(4))-&
                    0.5_wp*(q(nx+4-1,j,indvars(2))*q(nx+4-1,j,indvars(2))+&
                    q(nx+4-1,j,indvars(3))*q(nx+4-1,j,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4-1,j,indvars(4))-&
                    0.5_wp*(q(nx+4-1,j,indvars(2))*q(nx+4-1,j,indvars(2))+&
                    q(nx+4-1,j,indvars(3))*q(nx+4-1,j,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4-1,j,indvars(1)))**3/((q(nx+4-1,j,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4-1,j,indvars(4))-&
                    0.5_wp*(q(nx+4-1,j,indvars(2))*q(nx+4-1,j,indvars(2))+&
                    q(nx+4-1,j,indvars(3))*q(nx+4-1,j,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4-1,j,indvars(4))-&
                    0.5_wp*(q(nx+4-1,j,indvars(2))*q(nx+4-1,j,indvars(2))+&
                    q(nx+4-1,j,indvars(3))*q(nx+4-1,j,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4-1,j,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+4-1,j,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4-1,j,indvars(4))-&
                    0.5_wp*(q(nx+4-1,j,indvars(2))*q(nx+4-1,j,indvars(2))+&
                    q(nx+4-1,j,indvars(3))*q(nx+4-1,j,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4-1,j,indvars(4))-&
                    0.5_wp*(q(nx+4-1,j,indvars(2))*q(nx+4-1,j,indvars(2))+&
                    q(nx+4-1,j,indvars(3))*q(nx+4-1,j,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4-1,j,indvars(1))))*param_float(1 + 5)*(d1_stemp_dy_0_nxp4m1jk)*qst(nx+4-1,j,indvarsst(11)))

     enddo


!***********************************************************
!                                                           
! Start building layers for BC : imax None None ************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 2 None None ************************************
!                                                           
!***********************************************************


     do j=idloop(3),idloop(4) 


!***********************************************************
!                                                           
! building source terms in RHS for layer 2 None None stemp *
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (((0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))**2+(0.5_wp*(deltaxI*([v]_1x)-deltayI*([u]_1y)))**2)*2)**0.5
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_nxp4m2m1jk = q(nx+4-2-1,j,indvars(3))

d1_stemp_dx_0_nxp4m2p1jk = q(nx+4-2+1,j,indvars(3))

d1_stemp_dx_0_nxp4m2jk = -&
          0.5_wp*d1_stemp_dx_0_nxp4m2m1jk+&
          0.5_wp*d1_stemp_dx_0_nxp4m2p1jk

d1_stemp_dx_0_nxp4m2jk = d1_stemp_dx_0_nxp4m2jk*param_float(1)

d1_stemp_dy_0_nxp4m2jm1k = q(nx+4-2,j-1,indvars(2))

d1_stemp_dy_0_nxp4m2jp1k = q(nx+4-2,j+1,indvars(2))

d1_stemp_dy_0_nxp4m2jk = -&
          0.5_wp*d1_stemp_dy_0_nxp4m2jm1k+&
          0.5_wp*d1_stemp_dy_0_nxp4m2jp1k

d1_stemp_dy_0_nxp4m2jk = d1_stemp_dy_0_nxp4m2jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 2 None None stemp **************
!                                                           
!***********************************************************


qst(nx+4-2,j,indvarsst(4)) =  (((0.5_wp*(qst(nx+4-2,j,indvarsst(11))*(d1_stemp_dy_0_nxp4m2jk)-&
                    qst(nx+4-2,j,indvarsst(10))*(d1_stemp_dx_0_nxp4m2jk)))**2+&
                    (0.5_wp*(qst(nx+4-2,j,indvarsst(10))*(d1_stemp_dx_0_nxp4m2jk)-&
                    qst(nx+4-2,j,indvarsst(11))*(d1_stemp_dy_0_nxp4m2jk)))**2)*2)**0.5



!***********************************************************
!                                                           
! building source terms in RHS for layer 2 None None SS ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (stemp+ReI*nut/(k**2*eta**2))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 2 None None SS *****************
!                                                           
!***********************************************************


qst(nx+4-2,j,indvarsst(12)) =  (qst(nx+4-2,j,indvarsst(4))+&
                    param_float(1 + 5)*q(nx+4-2,j,indvars(5))/(param_float(9 + 5)**2*qst(nx+4-2,j,indvarsst(2))**2))



!***********************************************************
!                                                           
! building source terms in RHS for layer 2 None None tau_wall 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (visc_t*([u]_1y)*deltayI)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 2 None None tau_wall ***********
!                                                           
!***********************************************************


qst(nx+4-2,j,indvarsst(13)) =  (((1+&
                    param_float(21 + 5))/(((q(nx+4-2,j,indvars(4))-&
                    0.5_wp*(q(nx+4-2,j,indvars(2))*q(nx+4-2,j,indvars(2))+&
                    q(nx+4-2,j,indvars(3))*q(nx+4-2,j,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4-2,j,indvars(4))-&
                    0.5_wp*(q(nx+4-2,j,indvars(2))*q(nx+4-2,j,indvars(2))+&
                    q(nx+4-2,j,indvars(3))*q(nx+4-2,j,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+4-2,j,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4-2,j,indvars(4))-&
                    0.5_wp*(q(nx+4-2,j,indvars(2))*q(nx+4-2,j,indvars(2))+&
                    q(nx+4-2,j,indvars(3))*q(nx+4-2,j,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4-2,j,indvars(4))-&
                    0.5_wp*(q(nx+4-2,j,indvars(2))*q(nx+4-2,j,indvars(2))+&
                    q(nx+4-2,j,indvars(3))*q(nx+4-2,j,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4-2,j,indvars(1)))**3/((q(nx+4-2,j,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4-2,j,indvars(4))-&
                    0.5_wp*(q(nx+4-2,j,indvars(2))*q(nx+4-2,j,indvars(2))+&
                    q(nx+4-2,j,indvars(3))*q(nx+4-2,j,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4-2,j,indvars(4))-&
                    0.5_wp*(q(nx+4-2,j,indvars(2))*q(nx+4-2,j,indvars(2))+&
                    q(nx+4-2,j,indvars(3))*q(nx+4-2,j,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4-2,j,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+4-2,j,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4-2,j,indvars(4))-&
                    0.5_wp*(q(nx+4-2,j,indvars(2))*q(nx+4-2,j,indvars(2))+&
                    q(nx+4-2,j,indvars(3))*q(nx+4-2,j,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4-2,j,indvars(4))-&
                    0.5_wp*(q(nx+4-2,j,indvars(2))*q(nx+4-2,j,indvars(2))+&
                    q(nx+4-2,j,indvars(3))*q(nx+4-2,j,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4-2,j,indvars(1))))*param_float(1 + 5)*(d1_stemp_dy_0_nxp4m2jk)*qst(nx+4-2,j,indvarsst(11)))

     enddo


!***********************************************************
!                                                           
! Start building layers for BC : imax None None ************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 3 None None ************************************
!                                                           
!***********************************************************


     do j=idloop(3),idloop(4) 


!***********************************************************
!                                                           
! building source terms in RHS for layer 3 None None stemp *
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (((0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))**2+(0.5_wp*(deltaxI*([v]_1x)-deltayI*([u]_1y)))**2)*2)**0.5
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_nxp4m3m1jk = q(nx+4-3-1,j,indvars(3))

d1_stemp_dx_0_nxp4m3p1jk = q(nx+4-3+1,j,indvars(3))

d1_stemp_dx_0_nxp4m3jk = -&
          0.5_wp*d1_stemp_dx_0_nxp4m3m1jk+&
          0.5_wp*d1_stemp_dx_0_nxp4m3p1jk

d1_stemp_dx_0_nxp4m3jk = d1_stemp_dx_0_nxp4m3jk*param_float(1)

d1_stemp_dy_0_nxp4m3jm1k = q(nx+4-3,j-1,indvars(2))

d1_stemp_dy_0_nxp4m3jp1k = q(nx+4-3,j+1,indvars(2))

d1_stemp_dy_0_nxp4m3jk = -&
          0.5_wp*d1_stemp_dy_0_nxp4m3jm1k+&
          0.5_wp*d1_stemp_dy_0_nxp4m3jp1k

d1_stemp_dy_0_nxp4m3jk = d1_stemp_dy_0_nxp4m3jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 3 None None stemp **************
!                                                           
!***********************************************************


qst(nx+4-3,j,indvarsst(4)) =  (((0.5_wp*(qst(nx+4-3,j,indvarsst(11))*(d1_stemp_dy_0_nxp4m3jk)-&
                    qst(nx+4-3,j,indvarsst(10))*(d1_stemp_dx_0_nxp4m3jk)))**2+&
                    (0.5_wp*(qst(nx+4-3,j,indvarsst(10))*(d1_stemp_dx_0_nxp4m3jk)-&
                    qst(nx+4-3,j,indvarsst(11))*(d1_stemp_dy_0_nxp4m3jk)))**2)*2)**0.5



!***********************************************************
!                                                           
! building source terms in RHS for layer 3 None None SS ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (stemp+ReI*nut/(k**2*eta**2))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 3 None None SS *****************
!                                                           
!***********************************************************


qst(nx+4-3,j,indvarsst(12)) =  (qst(nx+4-3,j,indvarsst(4))+&
                    param_float(1 + 5)*q(nx+4-3,j,indvars(5))/(param_float(9 + 5)**2*qst(nx+4-3,j,indvarsst(2))**2))



!***********************************************************
!                                                           
! building source terms in RHS for layer 3 None None tau_wall 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (visc_t*([u]_1y)*deltayI)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 3 None None tau_wall ***********
!                                                           
!***********************************************************


qst(nx+4-3,j,indvarsst(13)) =  (((1+&
                    param_float(21 + 5))/(((q(nx+4-3,j,indvars(4))-&
                    0.5_wp*(q(nx+4-3,j,indvars(2))*q(nx+4-3,j,indvars(2))+&
                    q(nx+4-3,j,indvars(3))*q(nx+4-3,j,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4-3,j,indvars(4))-&
                    0.5_wp*(q(nx+4-3,j,indvars(2))*q(nx+4-3,j,indvars(2))+&
                    q(nx+4-3,j,indvars(3))*q(nx+4-3,j,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(nx+4-3,j,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4-3,j,indvars(4))-&
                    0.5_wp*(q(nx+4-3,j,indvars(2))*q(nx+4-3,j,indvars(2))+&
                    q(nx+4-3,j,indvars(3))*q(nx+4-3,j,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4-3,j,indvars(4))-&
                    0.5_wp*(q(nx+4-3,j,indvars(2))*q(nx+4-3,j,indvars(2))+&
                    q(nx+4-3,j,indvars(3))*q(nx+4-3,j,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4-3,j,indvars(1)))**3/((q(nx+4-3,j,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4-3,j,indvars(4))-&
                    0.5_wp*(q(nx+4-3,j,indvars(2))*q(nx+4-3,j,indvars(2))+&
                    q(nx+4-3,j,indvars(3))*q(nx+4-3,j,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4-3,j,indvars(4))-&
                    0.5_wp*(q(nx+4-3,j,indvars(2))*q(nx+4-3,j,indvars(2))+&
                    q(nx+4-3,j,indvars(3))*q(nx+4-3,j,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4-3,j,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(nx+4-3,j,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(nx+4-3,j,indvars(4))-&
                    0.5_wp*(q(nx+4-3,j,indvars(2))*q(nx+4-3,j,indvars(2))+&
                    q(nx+4-3,j,indvars(3))*q(nx+4-3,j,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(nx+4-3,j,indvars(4))-&
                    0.5_wp*(q(nx+4-3,j,indvars(2))*q(nx+4-3,j,indvars(2))+&
                    q(nx+4-3,j,indvars(3))*q(nx+4-3,j,indvars(3)))))/param_float(4 + 5)**1.5*q(nx+4-3,j,indvars(1))))*param_float(1 + 5)*(d1_stemp_dy_0_nxp4m3jk)*qst(nx+4-3,j,indvarsst(11)))

     enddo
