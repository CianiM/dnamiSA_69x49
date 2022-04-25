

!***********************************************************
!                                                           
! Start building layers for BC : None j1 None **************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: None 0 None ************************************
!                                                           
!***********************************************************


 
      do i=idloop(1),idloop(2) 


!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None stemp *
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (((0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))**2+(0.5_wp*(deltaxI*([v]_1x)-deltayI*([u]_1y)))**2)*2)**0.5
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_im11m4p0k = q(i-1,1-4+0,indvars(3))

d1_stemp_dx_0_ip11m4p0k = q(i+1,1-4+0,indvars(3))

d1_stemp_dx_0_i1m4p0k = -&
          0.5_wp*d1_stemp_dx_0_im11m4p0k+&
          0.5_wp*d1_stemp_dx_0_ip11m4p0k

d1_stemp_dx_0_i1m4p0k = d1_stemp_dx_0_i1m4p0k*param_float(1)

d1_stemp_dy_0_i1m4p0p0k = q(i,1-4+0+0,indvars(2))

d1_stemp_dy_0_i1m4p0p1k = q(i,1-4+0+1,indvars(2))

d1_stemp_dy_0_i1m4p0p2k = q(i,1-4+0+2,indvars(2))

d1_stemp_dy_0_i1m4p0k = -&
          1.5_wp*d1_stemp_dy_0_i1m4p0p0k+&
          2.0_wp*d1_stemp_dy_0_i1m4p0p1k-&
          0.5_wp*d1_stemp_dy_0_i1m4p0p2k

d1_stemp_dy_0_i1m4p0k = d1_stemp_dy_0_i1m4p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None stemp **************
!                                                           
!***********************************************************


qst(i,1-4+0,indvarsst(4)) =  (((0.5_wp*(qst(i,1-4+0,indvarsst(11))*(d1_stemp_dy_0_i1m4p0k)-&
                    qst(i,1-4+0,indvarsst(10))*(d1_stemp_dx_0_i1m4p0k)))**2+&
                    (0.5_wp*(qst(i,1-4+0,indvarsst(10))*(d1_stemp_dx_0_i1m4p0k)-&
                    qst(i,1-4+0,indvarsst(11))*(d1_stemp_dy_0_i1m4p0k)))**2)*2)**0.5



!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None SS ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (stemp+ReI*nut/(k**2*eta**2))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None SS *****************
!                                                           
!***********************************************************


qst(i,1-4+0,indvarsst(12)) =  (qst(i,1-4+0,indvarsst(4))+&
                    param_float(1 + 5)*q(i,1-4+0,indvars(5))/(param_float(9 + 5)**2*qst(i,1-4+0,indvarsst(2))**2))



!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None tau_wall 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (visc_t*([u]_1y)*deltayI)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None tau_wall ***********
!                                                           
!***********************************************************


qst(i,1-4+0,indvarsst(13)) =  (((1+&
                    param_float(21 + 5))/(((q(i,1-4+0,indvars(4))-&
                    0.5_wp*(q(i,1-4+0,indvars(2))*q(i,1-4+0,indvars(2))+&
                    q(i,1-4+0,indvars(3))*q(i,1-4+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,1-4+0,indvars(4))-&
                    0.5_wp*(q(i,1-4+0,indvars(2))*q(i,1-4+0,indvars(2))+&
                    q(i,1-4+0,indvars(3))*q(i,1-4+0,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(i,1-4+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,1-4+0,indvars(4))-&
                    0.5_wp*(q(i,1-4+0,indvars(2))*q(i,1-4+0,indvars(2))+&
                    q(i,1-4+0,indvars(3))*q(i,1-4+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,1-4+0,indvars(4))-&
                    0.5_wp*(q(i,1-4+0,indvars(2))*q(i,1-4+0,indvars(2))+&
                    q(i,1-4+0,indvars(3))*q(i,1-4+0,indvars(3)))))/param_float(4 + 5)**1.5*q(i,1-4+0,indvars(1)))**3/((q(i,1-4+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,1-4+0,indvars(4))-&
                    0.5_wp*(q(i,1-4+0,indvars(2))*q(i,1-4+0,indvars(2))+&
                    q(i,1-4+0,indvars(3))*q(i,1-4+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,1-4+0,indvars(4))-&
                    0.5_wp*(q(i,1-4+0,indvars(2))*q(i,1-4+0,indvars(2))+&
                    q(i,1-4+0,indvars(3))*q(i,1-4+0,indvars(3)))))/param_float(4 + 5)**1.5*q(i,1-4+0,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(i,1-4+0,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,1-4+0,indvars(4))-&
                    0.5_wp*(q(i,1-4+0,indvars(2))*q(i,1-4+0,indvars(2))+&
                    q(i,1-4+0,indvars(3))*q(i,1-4+0,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,1-4+0,indvars(4))-&
                    0.5_wp*(q(i,1-4+0,indvars(2))*q(i,1-4+0,indvars(2))+&
                    q(i,1-4+0,indvars(3))*q(i,1-4+0,indvars(3)))))/param_float(4 + 5)**1.5*q(i,1-4+0,indvars(1))))*param_float(1 + 5)*(d1_stemp_dy_0_i1m4p0k)*qst(i,1-4+0,indvarsst(11)))

   enddo


!***********************************************************
!                                                           
! Start building layers for BC : None j1 None **************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: None 1 None ************************************
!                                                           
!***********************************************************


 
      do i=idloop(1),idloop(2) 


!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None stemp *
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (((0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))**2+(0.5_wp*(deltaxI*([v]_1x)-deltayI*([u]_1y)))**2)*2)**0.5
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_im11m4p1k = q(i-1,1-4+1,indvars(3))

d1_stemp_dx_0_ip11m4p1k = q(i+1,1-4+1,indvars(3))

d1_stemp_dx_0_i1m4p1k = -&
          0.5_wp*d1_stemp_dx_0_im11m4p1k+&
          0.5_wp*d1_stemp_dx_0_ip11m4p1k

d1_stemp_dx_0_i1m4p1k = d1_stemp_dx_0_i1m4p1k*param_float(1)

d1_stemp_dy_0_i1m4p1m1k = q(i,1-4+1-1,indvars(2))

d1_stemp_dy_0_i1m4p1p1k = q(i,1-4+1+1,indvars(2))

d1_stemp_dy_0_i1m4p1k = -&
          0.5_wp*d1_stemp_dy_0_i1m4p1m1k+&
          0.5_wp*d1_stemp_dy_0_i1m4p1p1k

d1_stemp_dy_0_i1m4p1k = d1_stemp_dy_0_i1m4p1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None stemp **************
!                                                           
!***********************************************************


qst(i,1-4+1,indvarsst(4)) =  (((0.5_wp*(qst(i,1-4+1,indvarsst(11))*(d1_stemp_dy_0_i1m4p1k)-&
                    qst(i,1-4+1,indvarsst(10))*(d1_stemp_dx_0_i1m4p1k)))**2+&
                    (0.5_wp*(qst(i,1-4+1,indvarsst(10))*(d1_stemp_dx_0_i1m4p1k)-&
                    qst(i,1-4+1,indvarsst(11))*(d1_stemp_dy_0_i1m4p1k)))**2)*2)**0.5



!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None SS ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (stemp+ReI*nut/(k**2*eta**2))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None SS *****************
!                                                           
!***********************************************************


qst(i,1-4+1,indvarsst(12)) =  (qst(i,1-4+1,indvarsst(4))+&
                    param_float(1 + 5)*q(i,1-4+1,indvars(5))/(param_float(9 + 5)**2*qst(i,1-4+1,indvarsst(2))**2))



!***********************************************************
!                                                           
! building source terms in RHS for layer None 1 None tau_wall 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (visc_t*([u]_1y)*deltayI)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 1 None tau_wall ***********
!                                                           
!***********************************************************


qst(i,1-4+1,indvarsst(13)) =  (((1+&
                    param_float(21 + 5))/(((q(i,1-4+1,indvars(4))-&
                    0.5_wp*(q(i,1-4+1,indvars(2))*q(i,1-4+1,indvars(2))+&
                    q(i,1-4+1,indvars(3))*q(i,1-4+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,1-4+1,indvars(4))-&
                    0.5_wp*(q(i,1-4+1,indvars(2))*q(i,1-4+1,indvars(2))+&
                    q(i,1-4+1,indvars(3))*q(i,1-4+1,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(i,1-4+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,1-4+1,indvars(4))-&
                    0.5_wp*(q(i,1-4+1,indvars(2))*q(i,1-4+1,indvars(2))+&
                    q(i,1-4+1,indvars(3))*q(i,1-4+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,1-4+1,indvars(4))-&
                    0.5_wp*(q(i,1-4+1,indvars(2))*q(i,1-4+1,indvars(2))+&
                    q(i,1-4+1,indvars(3))*q(i,1-4+1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,1-4+1,indvars(1)))**3/((q(i,1-4+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,1-4+1,indvars(4))-&
                    0.5_wp*(q(i,1-4+1,indvars(2))*q(i,1-4+1,indvars(2))+&
                    q(i,1-4+1,indvars(3))*q(i,1-4+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,1-4+1,indvars(4))-&
                    0.5_wp*(q(i,1-4+1,indvars(2))*q(i,1-4+1,indvars(2))+&
                    q(i,1-4+1,indvars(3))*q(i,1-4+1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,1-4+1,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(i,1-4+1,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,1-4+1,indvars(4))-&
                    0.5_wp*(q(i,1-4+1,indvars(2))*q(i,1-4+1,indvars(2))+&
                    q(i,1-4+1,indvars(3))*q(i,1-4+1,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,1-4+1,indvars(4))-&
                    0.5_wp*(q(i,1-4+1,indvars(2))*q(i,1-4+1,indvars(2))+&
                    q(i,1-4+1,indvars(3))*q(i,1-4+1,indvars(3)))))/param_float(4 + 5)**1.5*q(i,1-4+1,indvars(1))))*param_float(1 + 5)*(d1_stemp_dy_0_i1m4p1k)*qst(i,1-4+1,indvarsst(11)))

   enddo


!***********************************************************
!                                                           
! Start building layers for BC : None j1 None **************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: None 2 None ************************************
!                                                           
!***********************************************************


 
      do i=idloop(1),idloop(2) 


!***********************************************************
!                                                           
! building source terms in RHS for layer None 2 None stemp *
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (((0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))**2+(0.5_wp*(deltaxI*([v]_1x)-deltayI*([u]_1y)))**2)*2)**0.5
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_im11m4p2k = q(i-1,1-4+2,indvars(3))

d1_stemp_dx_0_ip11m4p2k = q(i+1,1-4+2,indvars(3))

d1_stemp_dx_0_i1m4p2k = -&
          0.5_wp*d1_stemp_dx_0_im11m4p2k+&
          0.5_wp*d1_stemp_dx_0_ip11m4p2k

d1_stemp_dx_0_i1m4p2k = d1_stemp_dx_0_i1m4p2k*param_float(1)

d1_stemp_dy_0_i1m4p2m1k = q(i,1-4+2-1,indvars(2))

d1_stemp_dy_0_i1m4p2p1k = q(i,1-4+2+1,indvars(2))

d1_stemp_dy_0_i1m4p2k = -&
          0.5_wp*d1_stemp_dy_0_i1m4p2m1k+&
          0.5_wp*d1_stemp_dy_0_i1m4p2p1k

d1_stemp_dy_0_i1m4p2k = d1_stemp_dy_0_i1m4p2k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None 2 None stemp **************
!                                                           
!***********************************************************


qst(i,1-4+2,indvarsst(4)) =  (((0.5_wp*(qst(i,1-4+2,indvarsst(11))*(d1_stemp_dy_0_i1m4p2k)-&
                    qst(i,1-4+2,indvarsst(10))*(d1_stemp_dx_0_i1m4p2k)))**2+&
                    (0.5_wp*(qst(i,1-4+2,indvarsst(10))*(d1_stemp_dx_0_i1m4p2k)-&
                    qst(i,1-4+2,indvarsst(11))*(d1_stemp_dy_0_i1m4p2k)))**2)*2)**0.5



!***********************************************************
!                                                           
! building source terms in RHS for layer None 2 None SS ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (stemp+ReI*nut/(k**2*eta**2))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 2 None SS *****************
!                                                           
!***********************************************************


qst(i,1-4+2,indvarsst(12)) =  (qst(i,1-4+2,indvarsst(4))+&
                    param_float(1 + 5)*q(i,1-4+2,indvars(5))/(param_float(9 + 5)**2*qst(i,1-4+2,indvarsst(2))**2))



!***********************************************************
!                                                           
! building source terms in RHS for layer None 2 None tau_wall 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (visc_t*([u]_1y)*deltayI)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 2 None tau_wall ***********
!                                                           
!***********************************************************


qst(i,1-4+2,indvarsst(13)) =  (((1+&
                    param_float(21 + 5))/(((q(i,1-4+2,indvars(4))-&
                    0.5_wp*(q(i,1-4+2,indvars(2))*q(i,1-4+2,indvars(2))+&
                    q(i,1-4+2,indvars(3))*q(i,1-4+2,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,1-4+2,indvars(4))-&
                    0.5_wp*(q(i,1-4+2,indvars(2))*q(i,1-4+2,indvars(2))+&
                    q(i,1-4+2,indvars(3))*q(i,1-4+2,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(i,1-4+2,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,1-4+2,indvars(4))-&
                    0.5_wp*(q(i,1-4+2,indvars(2))*q(i,1-4+2,indvars(2))+&
                    q(i,1-4+2,indvars(3))*q(i,1-4+2,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,1-4+2,indvars(4))-&
                    0.5_wp*(q(i,1-4+2,indvars(2))*q(i,1-4+2,indvars(2))+&
                    q(i,1-4+2,indvars(3))*q(i,1-4+2,indvars(3)))))/param_float(4 + 5)**1.5*q(i,1-4+2,indvars(1)))**3/((q(i,1-4+2,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,1-4+2,indvars(4))-&
                    0.5_wp*(q(i,1-4+2,indvars(2))*q(i,1-4+2,indvars(2))+&
                    q(i,1-4+2,indvars(3))*q(i,1-4+2,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,1-4+2,indvars(4))-&
                    0.5_wp*(q(i,1-4+2,indvars(2))*q(i,1-4+2,indvars(2))+&
                    q(i,1-4+2,indvars(3))*q(i,1-4+2,indvars(3)))))/param_float(4 + 5)**1.5*q(i,1-4+2,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(i,1-4+2,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,1-4+2,indvars(4))-&
                    0.5_wp*(q(i,1-4+2,indvars(2))*q(i,1-4+2,indvars(2))+&
                    q(i,1-4+2,indvars(3))*q(i,1-4+2,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,1-4+2,indvars(4))-&
                    0.5_wp*(q(i,1-4+2,indvars(2))*q(i,1-4+2,indvars(2))+&
                    q(i,1-4+2,indvars(3))*q(i,1-4+2,indvars(3)))))/param_float(4 + 5)**1.5*q(i,1-4+2,indvars(1))))*param_float(1 + 5)*(d1_stemp_dy_0_i1m4p2k)*qst(i,1-4+2,indvarsst(11)))

   enddo


!***********************************************************
!                                                           
! Start building layers for BC : None j1 None **************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: None 3 None ************************************
!                                                           
!***********************************************************


 
      do i=idloop(1),idloop(2) 


!***********************************************************
!                                                           
! building source terms in RHS for layer None 3 None stemp *
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (((0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))**2+(0.5_wp*(deltaxI*([v]_1x)-deltayI*([u]_1y)))**2)*2)**0.5
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_im11m4p3k = q(i-1,1-4+3,indvars(3))

d1_stemp_dx_0_ip11m4p3k = q(i+1,1-4+3,indvars(3))

d1_stemp_dx_0_i1m4p3k = -&
          0.5_wp*d1_stemp_dx_0_im11m4p3k+&
          0.5_wp*d1_stemp_dx_0_ip11m4p3k

d1_stemp_dx_0_i1m4p3k = d1_stemp_dx_0_i1m4p3k*param_float(1)

d1_stemp_dy_0_i1m4p3m1k = q(i,1-4+3-1,indvars(2))

d1_stemp_dy_0_i1m4p3p1k = q(i,1-4+3+1,indvars(2))

d1_stemp_dy_0_i1m4p3k = -&
          0.5_wp*d1_stemp_dy_0_i1m4p3m1k+&
          0.5_wp*d1_stemp_dy_0_i1m4p3p1k

d1_stemp_dy_0_i1m4p3k = d1_stemp_dy_0_i1m4p3k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None 3 None stemp **************
!                                                           
!***********************************************************


qst(i,1-4+3,indvarsst(4)) =  (((0.5_wp*(qst(i,1-4+3,indvarsst(11))*(d1_stemp_dy_0_i1m4p3k)-&
                    qst(i,1-4+3,indvarsst(10))*(d1_stemp_dx_0_i1m4p3k)))**2+&
                    (0.5_wp*(qst(i,1-4+3,indvarsst(10))*(d1_stemp_dx_0_i1m4p3k)-&
                    qst(i,1-4+3,indvarsst(11))*(d1_stemp_dy_0_i1m4p3k)))**2)*2)**0.5



!***********************************************************
!                                                           
! building source terms in RHS for layer None 3 None SS ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (stemp+ReI*nut/(k**2*eta**2))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 3 None SS *****************
!                                                           
!***********************************************************


qst(i,1-4+3,indvarsst(12)) =  (qst(i,1-4+3,indvarsst(4))+&
                    param_float(1 + 5)*q(i,1-4+3,indvars(5))/(param_float(9 + 5)**2*qst(i,1-4+3,indvarsst(2))**2))



!***********************************************************
!                                                           
! building source terms in RHS for layer None 3 None tau_wall 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (visc_t*([u]_1y)*deltayI)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer None 3 None tau_wall ***********
!                                                           
!***********************************************************


qst(i,1-4+3,indvarsst(13)) =  (((1+&
                    param_float(21 + 5))/(((q(i,1-4+3,indvars(4))-&
                    0.5_wp*(q(i,1-4+3,indvars(2))*q(i,1-4+3,indvars(2))+&
                    q(i,1-4+3,indvars(3))*q(i,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,1-4+3,indvars(4))-&
                    0.5_wp*(q(i,1-4+3,indvars(2))*q(i,1-4+3,indvars(2))+&
                    q(i,1-4+3,indvars(3))*q(i,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(i,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,1-4+3,indvars(4))-&
                    0.5_wp*(q(i,1-4+3,indvars(2))*q(i,1-4+3,indvars(2))+&
                    q(i,1-4+3,indvars(3))*q(i,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,1-4+3,indvars(4))-&
                    0.5_wp*(q(i,1-4+3,indvars(2))*q(i,1-4+3,indvars(2))+&
                    q(i,1-4+3,indvars(3))*q(i,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(i,1-4+3,indvars(1)))**3/((q(i,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,1-4+3,indvars(4))-&
                    0.5_wp*(q(i,1-4+3,indvars(2))*q(i,1-4+3,indvars(2))+&
                    q(i,1-4+3,indvars(3))*q(i,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,1-4+3,indvars(4))-&
                    0.5_wp*(q(i,1-4+3,indvars(2))*q(i,1-4+3,indvars(2))+&
                    q(i,1-4+3,indvars(3))*q(i,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(i,1-4+3,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(i,1-4+3,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,1-4+3,indvars(4))-&
                    0.5_wp*(q(i,1-4+3,indvars(2))*q(i,1-4+3,indvars(2))+&
                    q(i,1-4+3,indvars(3))*q(i,1-4+3,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,1-4+3,indvars(4))-&
                    0.5_wp*(q(i,1-4+3,indvars(2))*q(i,1-4+3,indvars(2))+&
                    q(i,1-4+3,indvars(3))*q(i,1-4+3,indvars(3)))))/param_float(4 + 5)**1.5*q(i,1-4+3,indvars(1))))*param_float(1 + 5)*(d1_stemp_dy_0_i1m4p3k)*qst(i,1-4+3,indvarsst(11)))

   enddo
