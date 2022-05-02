

!***********************************************************
!                                                           
! Start building layers for BC : i1 j1 None ****************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 0 0 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None stemp ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ((dabs(0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_1m2p0p01m2p0k = q(1-2+0+0,1-2+0,indvars(3))

d1_stemp_dx_0_1m2p0p11m2p0k = q(1-2+0+1,1-2+0,indvars(3))

d1_stemp_dx_0_1m2p0p21m2p0k = q(1-2+0+2,1-2+0,indvars(3))

d1_stemp_dx_0_1m2p01m2p0k = -&
          1.5_wp*d1_stemp_dx_0_1m2p0p01m2p0k+&
          2.0_wp*d1_stemp_dx_0_1m2p0p11m2p0k-&
          0.5_wp*d1_stemp_dx_0_1m2p0p21m2p0k

d1_stemp_dx_0_1m2p01m2p0k = d1_stemp_dx_0_1m2p01m2p0k*param_float(1)

d1_stemp_dy_0_1m2p01m2p0p0k = q(1-2+0,1-2+0+0,indvars(2))

d1_stemp_dy_0_1m2p01m2p0p1k = q(1-2+0,1-2+0+1,indvars(2))

d1_stemp_dy_0_1m2p01m2p0p2k = q(1-2+0,1-2+0+2,indvars(2))

d1_stemp_dy_0_1m2p01m2p0k = -&
          1.5_wp*d1_stemp_dy_0_1m2p01m2p0p0k+&
          2.0_wp*d1_stemp_dy_0_1m2p01m2p0p1k-&
          0.5_wp*d1_stemp_dy_0_1m2p01m2p0p2k

d1_stemp_dy_0_1m2p01m2p0k = d1_stemp_dy_0_1m2p01m2p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None stemp *****************
!                                                           
!***********************************************************


qst(1-2+0,1-2+0,indvarsst(4)) =  ((dabs(0.5_wp*(qst(1-2+0,1-2+0,indvarsst(11))*(d1_stemp_dy_0_1m2p01m2p0k)-&
                    qst(1-2+0,1-2+0,indvarsst(10))*(d1_stemp_dx_0_1m2p01m2p0k)))))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None SA *******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (stemp+fv2*ReI*nut/(k**2.0_wp*eta**2.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None SA ********************
!                                                           
!***********************************************************


qst(1-2+0,1-2+0,indvarsst(12)) =  (qst(1-2+0,1-2+0,indvarsst(4))+&
                    (1.0_wp-&
                    (q(1-2+0,1-2+0,indvars(5))/1.0_wp)/(1.0_wp+&
                    (q(1-2+0,1-2+0,indvars(5))/1.0_wp)*((q(1-2+0,1-2+0,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+0,1-2+0,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))))*param_float(1 + 5)*q(1-2+0,1-2+0,indvars(5))/(param_float(9 + 5)**2.0_wp*qst(1-2+0,1-2+0,indvarsst(2))**2.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None SS *******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (max(SA,0.3*stemp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None SS ********************
!                                                           
!***********************************************************


qst(1-2+0,1-2+0,indvarsst(13)) =  (max(qst(1-2+0,1-2+0,indvarsst(12)),0.3*qst(1-2+0,1-2+0,indvarsst(4))))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None r ********
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (min(ReI*(nut/(SS*k**2.0_wp*eta**2.0_wp)),10.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None r *********************
!                                                           
!***********************************************************


qst(1-2+0,1-2+0,indvarsst(14)) =  (min(param_float(1 + 5)*(q(1-2+0,1-2+0,indvars(5))/(qst(1-2+0,1-2+0,indvarsst(13))*param_float(9 + 5)**2.0_wp*qst(1-2+0,1-2+0,indvarsst(2))**2.0_wp)),10.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None g ********
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (r+Cw2*(r**6.0_wp-r))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None g *********************
!                                                           
!***********************************************************


qst(1-2+0,1-2+0,indvarsst(15)) =  (qst(1-2+0,1-2+0,indvarsst(14))+&
                    param_float(11 + 5)*(qst(1-2+0,1-2+0,indvarsst(14))**6.0_wp-&
                    qst(1-2+0,1-2+0,indvarsst(14))))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None fw *******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (g*((1.0_wp+Cw3**6.0_wp)/(g**6.0_wp+Cw3**6.0_wp))**(1.0_wp/6.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None fw ********************
!                                                           
!***********************************************************


qst(1-2+0,1-2+0,indvarsst(16)) =  (qst(1-2+0,1-2+0,indvarsst(15))*((1.0_wp+&
                    param_float(12 + 5)**6.0_wp)/(qst(1-2+0,1-2+0,indvarsst(15))**6.0_wp+&
                    param_float(12 + 5)**6.0_wp))**(1.0_wp/6.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None tau_wall *
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (visc_t*([u]_1y)*deltayI)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None tau_wall **************
!                                                           
!***********************************************************


qst(1-2+0,1-2+0,indvarsst(17)) =  ((1.0_wp)*(1.0_wp+&
                    ((q(1-2+0,1-2+0,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+0,1-2+0,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+0,1-2+0,indvars(5))/1.0_wp))*param_float(1 + 5)*(d1_stemp_dy_0_1m2p01m2p0k)*qst(1-2+0,1-2+0,indvarsst(11)))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None visc_SA **
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! nut*rho
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None visc_SA ***************
!                                                           
!***********************************************************


qst(1-2+0,1-2+0,indvarsst(18)) =  q(1-2+0,1-2+0,indvars(5))*q(1-2+0,1-2+0,indvars(1))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None visc_turb 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! nut*rho*fv1
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None visc_turb *************
!                                                           
!***********************************************************


qst(1-2+0,1-2+0,indvarsst(19)) =  q(1-2+0,1-2+0,indvars(5))*q(1-2+0,1-2+0,indvars(1))*((q(1-2+0,1-2+0,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+0,1-2+0,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None Pressure *
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (gamma_m1)*rho*(e)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None Pressure **************
!                                                           
!***********************************************************


qst(1-2+0,1-2+0,indvarsst(20)) =  (param_float(3 + 5))*q(1-2+0,1-2+0,indvars(1))*((q(1-2+0,1-2+0,indvars(4))-&
                    0.5_wp*(q(1-2+0,1-2+0,indvars(2))*q(1-2+0,1-2+0,indvars(2))+&
                    q(1-2+0,1-2+0,indvars(3))*q(1-2+0,1-2+0,indvars(3)))))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None chi_coeff 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! nut/visc*rho
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None chi_coeff *************
!                                                           
!***********************************************************


qst(1-2+0,1-2+0,indvarsst(21)) =  q(1-2+0,1-2+0,indvars(5))/1.0_wp*q(1-2+0,1-2+0,indvars(1))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 0 None Production 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! Cb1*(1.0_wp-ft2)*SS*rho*nut
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 0 None Production ************
!                                                           
!***********************************************************


qst(1-2+0,1-2+0,indvarsst(22)) =  param_float(6 + 5)*(1.0_wp-&
                    (param_float(16 + 5)*exp(-&
                    param_float(17 + 5)*(q(1-2+0,1-2+0,indvars(5))/1.0_wp)**2.0_wp)))*qst(1-2+0,1-2+0,indvarsst(13))*q(1-2+0,1-2+0,indvars(1))*q(1-2+0,1-2+0,indvars(5))

