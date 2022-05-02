

!***********************************************************
!                                                           
! Start building layers for BC : i1 jmax None **************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 1 1 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None stemp ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! ((dabs(0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_1m2p1m1nyp2m1k = q(1-2+1-1,ny+2-1,indvars(3))

d1_stemp_dx_0_1m2p1p1nyp2m1k = q(1-2+1+1,ny+2-1,indvars(3))

d1_stemp_dx_0_1m2p1nyp2m1k = -&
          0.5_wp*d1_stemp_dx_0_1m2p1m1nyp2m1k+&
          0.5_wp*d1_stemp_dx_0_1m2p1p1nyp2m1k

d1_stemp_dx_0_1m2p1nyp2m1k = d1_stemp_dx_0_1m2p1nyp2m1k*param_float(1)

d1_stemp_dy_0_1m2p1nyp2m1m1k = q(1-2+1,ny+2-1-1,indvars(2))

d1_stemp_dy_0_1m2p1nyp2m1p1k = q(1-2+1,ny+2-1+1,indvars(2))

d1_stemp_dy_0_1m2p1nyp2m1k = -&
          0.5_wp*d1_stemp_dy_0_1m2p1nyp2m1m1k+&
          0.5_wp*d1_stemp_dy_0_1m2p1nyp2m1p1k

d1_stemp_dy_0_1m2p1nyp2m1k = d1_stemp_dy_0_1m2p1nyp2m1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None stemp *****************
!                                                           
!***********************************************************


qst(1-2+1,ny+2-1,indvarsst(4)) =  ((dabs(0.5_wp*(qst(1-2+1,ny+2-1,indvarsst(11))*(d1_stemp_dy_0_1m2p1nyp2m1k)-&
                    qst(1-2+1,ny+2-1,indvarsst(10))*(d1_stemp_dx_0_1m2p1nyp2m1k)))))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None SA *******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (stemp+fv2*ReI*nut/(k**2.0_wp*eta**2.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None SA ********************
!                                                           
!***********************************************************


qst(1-2+1,ny+2-1,indvarsst(12)) =  (qst(1-2+1,ny+2-1,indvarsst(4))+&
                    (1.0_wp-&
                    (q(1-2+1,ny+2-1,indvars(5))/1.0_wp)/(1.0_wp+&
                    (q(1-2+1,ny+2-1,indvars(5))/1.0_wp)*((q(1-2+1,ny+2-1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1,ny+2-1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))))*param_float(1 + 5)*q(1-2+1,ny+2-1,indvars(5))/(param_float(9 + 5)**2.0_wp*qst(1-2+1,ny+2-1,indvarsst(2))**2.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None SS *******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (max(SA,0.3*stemp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None SS ********************
!                                                           
!***********************************************************


qst(1-2+1,ny+2-1,indvarsst(13)) =  (max(qst(1-2+1,ny+2-1,indvarsst(12)),0.3*qst(1-2+1,ny+2-1,indvarsst(4))))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None r ********
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (min(ReI*(nut/(SS*k**2.0_wp*eta**2.0_wp)),10.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None r *********************
!                                                           
!***********************************************************


qst(1-2+1,ny+2-1,indvarsst(14)) =  (min(param_float(1 + 5)*(q(1-2+1,ny+2-1,indvars(5))/(qst(1-2+1,ny+2-1,indvarsst(13))*param_float(9 + 5)**2.0_wp*qst(1-2+1,ny+2-1,indvarsst(2))**2.0_wp)),10.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None g ********
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (r+Cw2*(r**6.0_wp-r))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None g *********************
!                                                           
!***********************************************************


qst(1-2+1,ny+2-1,indvarsst(15)) =  (qst(1-2+1,ny+2-1,indvarsst(14))+&
                    param_float(11 + 5)*(qst(1-2+1,ny+2-1,indvarsst(14))**6.0_wp-&
                    qst(1-2+1,ny+2-1,indvarsst(14))))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None fw *******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (g*((1.0_wp+Cw3**6.0_wp)/(g**6.0_wp+Cw3**6.0_wp))**(1.0_wp/6.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None fw ********************
!                                                           
!***********************************************************


qst(1-2+1,ny+2-1,indvarsst(16)) =  (qst(1-2+1,ny+2-1,indvarsst(15))*((1.0_wp+&
                    param_float(12 + 5)**6.0_wp)/(qst(1-2+1,ny+2-1,indvarsst(15))**6.0_wp+&
                    param_float(12 + 5)**6.0_wp))**(1.0_wp/6.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None tau_wall *
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (visc_t*([u]_1y)*deltayI)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None tau_wall **************
!                                                           
!***********************************************************


qst(1-2+1,ny+2-1,indvarsst(17)) =  ((1.0_wp)*(1.0_wp+&
                    ((q(1-2+1,ny+2-1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1,ny+2-1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+1,ny+2-1,indvars(5))/1.0_wp))*param_float(1 + 5)*(d1_stemp_dy_0_1m2p1nyp2m1k)*qst(1-2+1,ny+2-1,indvarsst(11)))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None visc_SA **
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! nut*rho
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None visc_SA ***************
!                                                           
!***********************************************************


qst(1-2+1,ny+2-1,indvarsst(18)) =  q(1-2+1,ny+2-1,indvars(5))*q(1-2+1,ny+2-1,indvars(1))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None visc_turb 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! nut*rho*fv1
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None visc_turb *************
!                                                           
!***********************************************************


qst(1-2+1,ny+2-1,indvarsst(19)) =  q(1-2+1,ny+2-1,indvars(5))*q(1-2+1,ny+2-1,indvars(1))*((q(1-2+1,ny+2-1,indvars(5))/1.0_wp)**3.0_wp/((q(1-2+1,ny+2-1,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None Pressure *
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (gamma_m1)*rho*(e)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None Pressure **************
!                                                           
!***********************************************************


qst(1-2+1,ny+2-1,indvarsst(20)) =  (param_float(3 + 5))*q(1-2+1,ny+2-1,indvars(1))*((q(1-2+1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(1-2+1,ny+2-1,indvars(2))*q(1-2+1,ny+2-1,indvars(2))+&
                    q(1-2+1,ny+2-1,indvars(3))*q(1-2+1,ny+2-1,indvars(3)))))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None chi_coeff 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! nut/visc*rho
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None chi_coeff *************
!                                                           
!***********************************************************


qst(1-2+1,ny+2-1,indvarsst(21)) =  q(1-2+1,ny+2-1,indvars(5))/1.0_wp*q(1-2+1,ny+2-1,indvars(1))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None Production 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! Cb1*(1.0_wp-ft2)*SS*rho*nut
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None Production ************
!                                                           
!***********************************************************


qst(1-2+1,ny+2-1,indvarsst(22)) =  param_float(6 + 5)*(1.0_wp-&
                    (param_float(16 + 5)*exp(-&
                    param_float(17 + 5)*(q(1-2+1,ny+2-1,indvars(5))/1.0_wp)**2.0_wp)))*qst(1-2+1,ny+2-1,indvarsst(13))*q(1-2+1,ny+2-1,indvars(1))*q(1-2+1,ny+2-1,indvars(5))

