

!***********************************************************
!                                                           
! Start building layers for BC : i1 j1 None ****************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 0 1 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None stemp ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (2.0_wp*(dabs(0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_1m2p0p01m2p1k = q(1-2+0+0,1-2+1,indvars(3))

d1_stemp_dx_0_1m2p0p11m2p1k = q(1-2+0+1,1-2+1,indvars(3))

d1_stemp_dx_0_1m2p0p21m2p1k = q(1-2+0+2,1-2+1,indvars(3))

d1_stemp_dx_0_1m2p01m2p1k = -&
          1.5_wp*d1_stemp_dx_0_1m2p0p01m2p1k+&
          2.0_wp*d1_stemp_dx_0_1m2p0p11m2p1k-&
          0.5_wp*d1_stemp_dx_0_1m2p0p21m2p1k

d1_stemp_dx_0_1m2p01m2p1k = d1_stemp_dx_0_1m2p01m2p1k*param_float(1)

d1_stemp_dy_0_1m2p01m2p1m1k = q(1-2+0,1-2+1-1,indvars(2))

d1_stemp_dy_0_1m2p01m2p1p1k = q(1-2+0,1-2+1+1,indvars(2))

d1_stemp_dy_0_1m2p01m2p1k = -&
          0.5_wp*d1_stemp_dy_0_1m2p01m2p1m1k+&
          0.5_wp*d1_stemp_dy_0_1m2p01m2p1p1k

d1_stemp_dy_0_1m2p01m2p1k = d1_stemp_dy_0_1m2p01m2p1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None stemp *****************
!                                                           
!***********************************************************


qst(1-2+0,1-2+1,indvarsst(4)) =  (2.0_wp*(dabs(0.5_wp*(qst(1-2+0,1-2+1,indvarsst(11))*(d1_stemp_dy_0_1m2p01m2p1k)-&
                    qst(1-2+0,1-2+1,indvarsst(10))*(d1_stemp_dx_0_1m2p01m2p1k)))))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None SS *******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (stemp+fv2*ReI*nut/(k**2.0_wp*eta**2.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None SS ********************
!                                                           
!***********************************************************


qst(1-2+0,1-2+1,indvarsst(12)) =  (qst(1-2+0,1-2+1,indvarsst(4))+&
                    (1.0_wp-&
                    (q(1-2+0,1-2+1,indvars(5))/1.0_wp*q(1-2+0,1-2+1,indvars(1)))/(1.0_wp+&
                    (q(1-2+0,1-2+1,indvars(5))/1.0_wp*q(1-2+0,1-2+1,indvars(1)))*((q(1-2+0,1-2+1,indvars(5))/1.0_wp*q(1-2+0,1-2+1,indvars(1)))**3.0_wp/((q(1-2+0,1-2+1,indvars(5))/1.0_wp*q(1-2+0,1-2+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))))*param_float(1 + 5)*q(1-2+0,1-2+1,indvars(5))/(param_float(9 + 5)**2.0_wp*qst(1-2+0,1-2+1,indvarsst(2))**2.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None tau_wall *
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (visc_t*([u]_1y)*deltayI)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None tau_wall **************
!                                                           
!***********************************************************


qst(1-2+0,1-2+1,indvarsst(13)) =  ((1.0_wp)*(1.0_wp+&
                    ((q(1-2+0,1-2+1,indvars(5))/1.0_wp*q(1-2+0,1-2+1,indvars(1)))**3.0_wp/((q(1-2+0,1-2+1,indvars(5))/1.0_wp*q(1-2+0,1-2+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(1-2+0,1-2+1,indvars(5))/1.0_wp*q(1-2+0,1-2+1,indvars(1))))*param_float(1 + 5)*(d1_stemp_dy_0_1m2p01m2p1k)*qst(1-2+0,1-2+1,indvarsst(11)))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None visc_SA **
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! nut*rho
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None visc_SA ***************
!                                                           
!***********************************************************


qst(1-2+0,1-2+1,indvarsst(14)) =  q(1-2+0,1-2+1,indvars(5))*q(1-2+0,1-2+1,indvars(1))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None visc_turb 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! nut*rho*fv1
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None visc_turb *************
!                                                           
!***********************************************************


qst(1-2+0,1-2+1,indvarsst(15)) =  q(1-2+0,1-2+1,indvars(5))*q(1-2+0,1-2+1,indvars(1))*((q(1-2+0,1-2+1,indvars(5))/1.0_wp*q(1-2+0,1-2+1,indvars(1)))**3.0_wp/((q(1-2+0,1-2+1,indvars(5))/1.0_wp*q(1-2+0,1-2+1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 1 None Pressure *
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (gamma_m1)*rho*(e)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 1 None Pressure **************
!                                                           
!***********************************************************


qst(1-2+0,1-2+1,indvarsst(16)) =  (param_float(3 + 5))*q(1-2+0,1-2+1,indvars(1))*((q(1-2+0,1-2+1,indvars(4))-&
                    0.5_wp*(q(1-2+0,1-2+1,indvars(2))*q(1-2+0,1-2+1,indvars(2))+&
                    q(1-2+0,1-2+1,indvars(3))*q(1-2+0,1-2+1,indvars(3)))))

