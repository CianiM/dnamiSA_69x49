

!***********************************************************
!                                                           
! Start building layers for BC : imax j1 None **************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! BC layer: 2 2 None ***************************************
!                                                           
!***********************************************************




!***********************************************************
!                                                           
! building source terms in RHS for layer 2 2 None stemp ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (2.0_wp*(dabs(0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_nxp4m2m11m4p2k = q(nx+4-2-1,1-4+2,indvars(3))

d1_stemp_dx_0_nxp4m2p11m4p2k = q(nx+4-2+1,1-4+2,indvars(3))

d1_stemp_dx_0_nxp4m21m4p2k = -&
          0.5_wp*d1_stemp_dx_0_nxp4m2m11m4p2k+&
          0.5_wp*d1_stemp_dx_0_nxp4m2p11m4p2k

d1_stemp_dx_0_nxp4m21m4p2k = d1_stemp_dx_0_nxp4m21m4p2k*param_float(1)

d1_stemp_dy_0_nxp4m21m4p2m1k = q(nx+4-2,1-4+2-1,indvars(2))

d1_stemp_dy_0_nxp4m21m4p2p1k = q(nx+4-2,1-4+2+1,indvars(2))

d1_stemp_dy_0_nxp4m21m4p2k = -&
          0.5_wp*d1_stemp_dy_0_nxp4m21m4p2m1k+&
          0.5_wp*d1_stemp_dy_0_nxp4m21m4p2p1k

d1_stemp_dy_0_nxp4m21m4p2k = d1_stemp_dy_0_nxp4m21m4p2k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 2 2 None stemp *****************
!                                                           
!***********************************************************


qst(nx+4-2,1-4+2,indvarsst(4)) =  (2.0_wp*(dabs(0.5_wp*(qst(nx+4-2,1-4+2,indvarsst(11))*(d1_stemp_dy_0_nxp4m21m4p2k)-&
                    qst(nx+4-2,1-4+2,indvarsst(10))*(d1_stemp_dx_0_nxp4m21m4p2k)))))



!***********************************************************
!                                                           
! building source terms in RHS for layer 2 2 None SS *******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (stemp+fv2*ReI*nut/(k**2.0_wp*eta**2.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 2 2 None SS ********************
!                                                           
!***********************************************************


qst(nx+4-2,1-4+2,indvarsst(12)) =  (qst(nx+4-2,1-4+2,indvarsst(4))+&
                    (1.0_wp-&
                    (q(nx+4-2,1-4+2,indvars(5))/1.0_wp*q(nx+4-2,1-4+2,indvars(1)))/(1.0_wp+&
                    (q(nx+4-2,1-4+2,indvars(5))/1.0_wp*q(nx+4-2,1-4+2,indvars(1)))*((q(nx+4-2,1-4+2,indvars(5))/1.0_wp*q(nx+4-2,1-4+2,indvars(1)))**3.0_wp/((q(nx+4-2,1-4+2,indvars(5))/1.0_wp*q(nx+4-2,1-4+2,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))))*param_float(1 + 5)*q(nx+4-2,1-4+2,indvars(5))/(param_float(9 + 5)**2.0_wp*qst(nx+4-2,1-4+2,indvarsst(2))**2.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 2 2 None tau_wall *
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (visc_t*([u]_1y)*deltayI)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 2 2 None tau_wall **************
!                                                           
!***********************************************************


qst(nx+4-2,1-4+2,indvarsst(13)) =  ((1.0_wp)*(1.0_wp+&
                    ((q(nx+4-2,1-4+2,indvars(5))/1.0_wp*q(nx+4-2,1-4+2,indvars(1)))**3.0_wp/((q(nx+4-2,1-4+2,indvars(5))/1.0_wp*q(nx+4-2,1-4+2,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+4-2,1-4+2,indvars(5))/1.0_wp*q(nx+4-2,1-4+2,indvars(1))))*param_float(1 + 5)*(d1_stemp_dy_0_nxp4m21m4p2k)*qst(nx+4-2,1-4+2,indvarsst(11)))



!***********************************************************
!                                                           
! building source terms in RHS for layer 2 2 None visc_SA **
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! nut*rho
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 2 2 None visc_SA ***************
!                                                           
!***********************************************************


qst(nx+4-2,1-4+2,indvarsst(14)) =  q(nx+4-2,1-4+2,indvars(5))*q(nx+4-2,1-4+2,indvars(1))



!***********************************************************
!                                                           
! building source terms in RHS for layer 2 2 None visc_turb 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! nut*rho*fv1
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 2 2 None visc_turb *************
!                                                           
!***********************************************************


qst(nx+4-2,1-4+2,indvarsst(15)) =  q(nx+4-2,1-4+2,indvars(5))*q(nx+4-2,1-4+2,indvars(1))*((q(nx+4-2,1-4+2,indvars(5))/1.0_wp*q(nx+4-2,1-4+2,indvars(1)))**3.0_wp/((q(nx+4-2,1-4+2,indvars(5))/1.0_wp*q(nx+4-2,1-4+2,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 2 2 None Pressure *
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (gamma_m1)*rho*(e)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 2 2 None Pressure **************
!                                                           
!***********************************************************


qst(nx+4-2,1-4+2,indvarsst(16)) =  (param_float(3 + 5))*q(nx+4-2,1-4+2,indvars(1))*((q(nx+4-2,1-4+2,indvars(4))-&
                    0.5_wp*(q(nx+4-2,1-4+2,indvars(2))*q(nx+4-2,1-4+2,indvars(2))+&
                    q(nx+4-2,1-4+2,indvars(3))*q(nx+4-2,1-4+2,indvars(3)))))

