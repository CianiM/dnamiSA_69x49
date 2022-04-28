

!***********************************************************
!                                                           
! Start building layers for BC : imax jmax None ************
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
! (2.0_wp*(dabs(0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_nxp2m1m1nyp2m1k = q(nx+2-1-1,ny+2-1,indvars(3))

d1_stemp_dx_0_nxp2m1p1nyp2m1k = q(nx+2-1+1,ny+2-1,indvars(3))

d1_stemp_dx_0_nxp2m1nyp2m1k = -&
          0.5_wp*d1_stemp_dx_0_nxp2m1m1nyp2m1k+&
          0.5_wp*d1_stemp_dx_0_nxp2m1p1nyp2m1k

d1_stemp_dx_0_nxp2m1nyp2m1k = d1_stemp_dx_0_nxp2m1nyp2m1k*param_float(1)

d1_stemp_dy_0_nxp2m1nyp2m1m1k = q(nx+2-1,ny+2-1-1,indvars(2))

d1_stemp_dy_0_nxp2m1nyp2m1p1k = q(nx+2-1,ny+2-1+1,indvars(2))

d1_stemp_dy_0_nxp2m1nyp2m1k = -&
          0.5_wp*d1_stemp_dy_0_nxp2m1nyp2m1m1k+&
          0.5_wp*d1_stemp_dy_0_nxp2m1nyp2m1p1k

d1_stemp_dy_0_nxp2m1nyp2m1k = d1_stemp_dy_0_nxp2m1nyp2m1k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None stemp *****************
!                                                           
!***********************************************************


qst(nx+2-1,ny+2-1,indvarsst(4)) =  (2.0_wp*(dabs(0.5_wp*(qst(nx+2-1,ny+2-1,indvarsst(11))*(d1_stemp_dy_0_nxp2m1nyp2m1k)-&
                    qst(nx+2-1,ny+2-1,indvarsst(10))*(d1_stemp_dx_0_nxp2m1nyp2m1k)))))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 1 None SS *******
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (stemp+fv2*ReI*nut/(k**2.0_wp*eta**2.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 1 None SS ********************
!                                                           
!***********************************************************


qst(nx+2-1,ny+2-1,indvarsst(12)) =  (qst(nx+2-1,ny+2-1,indvarsst(4))+&
                    (1.0_wp-&
                    (q(nx+2-1,ny+2-1,indvars(5))/1.0_wp*q(nx+2-1,ny+2-1,indvars(1)))/(1.0_wp+&
                    (q(nx+2-1,ny+2-1,indvars(5))/1.0_wp*q(nx+2-1,ny+2-1,indvars(1)))*((q(nx+2-1,ny+2-1,indvars(5))/1.0_wp*q(nx+2-1,ny+2-1,indvars(1)))**3.0_wp/((q(nx+2-1,ny+2-1,indvars(5))/1.0_wp*q(nx+2-1,ny+2-1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))))*param_float(1 + 5)*q(nx+2-1,ny+2-1,indvars(5))/(param_float(9 + 5)**2.0_wp*qst(nx+2-1,ny+2-1,indvarsst(2))**2.0_wp))



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


qst(nx+2-1,ny+2-1,indvarsst(13)) =  ((1.0_wp)*(1.0_wp+&
                    ((q(nx+2-1,ny+2-1,indvars(5))/1.0_wp*q(nx+2-1,ny+2-1,indvars(1)))**3.0_wp/((q(nx+2-1,ny+2-1,indvars(5))/1.0_wp*q(nx+2-1,ny+2-1,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+2-1,ny+2-1,indvars(5))/1.0_wp*q(nx+2-1,ny+2-1,indvars(1))))*param_float(1 + 5)*(d1_stemp_dy_0_nxp2m1nyp2m1k)*qst(nx+2-1,ny+2-1,indvarsst(11)))



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


qst(nx+2-1,ny+2-1,indvarsst(14)) =  q(nx+2-1,ny+2-1,indvars(5))*q(nx+2-1,ny+2-1,indvars(1))



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


qst(nx+2-1,ny+2-1,indvarsst(15)) =  q(nx+2-1,ny+2-1,indvars(5))*q(nx+2-1,ny+2-1,indvars(1))*((q(nx+2-1,ny+2-1,indvars(5))/1.0_wp*q(nx+2-1,ny+2-1,indvars(1)))**3.0_wp/((q(nx+2-1,ny+2-1,indvars(5))/1.0_wp*q(nx+2-1,ny+2-1,indvars(1)))**3.0_wp+&
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


qst(nx+2-1,ny+2-1,indvarsst(16)) =  (param_float(3 + 5))*q(nx+2-1,ny+2-1,indvars(1))*((q(nx+2-1,ny+2-1,indvars(4))-&
                    0.5_wp*(q(nx+2-1,ny+2-1,indvars(2))*q(nx+2-1,ny+2-1,indvars(2))+&
                    q(nx+2-1,ny+2-1,indvars(3))*q(nx+2-1,ny+2-1,indvars(3)))))

