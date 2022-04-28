

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
! (2.0_wp*(dabs(0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_nxp2p0p0jk = q(nx+2+0+0,j,indvars(3))

d1_stemp_dx_0_nxp2p0m1jk = q(nx+2+0-1,j,indvars(3))

d1_stemp_dx_0_nxp2p0m2jk = q(nx+2+0-2,j,indvars(3))

d1_stemp_dx_0_nxp2p0jk = 1.5_wp*d1_stemp_dx_0_nxp2p0p0jk-&
          2.0_wp*d1_stemp_dx_0_nxp2p0m1jk+&
          0.5_wp*d1_stemp_dx_0_nxp2p0m2jk

d1_stemp_dx_0_nxp2p0jk = d1_stemp_dx_0_nxp2p0jk*param_float(1)

d1_stemp_dy_0_nxp2p0jm1k = q(nx+2+0,j-1,indvars(2))

d1_stemp_dy_0_nxp2p0jp1k = q(nx+2+0,j+1,indvars(2))

d1_stemp_dy_0_nxp2p0jk = -&
          0.5_wp*d1_stemp_dy_0_nxp2p0jm1k+&
          0.5_wp*d1_stemp_dy_0_nxp2p0jp1k

d1_stemp_dy_0_nxp2p0jk = d1_stemp_dy_0_nxp2p0jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None stemp **************
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(4)) =  (2.0_wp*(dabs(0.5_wp*(qst(nx+2+0,j,indvarsst(11))*(d1_stemp_dy_0_nxp2p0jk)-&
                    qst(nx+2+0,j,indvarsst(10))*(d1_stemp_dx_0_nxp2p0jk)))))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None SS ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (stemp+fv2*ReI*nut/(k**2.0_wp*eta**2.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None SS *****************
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(12)) =  (qst(nx+2+0,j,indvarsst(4))+&
                    (1.0_wp-&
                    (q(nx+2+0,j,indvars(5))/1.0_wp*q(nx+2+0,j,indvars(1)))/(1.0_wp+&
                    (q(nx+2+0,j,indvars(5))/1.0_wp*q(nx+2+0,j,indvars(1)))*((q(nx+2+0,j,indvars(5))/1.0_wp*q(nx+2+0,j,indvars(1)))**3.0_wp/((q(nx+2+0,j,indvars(5))/1.0_wp*q(nx+2+0,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))))*param_float(1 + 5)*q(nx+2+0,j,indvars(5))/(param_float(9 + 5)**2.0_wp*qst(nx+2+0,j,indvarsst(2))**2.0_wp))



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


qst(nx+2+0,j,indvarsst(13)) =  ((1.0_wp)*(1.0_wp+&
                    ((q(nx+2+0,j,indvars(5))/1.0_wp*q(nx+2+0,j,indvars(1)))**3.0_wp/((q(nx+2+0,j,indvars(5))/1.0_wp*q(nx+2+0,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+2+0,j,indvars(5))/1.0_wp*q(nx+2+0,j,indvars(1))))*param_float(1 + 5)*(d1_stemp_dy_0_nxp2p0jk)*qst(nx+2+0,j,indvarsst(11)))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None visc_SA 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! nut*rho
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None visc_SA ************
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(14)) =  q(nx+2+0,j,indvars(5))*q(nx+2+0,j,indvars(1))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None visc_turb 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! nut*rho*fv1
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None visc_turb **********
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(15)) =  q(nx+2+0,j,indvars(5))*q(nx+2+0,j,indvars(1))*((q(nx+2+0,j,indvars(5))/1.0_wp*q(nx+2+0,j,indvars(1)))**3.0_wp/((q(nx+2+0,j,indvars(5))/1.0_wp*q(nx+2+0,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 0 None None Pressure 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (gamma_m1)*rho*(e)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None Pressure ***********
!                                                           
!***********************************************************


qst(nx+2+0,j,indvarsst(16)) =  (param_float(3 + 5))*q(nx+2+0,j,indvars(1))*((q(nx+2+0,j,indvars(4))-&
                    0.5_wp*(q(nx+2+0,j,indvars(2))*q(nx+2+0,j,indvars(2))+&
                    q(nx+2+0,j,indvars(3))*q(nx+2+0,j,indvars(3)))))

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
! (2.0_wp*(dabs(0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_nxp2m1m1jk = q(nx+2-1-1,j,indvars(3))

d1_stemp_dx_0_nxp2m1p1jk = q(nx+2-1+1,j,indvars(3))

d1_stemp_dx_0_nxp2m1jk = -&
          0.5_wp*d1_stemp_dx_0_nxp2m1m1jk+&
          0.5_wp*d1_stemp_dx_0_nxp2m1p1jk

d1_stemp_dx_0_nxp2m1jk = d1_stemp_dx_0_nxp2m1jk*param_float(1)

d1_stemp_dy_0_nxp2m1jm1k = q(nx+2-1,j-1,indvars(2))

d1_stemp_dy_0_nxp2m1jp1k = q(nx+2-1,j+1,indvars(2))

d1_stemp_dy_0_nxp2m1jk = -&
          0.5_wp*d1_stemp_dy_0_nxp2m1jm1k+&
          0.5_wp*d1_stemp_dy_0_nxp2m1jp1k

d1_stemp_dy_0_nxp2m1jk = d1_stemp_dy_0_nxp2m1jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None stemp **************
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(4)) =  (2.0_wp*(dabs(0.5_wp*(qst(nx+2-1,j,indvarsst(11))*(d1_stemp_dy_0_nxp2m1jk)-&
                    qst(nx+2-1,j,indvarsst(10))*(d1_stemp_dx_0_nxp2m1jk)))))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None SS ****
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (stemp+fv2*ReI*nut/(k**2.0_wp*eta**2.0_wp))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None SS *****************
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(12)) =  (qst(nx+2-1,j,indvarsst(4))+&
                    (1.0_wp-&
                    (q(nx+2-1,j,indvars(5))/1.0_wp*q(nx+2-1,j,indvars(1)))/(1.0_wp+&
                    (q(nx+2-1,j,indvars(5))/1.0_wp*q(nx+2-1,j,indvars(1)))*((q(nx+2-1,j,indvars(5))/1.0_wp*q(nx+2-1,j,indvars(1)))**3.0_wp/((q(nx+2-1,j,indvars(5))/1.0_wp*q(nx+2-1,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))))*param_float(1 + 5)*q(nx+2-1,j,indvars(5))/(param_float(9 + 5)**2.0_wp*qst(nx+2-1,j,indvarsst(2))**2.0_wp))



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


qst(nx+2-1,j,indvarsst(13)) =  ((1.0_wp)*(1.0_wp+&
                    ((q(nx+2-1,j,indvars(5))/1.0_wp*q(nx+2-1,j,indvars(1)))**3.0_wp/((q(nx+2-1,j,indvars(5))/1.0_wp*q(nx+2-1,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(nx+2-1,j,indvars(5))/1.0_wp*q(nx+2-1,j,indvars(1))))*param_float(1 + 5)*(d1_stemp_dy_0_nxp2m1jk)*qst(nx+2-1,j,indvarsst(11)))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None visc_SA 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! nut*rho
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None visc_SA ************
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(14)) =  q(nx+2-1,j,indvars(5))*q(nx+2-1,j,indvars(1))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None visc_turb 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! nut*rho*fv1
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None visc_turb **********
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(15)) =  q(nx+2-1,j,indvars(5))*q(nx+2-1,j,indvars(1))*((q(nx+2-1,j,indvars(5))/1.0_wp*q(nx+2-1,j,indvars(1)))**3.0_wp/((q(nx+2-1,j,indvars(5))/1.0_wp*q(nx+2-1,j,indvars(1)))**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))



!***********************************************************
!                                                           
! building source terms in RHS for layer 1 None None Pressure 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (gamma_m1)*rho*(e)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update BC terms for layer 1 None None Pressure ***********
!                                                           
!***********************************************************


qst(nx+2-1,j,indvarsst(16)) =  (param_float(3 + 5))*q(nx+2-1,j,indvars(1))*((q(nx+2-1,j,indvars(4))-&
                    0.5_wp*(q(nx+2-1,j,indvars(2))*q(nx+2-1,j,indvars(2))+&
                    q(nx+2-1,j,indvars(3))*q(nx+2-1,j,indvars(3)))))

     enddo
