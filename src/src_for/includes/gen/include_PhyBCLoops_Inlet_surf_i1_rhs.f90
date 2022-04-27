

!***********************************************************
!                                                           
! Start building layers for BC : i1 None None **************
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
! building source terms in RHS for layer 0 None None d(rho)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (1.0_wp/c**2*((((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI)*(gamma_m1)+0.5_wp*(((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI+((c+u)*(c*rho*[u]_1x+1.0_wp*[p]_1x))*deltaxI)))+[rho*v]_1y
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_rhs_rho_dx_0_1m4p0p0jk = q(1-4+0+0,j,indvars(2))

d1_rhs_rho_dx_0_1m4p0p1jk = q(1-4+0+1,j,indvars(2))

d1_rhs_rho_dx_0_1m4p0p2jk = q(1-4+0+2,j,indvars(2))

d1_rhs_rho_dx_0_1m4p0jk = -&
          1.5_wp*d1_rhs_rho_dx_0_1m4p0p0jk+&
          2.0_wp*d1_rhs_rho_dx_0_1m4p0p1jk-&
          0.5_wp*d1_rhs_rho_dx_0_1m4p0p2jk

d1_rhs_rho_dx_0_1m4p0jk = d1_rhs_rho_dx_0_1m4p0jk*param_float(1)

d1_rhs_rho_dx_1_1m4p0p0jk = (param_float(3 + 5))*q(1-4+0+0,j,indvars(1))*((q(1-4+0+0,j,indvars(4))-&
                    0.5_wp*(q(1-4+0+0,j,indvars(2))*q(1-4+0+0,j,indvars(2))+&
                    q(1-4+0+0,j,indvars(3))*q(1-4+0+0,j,indvars(3)))))

d1_rhs_rho_dx_1_1m4p0p1jk = (param_float(3 + 5))*q(1-4+0+1,j,indvars(1))*((q(1-4+0+1,j,indvars(4))-&
                    0.5_wp*(q(1-4+0+1,j,indvars(2))*q(1-4+0+1,j,indvars(2))+&
                    q(1-4+0+1,j,indvars(3))*q(1-4+0+1,j,indvars(3)))))

d1_rhs_rho_dx_1_1m4p0p2jk = (param_float(3 + 5))*q(1-4+0+2,j,indvars(1))*((q(1-4+0+2,j,indvars(4))-&
                    0.5_wp*(q(1-4+0+2,j,indvars(2))*q(1-4+0+2,j,indvars(2))+&
                    q(1-4+0+2,j,indvars(3))*q(1-4+0+2,j,indvars(3)))))

d1_rhs_rho_dx_1_1m4p0jk = -&
          1.5_wp*d1_rhs_rho_dx_1_1m4p0p0jk+&
          2.0_wp*d1_rhs_rho_dx_1_1m4p0p1jk-&
          0.5_wp*d1_rhs_rho_dx_1_1m4p0p2jk

d1_rhs_rho_dx_1_1m4p0jk = d1_rhs_rho_dx_1_1m4p0jk*param_float(1)

d1_rhs_rho_dy_0_1m4p0jm1k = q(1-4+0,j-1,indvars(1))*q(1-4+0,j-1,indvars(3))

d1_rhs_rho_dy_0_1m4p0jp1k = q(1-4+0,j+1,indvars(1))*q(1-4+0,j+1,indvars(3))

d1_rhs_rho_dy_0_1m4p0jk = -&
          0.5_wp*d1_rhs_rho_dy_0_1m4p0jm1k+&
          0.5_wp*d1_rhs_rho_dy_0_1m4p0jp1k

d1_rhs_rho_dy_0_1m4p0jk = d1_rhs_rho_dy_0_1m4p0jk*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer 0 None None d(rho)/dt **********
!                                                           
!***********************************************************


rhs(1-4+0,j,indvars(1)) =   -  ( (1.0_wp/(param_float(23 + 5)*(param_float(3 + 5))*q(1-4+0,j,indvars(1))*((q(1-4+0,j,indvars(4))-&
                    0.5_wp*(q(1-4+0,j,indvars(2))*q(1-4+0,j,indvars(2))+&
                    q(1-4+0,j,indvars(3))*q(1-4+0,j,indvars(3)))))/q(1-4+0,j,indvars(1)))**0.5**2*(((((param_float(23 + 5)*(param_float(3 + 5))*q(1-4+0,j,indvars(1))*((q(1-4+0,j,indvars(4))-&
                    0.5_wp*(q(1-4+0,j,indvars(2))*q(1-4+0,j,indvars(2))+&
                    q(1-4+0,j,indvars(3))*q(1-4+0,j,indvars(3)))))/q(1-4+0,j,indvars(1)))**0.5+&
                    q(1-4+0,j,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(1-4+0,j,indvars(1))*((q(1-4+0,j,indvars(4))-&
                    0.5_wp*(q(1-4+0,j,indvars(2))*q(1-4+0,j,indvars(2))+&
                    q(1-4+0,j,indvars(3))*q(1-4+0,j,indvars(3)))))/q(1-4+0,j,indvars(1)))**0.5*q(1-4+0,j,indvars(1))*d1_rhs_rho_dx_0_1m4p0jk+&
                    1.0_wp*d1_rhs_rho_dx_1_1m4p0jk))*qst(1-4+0,j,indvarsst(10)))*(param_float(3 + 5))+&
                    0.5_wp*((((param_float(23 + 5)*(param_float(3 + 5))*q(1-4+0,j,indvars(1))*((q(1-4+0,j,indvars(4))-&
                    0.5_wp*(q(1-4+0,j,indvars(2))*q(1-4+0,j,indvars(2))+&
                    q(1-4+0,j,indvars(3))*q(1-4+0,j,indvars(3)))))/q(1-4+0,j,indvars(1)))**0.5+&
                    q(1-4+0,j,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(1-4+0,j,indvars(1))*((q(1-4+0,j,indvars(4))-&
                    0.5_wp*(q(1-4+0,j,indvars(2))*q(1-4+0,j,indvars(2))+&
                    q(1-4+0,j,indvars(3))*q(1-4+0,j,indvars(3)))))/q(1-4+0,j,indvars(1)))**0.5*q(1-4+0,j,indvars(1))*d1_rhs_rho_dx_0_1m4p0jk+&
                    1.0_wp*d1_rhs_rho_dx_1_1m4p0jk))*qst(1-4+0,j,indvarsst(10))+&
                    (((param_float(23 + 5)*(param_float(3 + 5))*q(1-4+0,j,indvars(1))*((q(1-4+0,j,indvars(4))-&
                    0.5_wp*(q(1-4+0,j,indvars(2))*q(1-4+0,j,indvars(2))+&
                    q(1-4+0,j,indvars(3))*q(1-4+0,j,indvars(3)))))/q(1-4+0,j,indvars(1)))**0.5+&
                    q(1-4+0,j,indvars(2)))*((param_float(23 + 5)*(param_float(3 + 5))*q(1-4+0,j,indvars(1))*((q(1-4+0,j,indvars(4))-&
                    0.5_wp*(q(1-4+0,j,indvars(2))*q(1-4+0,j,indvars(2))+&
                    q(1-4+0,j,indvars(3))*q(1-4+0,j,indvars(3)))))/q(1-4+0,j,indvars(1)))**0.5*q(1-4+0,j,indvars(1))*d1_rhs_rho_dx_0_1m4p0jk+&
                    1.0_wp*d1_rhs_rho_dx_1_1m4p0jk))*qst(1-4+0,j,indvarsst(10)))))+d1_rhs_rho_dy_0_1m4p0jk ) 

     enddo
