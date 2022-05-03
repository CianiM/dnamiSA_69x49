

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
! building source terms in RHS for layer None 0 None d(rho)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (deltaxI*([rho*u]_1x))*symm+(rho*([v]_1y)*deltayI)*wall
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_rhs_rho_dx_0_im11m2p0k = q(i-1,1-2+0,indvars(1))*q(i-1,1-2+0,indvars(2))

d1_rhs_rho_dx_0_ip11m2p0k = q(i+1,1-2+0,indvars(1))*q(i+1,1-2+0,indvars(2))

d1_rhs_rho_dx_0_i1m2p0k = -&
          0.5_wp*d1_rhs_rho_dx_0_im11m2p0k+&
          0.5_wp*d1_rhs_rho_dx_0_ip11m2p0k

d1_rhs_rho_dx_0_i1m2p0k = d1_rhs_rho_dx_0_i1m2p0k*param_float(1)

d1_rhs_rho_dy_0_i1m2p0p0k = q(i,1-2+0+0,indvars(3))

d1_rhs_rho_dy_0_i1m2p0p1k = q(i,1-2+0+1,indvars(3))

d1_rhs_rho_dy_0_i1m2p0p2k = q(i,1-2+0+2,indvars(3))

d1_rhs_rho_dy_0_i1m2p0k = -&
          1.5_wp*d1_rhs_rho_dy_0_i1m2p0p0k+&
          2.0_wp*d1_rhs_rho_dy_0_i1m2p0p1k-&
          0.5_wp*d1_rhs_rho_dy_0_i1m2p0p2k

d1_rhs_rho_dy_0_i1m2p0k = d1_rhs_rho_dy_0_i1m2p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None d(rho)/dt **********
!                                                           
!***********************************************************


rhs(i,1-2+0,indvars(1)) =   -  ( (qst(i,1-2+0,indvarsst(10))*(d1_rhs_rho_dx_0_i1m2p0k))*qst(i,1-2+0,indvarsst(5))+&
                    (q(i,1-2+0,indvars(1))*(d1_rhs_rho_dy_0_i1m2p0k)*qst(i,1-2+0,indvarsst(11)))*dabs(1.0_wp-&
                    qst(i,1-2+0,indvarsst(5))) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None d(rho u)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (deltaxI*([-4.0_wp/3.0_wp*visc_t*({u}_1x)*deltaxI]_1x)+deltaxI*([rho*u*u+p]_1x))*symm
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_rhs_u_dxdx_0_0_im11m2p0k_im1m11m2p0k = q(i-1-1,1-2+0,indvars(2))

d2_rhs_u_dxdx_0_0_im11m2p0k_im1p11m2p0k = q(i-1+1,1-2+0,indvars(2))

d2_rhs_u_dxdx_0_0_im11m2p0k = -&
          0.5_wp*d2_rhs_u_dxdx_0_0_im11m2p0k_im1m11m2p0k+&
          0.5_wp*d2_rhs_u_dxdx_0_0_im11m2p0k_im1p11m2p0k

d2_rhs_u_dxdx_0_0_im11m2p0k = d2_rhs_u_dxdx_0_0_im11m2p0k*param_float(1)

d2_rhs_u_dxdx_0_0_ip11m2p0k_ip1m11m2p0k = q(i+1-1,1-2+0,indvars(2))

d2_rhs_u_dxdx_0_0_ip11m2p0k_ip1p11m2p0k = q(i+1+1,1-2+0,indvars(2))

d2_rhs_u_dxdx_0_0_ip11m2p0k = -&
          0.5_wp*d2_rhs_u_dxdx_0_0_ip11m2p0k_ip1m11m2p0k+&
          0.5_wp*d2_rhs_u_dxdx_0_0_ip11m2p0k_ip1p11m2p0k

d2_rhs_u_dxdx_0_0_ip11m2p0k = d2_rhs_u_dxdx_0_0_ip11m2p0k*param_float(1)

d1_rhs_u_dx_1_im11m2p0k = q(i-1,1-2+0,indvars(1))*q(i-1,1-2+0,indvars(2))*q(i-1,1-2+0,indvars(2))+(param_float(3 + 5))*q(i-1,1-2+0,indvars(1))*((q(i-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(i-1,1-2+0,indvars(2))*q(i-1,1-2+0,indvars(2))+&
                    q(i-1,1-2+0,indvars(3))*q(i-1,1-2+0,indvars(3)))))

d1_rhs_u_dx_1_ip11m2p0k = q(i+1,1-2+0,indvars(1))*q(i+1,1-2+0,indvars(2))*q(i+1,1-2+0,indvars(2))+(param_float(3 + 5))*q(i+1,1-2+0,indvars(1))*((q(i+1,1-2+0,indvars(4))-&
                    0.5_wp*(q(i+1,1-2+0,indvars(2))*q(i+1,1-2+0,indvars(2))+&
                    q(i+1,1-2+0,indvars(3))*q(i+1,1-2+0,indvars(3)))))

d1_rhs_u_dx_1_i1m2p0k = -&
          0.5_wp*d1_rhs_u_dx_1_im11m2p0k+&
          0.5_wp*d1_rhs_u_dx_1_ip11m2p0k

d1_rhs_u_dx_1_i1m2p0k = d1_rhs_u_dx_1_i1m2p0k*param_float(1)

d1_rhs_u_dx_0_im11m2p0k = -4.0_wp/3.0_wp*(1.0_wp)*(1.0_wp+&
                    ((q(i-1,1-2+0,indvars(5))/1.0_wp)**3.0_wp/((q(i-1,1-2+0,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i-1,1-2+0,indvars(5))/1.0_wp))*param_float(1 + 5)*(d2_rhs_u_dxdx_0_0_im11m2p0k)*qst(i-1,1-2+0,indvarsst(10))

d1_rhs_u_dx_0_ip11m2p0k = -4.0_wp/3.0_wp*(1.0_wp)*(1.0_wp+&
                    ((q(i+1,1-2+0,indvars(5))/1.0_wp)**3.0_wp/((q(i+1,1-2+0,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i+1,1-2+0,indvars(5))/1.0_wp))*param_float(1 + 5)*(d2_rhs_u_dxdx_0_0_ip11m2p0k)*qst(i+1,1-2+0,indvarsst(10))

d1_rhs_u_dx_0_i1m2p0k = -&
          0.5_wp*d1_rhs_u_dx_0_im11m2p0k+&
          0.5_wp*d1_rhs_u_dx_0_ip11m2p0k

d1_rhs_u_dx_0_i1m2p0k = d1_rhs_u_dx_0_i1m2p0k*param_float(1)



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None d(rho u)/dt ********
!                                                           
!***********************************************************


rhs(i,1-2+0,indvars(2)) =   -  ( (qst(i,1-2+0,indvarsst(10))*(d1_rhs_u_dx_0_i1m2p0k)+&
                    qst(i,1-2+0,indvarsst(10))*(d1_rhs_u_dx_1_i1m2p0k))*qst(i,1-2+0,indvarsst(5)) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None d(rho v)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (deltaxI*([-visc_t*({v}_1x)*deltaxI]_1x)+deltaxI*([rho*u*v]_1x))*symm
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_rhs_v_dxdx_0_0_im11m2p0k_im1m11m2p0k = q(i-1-1,1-2+0,indvars(3))

d2_rhs_v_dxdx_0_0_im11m2p0k_im1p11m2p0k = q(i-1+1,1-2+0,indvars(3))

d2_rhs_v_dxdx_0_0_im11m2p0k = -&
          0.5_wp*d2_rhs_v_dxdx_0_0_im11m2p0k_im1m11m2p0k+&
          0.5_wp*d2_rhs_v_dxdx_0_0_im11m2p0k_im1p11m2p0k

d2_rhs_v_dxdx_0_0_im11m2p0k = d2_rhs_v_dxdx_0_0_im11m2p0k*param_float(1)

d2_rhs_v_dxdx_0_0_ip11m2p0k_ip1m11m2p0k = q(i+1-1,1-2+0,indvars(3))

d2_rhs_v_dxdx_0_0_ip11m2p0k_ip1p11m2p0k = q(i+1+1,1-2+0,indvars(3))

d2_rhs_v_dxdx_0_0_ip11m2p0k = -&
          0.5_wp*d2_rhs_v_dxdx_0_0_ip11m2p0k_ip1m11m2p0k+&
          0.5_wp*d2_rhs_v_dxdx_0_0_ip11m2p0k_ip1p11m2p0k

d2_rhs_v_dxdx_0_0_ip11m2p0k = d2_rhs_v_dxdx_0_0_ip11m2p0k*param_float(1)

d1_rhs_v_dx_1_im11m2p0k = q(i-1,1-2+0,indvars(1))*q(i-1,1-2+0,indvars(2))*q(i-1,1-2+0,indvars(3))

d1_rhs_v_dx_1_ip11m2p0k = q(i+1,1-2+0,indvars(1))*q(i+1,1-2+0,indvars(2))*q(i+1,1-2+0,indvars(3))

d1_rhs_v_dx_1_i1m2p0k = -&
          0.5_wp*d1_rhs_v_dx_1_im11m2p0k+&
          0.5_wp*d1_rhs_v_dx_1_ip11m2p0k

d1_rhs_v_dx_1_i1m2p0k = d1_rhs_v_dx_1_i1m2p0k*param_float(1)

d1_rhs_v_dx_0_im11m2p0k = -(1.0_wp)*(1.0_wp+&
                    ((q(i-1,1-2+0,indvars(5))/1.0_wp)**3.0_wp/((q(i-1,1-2+0,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i-1,1-2+0,indvars(5))/1.0_wp))*param_float(1 + 5)*(d2_rhs_v_dxdx_0_0_im11m2p0k)*qst(i-1,1-2+0,indvarsst(10))

d1_rhs_v_dx_0_ip11m2p0k = -(1.0_wp)*(1.0_wp+&
                    ((q(i+1,1-2+0,indvars(5))/1.0_wp)**3.0_wp/((q(i+1,1-2+0,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i+1,1-2+0,indvars(5))/1.0_wp))*param_float(1 + 5)*(d2_rhs_v_dxdx_0_0_ip11m2p0k)*qst(i+1,1-2+0,indvarsst(10))

d1_rhs_v_dx_0_i1m2p0k = -&
          0.5_wp*d1_rhs_v_dx_0_im11m2p0k+&
          0.5_wp*d1_rhs_v_dx_0_ip11m2p0k

d1_rhs_v_dx_0_i1m2p0k = d1_rhs_v_dx_0_i1m2p0k*param_float(1)



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None d(rho v)/dt ********
!                                                           
!***********************************************************


rhs(i,1-2+0,indvars(3)) =   -  ( (qst(i,1-2+0,indvarsst(10))*(d1_rhs_v_dx_0_i1m2p0k)+&
                    qst(i,1-2+0,indvarsst(10))*(d1_rhs_v_dx_1_i1m2p0k))*qst(i,1-2+0,indvarsst(5)) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None d(rho et)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (deltaxI*([-4.0_wp/3.0_wp*visc_t*({u}_1x)*deltaxI*u-visc_t*({v}_1x)*deltaxI*v-kappa*({T}_1x)*deltaxI]_1x)+deltaxI*([(rho*et+p)*u]_1x))*symm+(-visc_t*([u]_1y)*([u]_1y)*deltayI**2-4.0_wp/3.0_wp*visc_t*([v]_1y)*([v]_1y)*deltayI**2-kappa*([({T}_1y)*deltayI]_1y)-kappa*([({T}_1x)*deltaxI]_1x)+(rho*et+p)*([v]_1y)*deltayI)*wall
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_rhs_et_dxdx_0_2_im11m2p0k_im1m11m2p0k = ((q(i-1-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(i-1-1,1-2+0,indvars(2))*q(i-1-1,1-2+0,indvars(2))+&
                    q(i-1-1,1-2+0,indvars(3))*q(i-1-1,1-2+0,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dxdx_0_2_im11m2p0k_im1p11m2p0k = ((q(i-1+1,1-2+0,indvars(4))-&
                    0.5_wp*(q(i-1+1,1-2+0,indvars(2))*q(i-1+1,1-2+0,indvars(2))+&
                    q(i-1+1,1-2+0,indvars(3))*q(i-1+1,1-2+0,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dxdx_0_2_im11m2p0k = -&
          0.5_wp*d2_rhs_et_dxdx_0_2_im11m2p0k_im1m11m2p0k+&
          0.5_wp*d2_rhs_et_dxdx_0_2_im11m2p0k_im1p11m2p0k

d2_rhs_et_dxdx_0_2_im11m2p0k = d2_rhs_et_dxdx_0_2_im11m2p0k*param_float(1)

d2_rhs_et_dxdx_0_2_ip11m2p0k_ip1m11m2p0k = ((q(i+1-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(i+1-1,1-2+0,indvars(2))*q(i+1-1,1-2+0,indvars(2))+&
                    q(i+1-1,1-2+0,indvars(3))*q(i+1-1,1-2+0,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dxdx_0_2_ip11m2p0k_ip1p11m2p0k = ((q(i+1+1,1-2+0,indvars(4))-&
                    0.5_wp*(q(i+1+1,1-2+0,indvars(2))*q(i+1+1,1-2+0,indvars(2))+&
                    q(i+1+1,1-2+0,indvars(3))*q(i+1+1,1-2+0,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dxdx_0_2_ip11m2p0k = -&
          0.5_wp*d2_rhs_et_dxdx_0_2_ip11m2p0k_ip1m11m2p0k+&
          0.5_wp*d2_rhs_et_dxdx_0_2_ip11m2p0k_ip1p11m2p0k

d2_rhs_et_dxdx_0_2_ip11m2p0k = d2_rhs_et_dxdx_0_2_ip11m2p0k*param_float(1)

d1_rhs_et_dx_1_im11m2p0k = (q(i-1,1-2+0,indvars(1))*q(i-1,1-2+0,indvars(4))+&
                    (param_float(3 + 5))*q(i-1,1-2+0,indvars(1))*((q(i-1,1-2+0,indvars(4))-&
                    0.5_wp*(q(i-1,1-2+0,indvars(2))*q(i-1,1-2+0,indvars(2))+&
                    q(i-1,1-2+0,indvars(3))*q(i-1,1-2+0,indvars(3))))))*q(i-1,1-2+0,indvars(2))

d1_rhs_et_dx_1_ip11m2p0k = (q(i+1,1-2+0,indvars(1))*q(i+1,1-2+0,indvars(4))+&
                    (param_float(3 + 5))*q(i+1,1-2+0,indvars(1))*((q(i+1,1-2+0,indvars(4))-&
                    0.5_wp*(q(i+1,1-2+0,indvars(2))*q(i+1,1-2+0,indvars(2))+&
                    q(i+1,1-2+0,indvars(3))*q(i+1,1-2+0,indvars(3))))))*q(i+1,1-2+0,indvars(2))

d1_rhs_et_dx_1_i1m2p0k = -&
          0.5_wp*d1_rhs_et_dx_1_im11m2p0k+&
          0.5_wp*d1_rhs_et_dx_1_ip11m2p0k

d1_rhs_et_dx_1_i1m2p0k = d1_rhs_et_dx_1_i1m2p0k*param_float(1)

d1_rhs_et_dx_0_im11m2p0k = -4.0_wp/3.0_wp*(1.0_wp)*(1.0_wp+&
                    ((q(i-1,1-2+0,indvars(5))/1.0_wp)**3.0_wp/((q(i-1,1-2+0,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i-1,1-2+0,indvars(5))/1.0_wp))*param_float(1 + 5)*(d2_rhs_u_dxdx_0_0_im11m2p0k)*qst(i-1,1-2+0,indvarsst(10))*q(i-1,1-2+0,indvars(2))-&
                    (1.0_wp)*(1.0_wp+&
                    ((q(i-1,1-2+0,indvars(5))/1.0_wp)**3.0_wp/((q(i-1,1-2+0,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i-1,1-2+0,indvars(5))/1.0_wp))*param_float(1 + 5)*(d2_rhs_v_dxdx_0_0_im11m2p0k)*qst(i-1,1-2+0,indvarsst(10))*q(i-1,1-2+0,indvars(3))-&
                    param_float(2 + 5)*(d2_rhs_et_dxdx_0_2_im11m2p0k)*qst(i-1,1-2+0,indvarsst(10))

d1_rhs_et_dx_0_ip11m2p0k = -4.0_wp/3.0_wp*(1.0_wp)*(1.0_wp+&
                    ((q(i+1,1-2+0,indvars(5))/1.0_wp)**3.0_wp/((q(i+1,1-2+0,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i+1,1-2+0,indvars(5))/1.0_wp))*param_float(1 + 5)*(d2_rhs_u_dxdx_0_0_ip11m2p0k)*qst(i+1,1-2+0,indvarsst(10))*q(i+1,1-2+0,indvars(2))-&
                    (1.0_wp)*(1.0_wp+&
                    ((q(i+1,1-2+0,indvars(5))/1.0_wp)**3.0_wp/((q(i+1,1-2+0,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i+1,1-2+0,indvars(5))/1.0_wp))*param_float(1 + 5)*(d2_rhs_v_dxdx_0_0_ip11m2p0k)*qst(i+1,1-2+0,indvarsst(10))*q(i+1,1-2+0,indvars(3))-&
                    param_float(2 + 5)*(d2_rhs_et_dxdx_0_2_ip11m2p0k)*qst(i+1,1-2+0,indvarsst(10))

d1_rhs_et_dx_0_i1m2p0k = -&
          0.5_wp*d1_rhs_et_dx_0_im11m2p0k+&
          0.5_wp*d1_rhs_et_dx_0_ip11m2p0k

d1_rhs_et_dx_0_i1m2p0k = d1_rhs_et_dx_0_i1m2p0k*param_float(1)

d1_rhs_et_dx_2_im11m2p0k = (d2_rhs_et_dxdx_0_2_im11m2p0k)*qst(i-1,1-2+0,indvarsst(10))

d1_rhs_et_dx_2_ip11m2p0k = (d2_rhs_et_dxdx_0_2_ip11m2p0k)*qst(i+1,1-2+0,indvarsst(10))

d1_rhs_et_dx_2_i1m2p0k = -&
          0.5_wp*d1_rhs_et_dx_2_im11m2p0k+&
          0.5_wp*d1_rhs_et_dx_2_ip11m2p0k

d1_rhs_et_dx_2_i1m2p0k = d1_rhs_et_dx_2_i1m2p0k*param_float(1)

d1_rhs_et_dy_0_i1m2p0p0k = q(i,1-2+0+0,indvars(2))

d1_rhs_et_dy_0_i1m2p0p1k = q(i,1-2+0+1,indvars(2))

d1_rhs_et_dy_0_i1m2p0p2k = q(i,1-2+0+2,indvars(2))

d1_rhs_et_dy_0_i1m2p0k = -&
          1.5_wp*d1_rhs_et_dy_0_i1m2p0p0k+&
          2.0_wp*d1_rhs_et_dy_0_i1m2p0p1k-&
          0.5_wp*d1_rhs_et_dy_0_i1m2p0p2k

d1_rhs_et_dy_0_i1m2p0k = d1_rhs_et_dy_0_i1m2p0k*param_float(2)

d2_rhs_et_dydy_4_0_i1m2p0p0k_i1m2p0p0p0k = ((q(i,1-2+0+0+0,indvars(4))-&
                    0.5_wp*(q(i,1-2+0+0+0,indvars(2))*q(i,1-2+0+0+0,indvars(2))+&
                    q(i,1-2+0+0+0,indvars(3))*q(i,1-2+0+0+0,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_4_0_i1m2p0p0k_i1m2p0p0p1k = ((q(i,1-2+0+0+1,indvars(4))-&
                    0.5_wp*(q(i,1-2+0+0+1,indvars(2))*q(i,1-2+0+0+1,indvars(2))+&
                    q(i,1-2+0+0+1,indvars(3))*q(i,1-2+0+0+1,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_4_0_i1m2p0p0k_i1m2p0p0p2k = ((q(i,1-2+0+0+2,indvars(4))-&
                    0.5_wp*(q(i,1-2+0+0+2,indvars(2))*q(i,1-2+0+0+2,indvars(2))+&
                    q(i,1-2+0+0+2,indvars(3))*q(i,1-2+0+0+2,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_4_0_i1m2p0p0k = -&
          1.5_wp*d2_rhs_et_dydy_4_0_i1m2p0p0k_i1m2p0p0p0k+&
          2.0_wp*d2_rhs_et_dydy_4_0_i1m2p0p0k_i1m2p0p0p1k-&
          0.5_wp*d2_rhs_et_dydy_4_0_i1m2p0p0k_i1m2p0p0p2k

d2_rhs_et_dydy_4_0_i1m2p0p0k = d2_rhs_et_dydy_4_0_i1m2p0p0k*param_float(2)

d2_rhs_et_dydy_4_0_i1m2p0p1k_i1m2p0p1m1k = ((q(i,1-2+0+1-1,indvars(4))-&
                    0.5_wp*(q(i,1-2+0+1-1,indvars(2))*q(i,1-2+0+1-1,indvars(2))+&
                    q(i,1-2+0+1-1,indvars(3))*q(i,1-2+0+1-1,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_4_0_i1m2p0p1k_i1m2p0p1p1k = ((q(i,1-2+0+1+1,indvars(4))-&
                    0.5_wp*(q(i,1-2+0+1+1,indvars(2))*q(i,1-2+0+1+1,indvars(2))+&
                    q(i,1-2+0+1+1,indvars(3))*q(i,1-2+0+1+1,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_4_0_i1m2p0p1k = -&
          0.5_wp*d2_rhs_et_dydy_4_0_i1m2p0p1k_i1m2p0p1m1k+&
          0.5_wp*d2_rhs_et_dydy_4_0_i1m2p0p1k_i1m2p0p1p1k

d2_rhs_et_dydy_4_0_i1m2p0p1k = d2_rhs_et_dydy_4_0_i1m2p0p1k*param_float(2)

d2_rhs_et_dydy_4_0_i1m2p0p2k_i1m2p0p2m1k = ((q(i,1-2+0+2-1,indvars(4))-&
                    0.5_wp*(q(i,1-2+0+2-1,indvars(2))*q(i,1-2+0+2-1,indvars(2))+&
                    q(i,1-2+0+2-1,indvars(3))*q(i,1-2+0+2-1,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_4_0_i1m2p0p2k_i1m2p0p2p1k = ((q(i,1-2+0+2+1,indvars(4))-&
                    0.5_wp*(q(i,1-2+0+2+1,indvars(2))*q(i,1-2+0+2+1,indvars(2))+&
                    q(i,1-2+0+2+1,indvars(3))*q(i,1-2+0+2+1,indvars(3)))))/(param_float(4 + 5))

d2_rhs_et_dydy_4_0_i1m2p0p2k = -&
          0.5_wp*d2_rhs_et_dydy_4_0_i1m2p0p2k_i1m2p0p2m1k+&
          0.5_wp*d2_rhs_et_dydy_4_0_i1m2p0p2k_i1m2p0p2p1k

d2_rhs_et_dydy_4_0_i1m2p0p2k = d2_rhs_et_dydy_4_0_i1m2p0p2k*param_float(2)

d1_rhs_et_dy_4_i1m2p0p0k = (d2_rhs_et_dydy_4_0_i1m2p0p0k)*qst(i,1-2+0+0,indvarsst(11))

d1_rhs_et_dy_4_i1m2p0p1k = (d2_rhs_et_dydy_4_0_i1m2p0p1k)*qst(i,1-2+0+1,indvarsst(11))

d1_rhs_et_dy_4_i1m2p0p2k = (d2_rhs_et_dydy_4_0_i1m2p0p2k)*qst(i,1-2+0+2,indvarsst(11))

d1_rhs_et_dy_4_i1m2p0k = -&
          1.5_wp*d1_rhs_et_dy_4_i1m2p0p0k+&
          2.0_wp*d1_rhs_et_dy_4_i1m2p0p1k-&
          0.5_wp*d1_rhs_et_dy_4_i1m2p0p2k

d1_rhs_et_dy_4_i1m2p0k = d1_rhs_et_dy_4_i1m2p0k*param_float(2)



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None d(rho et)/dt *******
!                                                           
!***********************************************************


rhs(i,1-2+0,indvars(4)) =   -  ( (qst(i,1-2+0,indvarsst(10))*(d1_rhs_et_dx_0_i1m2p0k)+&
                    qst(i,1-2+0,indvarsst(10))*(d1_rhs_et_dx_1_i1m2p0k))*qst(i,1-2+0,indvarsst(5))+&
                    (-&
                    (1.0_wp)*(1.0_wp+&
                    ((q(i,1-2+0,indvars(5))/1.0_wp)**3.0_wp/((q(i,1-2+0,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i,1-2+0,indvars(5))/1.0_wp))*param_float(1 + 5)*(d1_rhs_et_dy_0_i1m2p0k)*(d1_rhs_et_dy_0_i1m2p0k)*qst(i,1-2+0,indvarsst(11))**2-&
                    4.0_wp/3.0_wp*(1.0_wp)*(1.0_wp+&
                    ((q(i,1-2+0,indvars(5))/1.0_wp)**3.0_wp/((q(i,1-2+0,indvars(5))/1.0_wp)**3.0_wp+&
                    param_float(13 + 5)**3.0_wp))*(q(i,1-2+0,indvars(5))/1.0_wp))*param_float(1 + 5)*(d1_rhs_rho_dy_0_i1m2p0k)*(d1_rhs_rho_dy_0_i1m2p0k)*qst(i,1-2+0,indvarsst(11))**2-&
                    param_float(2 + 5)*(d1_rhs_et_dy_4_i1m2p0k)-&
                    param_float(2 + 5)*(d1_rhs_et_dx_2_i1m2p0k)+&
                    (q(i,1-2+0,indvars(1))*q(i,1-2+0,indvars(4))+&
                    (param_float(3 + 5))*q(i,1-2+0,indvars(1))*((q(i,1-2+0,indvars(4))-&
                    0.5_wp*(q(i,1-2+0,indvars(2))*q(i,1-2+0,indvars(2))+&
                    q(i,1-2+0,indvars(3))*q(i,1-2+0,indvars(3))))))*(d1_rhs_rho_dy_0_i1m2p0k)*qst(i,1-2+0,indvarsst(11)))*dabs(1.0_wp-&
                    qst(i,1-2+0,indvarsst(5))) ) 



!***********************************************************
!                                                           
! building source terms in RHS for layer None 0 None d(rho nut)/dt 
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (deltaxI*([-ReI*(1.0_wp+chi)*sigmaI*deltaxI*({nut}_1x)]_1x)-ReI*Cb2*sigmaI*((deltaxI)**2*([rho*nut]_1x)*([nut]_1x))-Cb1*(1-ft2)*SS*rho*nut+ReI*(Cw1*fw-Cb1/k**2*ft2)*rho*(nut/eta)**2+deltaxI*([(rho*nut)*u]_1x))*symm
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d2_rhs_nut_dxdx_0_0_im11m2p0k_im1m11m2p0k = q(i-1-1,1-2+0,indvars(5))

d2_rhs_nut_dxdx_0_0_im11m2p0k_im1p11m2p0k = q(i-1+1,1-2+0,indvars(5))

d2_rhs_nut_dxdx_0_0_im11m2p0k = -&
          0.5_wp*d2_rhs_nut_dxdx_0_0_im11m2p0k_im1m11m2p0k+&
          0.5_wp*d2_rhs_nut_dxdx_0_0_im11m2p0k_im1p11m2p0k

d2_rhs_nut_dxdx_0_0_im11m2p0k = d2_rhs_nut_dxdx_0_0_im11m2p0k*param_float(1)

d2_rhs_nut_dxdx_0_0_ip11m2p0k_ip1m11m2p0k = q(i+1-1,1-2+0,indvars(5))

d2_rhs_nut_dxdx_0_0_ip11m2p0k_ip1p11m2p0k = q(i+1+1,1-2+0,indvars(5))

d2_rhs_nut_dxdx_0_0_ip11m2p0k = -&
          0.5_wp*d2_rhs_nut_dxdx_0_0_ip11m2p0k_ip1m11m2p0k+&
          0.5_wp*d2_rhs_nut_dxdx_0_0_ip11m2p0k_ip1p11m2p0k

d2_rhs_nut_dxdx_0_0_ip11m2p0k = d2_rhs_nut_dxdx_0_0_ip11m2p0k*param_float(1)

d1_rhs_nut_dx_1_im11m2p0k = q(i-1,1-2+0,indvars(1))*q(i-1,1-2+0,indvars(5))

d1_rhs_nut_dx_1_ip11m2p0k = q(i+1,1-2+0,indvars(1))*q(i+1,1-2+0,indvars(5))

d1_rhs_nut_dx_1_i1m2p0k = -&
          0.5_wp*d1_rhs_nut_dx_1_im11m2p0k+&
          0.5_wp*d1_rhs_nut_dx_1_ip11m2p0k

d1_rhs_nut_dx_1_i1m2p0k = d1_rhs_nut_dx_1_i1m2p0k*param_float(1)

d1_rhs_nut_dx_2_im11m2p0k = q(i-1,1-2+0,indvars(5))

d1_rhs_nut_dx_2_ip11m2p0k = q(i+1,1-2+0,indvars(5))

d1_rhs_nut_dx_2_i1m2p0k = -&
          0.5_wp*d1_rhs_nut_dx_2_im11m2p0k+&
          0.5_wp*d1_rhs_nut_dx_2_ip11m2p0k

d1_rhs_nut_dx_2_i1m2p0k = d1_rhs_nut_dx_2_i1m2p0k*param_float(1)

d1_rhs_nut_dx_3_im11m2p0k = (q(i-1,1-2+0,indvars(1))*q(i-1,1-2+0,indvars(5)))*q(i-1,1-2+0,indvars(2))

d1_rhs_nut_dx_3_ip11m2p0k = (q(i+1,1-2+0,indvars(1))*q(i+1,1-2+0,indvars(5)))*q(i+1,1-2+0,indvars(2))

d1_rhs_nut_dx_3_i1m2p0k = -&
          0.5_wp*d1_rhs_nut_dx_3_im11m2p0k+&
          0.5_wp*d1_rhs_nut_dx_3_ip11m2p0k

d1_rhs_nut_dx_3_i1m2p0k = d1_rhs_nut_dx_3_i1m2p0k*param_float(1)

d1_rhs_nut_dx_0_im11m2p0k = -param_float(1 + 5)*(1.0_wp+&
                    (q(i-1,1-2+0,indvars(5))/1.0_wp))*param_float(18 + 5)*qst(i-1,1-2+0,indvarsst(10))*(d2_rhs_nut_dxdx_0_0_im11m2p0k)

d1_rhs_nut_dx_0_ip11m2p0k = -param_float(1 + 5)*(1.0_wp+&
                    (q(i+1,1-2+0,indvars(5))/1.0_wp))*param_float(18 + 5)*qst(i+1,1-2+0,indvarsst(10))*(d2_rhs_nut_dxdx_0_0_ip11m2p0k)

d1_rhs_nut_dx_0_i1m2p0k = -&
          0.5_wp*d1_rhs_nut_dx_0_im11m2p0k+&
          0.5_wp*d1_rhs_nut_dx_0_ip11m2p0k

d1_rhs_nut_dx_0_i1m2p0k = d1_rhs_nut_dx_0_i1m2p0k*param_float(1)



!***********************************************************
!                                                           
! Update BC terms for layer None 0 None d(rho nut)/dt ******
!                                                           
!***********************************************************


rhs(i,1-2+0,indvars(5)) =   -  ( (qst(i,1-2+0,indvarsst(10))*(d1_rhs_nut_dx_0_i1m2p0k)-&
                    param_float(1 + 5)*param_float(7 + 5)*param_float(18 + 5)*((qst(i,1-2+0,indvarsst(10)))**2*(d1_rhs_nut_dx_1_i1m2p0k)*(d1_rhs_nut_dx_2_i1m2p0k))-&
                    param_float(6 + 5)*(1-&
                    (param_float(16 + 5)*exp(-&
                    param_float(17 + 5)*(q(i,1-2+0,indvars(5))/1.0_wp)**2.0_wp)))*qst(i,1-2+0,indvarsst(13))*q(i,1-2+0,indvars(1))*q(i,1-2+0,indvars(5))+&
                    param_float(1 + 5)*(param_float(10 + 5)*qst(i,1-2+0,indvarsst(16))-&
                    param_float(6 + 5)/param_float(9 + 5)**2*(param_float(16 + 5)*exp(-&
                    param_float(17 + 5)*(q(i,1-2+0,indvars(5))/1.0_wp)**2.0_wp)))*q(i,1-2+0,indvars(1))*(q(i,1-2+0,indvars(5))/qst(i,1-2+0,indvarsst(2)))**2+&
                    qst(i,1-2+0,indvarsst(10))*(d1_rhs_nut_dx_3_i1m2p0k))*qst(i,1-2+0,indvarsst(5)) ) 

   enddo
