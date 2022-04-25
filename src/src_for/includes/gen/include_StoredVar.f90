

!***********************************************************
!                                                           
! Start computing stored variables *************************
!                                                           
!***********************************************************


     do j=idloop(3),idloop(4) 
 
      do i=idloop(1),idloop(2) 


!***********************************************************
!                                                           
! building stored variable stemp ***************************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (((0.5_wp*(deltayI*([u]_1y)-deltaxI*([v]_1x)))**2+(0.5_wp*(deltaxI*([v]_1x)-deltayI*([u]_1y)))**2)*2)**0.5
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

d1_stemp_dx_0_im2jk = q(i-2,j,indvars(3))

d1_stemp_dx_0_im1jk = q(i-1,j,indvars(3))

d1_stemp_dx_0_ip1jk = q(i+1,j,indvars(3))

d1_stemp_dx_0_ip2jk = q(i+2,j,indvars(3))

d1_stemp_dx_0_ijk = 0.08333333333333333_wp*d1_stemp_dx_0_im2jk-&
          0.666666666666667_wp*d1_stemp_dx_0_im1jk+&
          0.666666666666667_wp*d1_stemp_dx_0_ip1jk-&
          0.08333333333333333_wp*d1_stemp_dx_0_ip2jk

d1_stemp_dx_0_ijk = d1_stemp_dx_0_ijk*param_float(1)

d1_stemp_dy_0_ijm2k = q(i,j-2,indvars(2))

d1_stemp_dy_0_ijm1k = q(i,j-1,indvars(2))

d1_stemp_dy_0_ijp1k = q(i,j+1,indvars(2))

d1_stemp_dy_0_ijp2k = q(i,j+2,indvars(2))

d1_stemp_dy_0_ijk = 0.08333333333333333_wp*d1_stemp_dy_0_ijm2k-&
          0.666666666666667_wp*d1_stemp_dy_0_ijm1k+&
          0.666666666666667_wp*d1_stemp_dy_0_ijp1k-&
          0.08333333333333333_wp*d1_stemp_dy_0_ijp2k

d1_stemp_dy_0_ijk = d1_stemp_dy_0_ijk*param_float(2)



!***********************************************************
!                                                           
! Update stored variables (( (  0.5_wp*( deltayI*( [u]_1y ) - deltaxI*( [v]_1x ) )  )**2 + (  0.5_wp*( deltaxI*( [v]_1x ) - deltayI*( [u]_1y ) )  )**2 )*2)**0.5 
!                                                           
!***********************************************************


qst(i,j,indvarsst(4)) =  (((0.5_wp*(qst(i,j,indvarsst(11))*(d1_stemp_dy_0_ijk)-&
                    qst(i,j,indvarsst(10))*(d1_stemp_dx_0_ijk)))**2+&
                    (0.5_wp*(qst(i,j,indvarsst(10))*(d1_stemp_dx_0_ijk)-&
                    qst(i,j,indvarsst(11))*(d1_stemp_dy_0_ijk)))**2)*2)**0.5



!***********************************************************
!                                                           
! building stored variable SS ******************************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (stemp+ReI*nut/(k**2*eta**2))
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update stored variables  ( stemp+ReI*nut/(k**2*eta**2) )  
!                                                           
!***********************************************************


qst(i,j,indvarsst(12)) =  (qst(i,j,indvarsst(4))+&
                    param_float(1 + 5)*q(i,j,indvars(5))/(param_float(9 + 5)**2*qst(i,j,indvarsst(2))**2))



!***********************************************************
!                                                           
! building stored variable tau_wall ************************
!                                                           
!***********************************************************


!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!
! (visc_t*([u]_1y)*deltayI)
!
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



!***********************************************************
!                                                           
! Update stored variables  ( visc_t*( [u]_1y )*deltayI ) ***
!                                                           
!***********************************************************


qst(i,j,indvarsst(13)) =  (((1+&
                    param_float(21 + 5))/(((q(i,j,indvars(4))-&
                    0.5_wp*(q(i,j,indvars(2))*q(i,j,indvars(2))+&
                    q(i,j,indvars(3))*q(i,j,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,j,indvars(4))-&
                    0.5_wp*(q(i,j,indvars(2))*q(i,j,indvars(2))+&
                    q(i,j,indvars(3))*q(i,j,indvars(3)))))/param_float(4 + 5)**1.5)*(1+&
                    ((q(i,j,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,j,indvars(4))-&
                    0.5_wp*(q(i,j,indvars(2))*q(i,j,indvars(2))+&
                    q(i,j,indvars(3))*q(i,j,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,j,indvars(4))-&
                    0.5_wp*(q(i,j,indvars(2))*q(i,j,indvars(2))+&
                    q(i,j,indvars(3))*q(i,j,indvars(3)))))/param_float(4 + 5)**1.5*q(i,j,indvars(1)))**3/((q(i,j,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,j,indvars(4))-&
                    0.5_wp*(q(i,j,indvars(2))*q(i,j,indvars(2))+&
                    q(i,j,indvars(3))*q(i,j,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,j,indvars(4))-&
                    0.5_wp*(q(i,j,indvars(2))*q(i,j,indvars(2))+&
                    q(i,j,indvars(3))*q(i,j,indvars(3)))))/param_float(4 + 5)**1.5*q(i,j,indvars(1)))**3+&
                    param_float(13 + 5)**3))*(q(i,j,indvars(5))/(1+&
                    param_float(21 + 5))/(((q(i,j,indvars(4))-&
                    0.5_wp*(q(i,j,indvars(2))*q(i,j,indvars(2))+&
                    q(i,j,indvars(3))*q(i,j,indvars(3)))))/param_float(4 + 5)+&
                    param_float(21 + 5))*((q(i,j,indvars(4))-&
                    0.5_wp*(q(i,j,indvars(2))*q(i,j,indvars(2))+&
                    q(i,j,indvars(3))*q(i,j,indvars(3)))))/param_float(4 + 5)**1.5*q(i,j,indvars(1))))*param_float(1 + 5)*(d1_stemp_dy_0_ijk)*qst(i,j,indvarsst(11)))

     enddo
   enddo
