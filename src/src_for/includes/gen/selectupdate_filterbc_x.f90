
 CASE (1)

   idloop(3) = idarray(3)

 CASE (2)

   idloop(4) = idarray(4)

 CASE (3)

   idloop(3) = idarray(3)

 CASE (4)

   idloop(4) = idarray(4)

 CASE (5)

        do j=idloop(3),idloop(4) 
   
   q(1-2+0,j,1) = q(1-2+0,j,1) - param_float(5)*q2(1-2+0,j,1)
   q(1-2+0,j,2) = q(1-2+0,j,2) - param_float(5)*q2(1-2+0,j,2)
   q(1-2+0,j,3) = q(1-2+0,j,3) - param_float(5)*q2(1-2+0,j,3)
   q(1-2+0,j,4) = q(1-2+0,j,4) - param_float(5)*q2(1-2+0,j,4)
   q(1-2+0,j,5) = q(1-2+0,j,5) - param_float(5)*q2(1-2+0,j,5)
   
   q(1-2+1,j,1) = q(1-2+1,j,1) - param_float(5)*q2(1-2+1,j,1)
   q(1-2+1,j,2) = q(1-2+1,j,2) - param_float(5)*q2(1-2+1,j,2)
   q(1-2+1,j,3) = q(1-2+1,j,3) - param_float(5)*q2(1-2+1,j,3)
   q(1-2+1,j,4) = q(1-2+1,j,4) - param_float(5)*q2(1-2+1,j,4)
   q(1-2+1,j,5) = q(1-2+1,j,5) - param_float(5)*q2(1-2+1,j,5)
   
        enddo

 CASE (6)

        do j=idloop(3),idloop(4) 
   
   q(nx+2+0,j,1) = q(nx+2+0,j,1) - param_float(5)*q2(nx+2+0,j,1)
   q(nx+2+0,j,2) = q(nx+2+0,j,2) - param_float(5)*q2(nx+2+0,j,2)
   q(nx+2+0,j,3) = q(nx+2+0,j,3) - param_float(5)*q2(nx+2+0,j,3)
   q(nx+2+0,j,4) = q(nx+2+0,j,4) - param_float(5)*q2(nx+2+0,j,4)
   q(nx+2+0,j,5) = q(nx+2+0,j,5) - param_float(5)*q2(nx+2+0,j,5)
   
   q(nx+2-1,j,1) = q(nx+2-1,j,1) - param_float(5)*q2(nx+2-1,j,1)
   q(nx+2-1,j,2) = q(nx+2-1,j,2) - param_float(5)*q2(nx+2-1,j,2)
   q(nx+2-1,j,3) = q(nx+2-1,j,3) - param_float(5)*q2(nx+2-1,j,3)
   q(nx+2-1,j,4) = q(nx+2-1,j,4) - param_float(5)*q2(nx+2-1,j,4)
   q(nx+2-1,j,5) = q(nx+2-1,j,5) - param_float(5)*q2(nx+2-1,j,5)
   
        enddo
