
 CASE (1)

   idloop(1) = idarray(1)

 CASE (2)

   idloop(1) = idarray(1)

 CASE (3)

   idloop(2) = idarray(2)

 CASE (4)

   idloop(2) = idarray(2)

 CASE (7)

    
         do i=idloop(1),idloop(2) 
   
   q(i,1-4+0,1) = q(i,1-4+0,1) - param_float(5)*q2(i,1-4+0,1)
   q(i,1-4+0,2) = q(i,1-4+0,2) - param_float(5)*q2(i,1-4+0,2)
   q(i,1-4+0,3) = q(i,1-4+0,3) - param_float(5)*q2(i,1-4+0,3)
   q(i,1-4+0,4) = q(i,1-4+0,4) - param_float(5)*q2(i,1-4+0,4)
   q(i,1-4+0,5) = q(i,1-4+0,5) - param_float(5)*q2(i,1-4+0,5)
   
   q(i,1-4+1,1) = q(i,1-4+1,1) - param_float(5)*q2(i,1-4+1,1)
   q(i,1-4+1,2) = q(i,1-4+1,2) - param_float(5)*q2(i,1-4+1,2)
   q(i,1-4+1,3) = q(i,1-4+1,3) - param_float(5)*q2(i,1-4+1,3)
   q(i,1-4+1,4) = q(i,1-4+1,4) - param_float(5)*q2(i,1-4+1,4)
   q(i,1-4+1,5) = q(i,1-4+1,5) - param_float(5)*q2(i,1-4+1,5)
   
   q(i,1-4+2,1) = q(i,1-4+2,1) - param_float(5)*q2(i,1-4+2,1)
   q(i,1-4+2,2) = q(i,1-4+2,2) - param_float(5)*q2(i,1-4+2,2)
   q(i,1-4+2,3) = q(i,1-4+2,3) - param_float(5)*q2(i,1-4+2,3)
   q(i,1-4+2,4) = q(i,1-4+2,4) - param_float(5)*q2(i,1-4+2,4)
   q(i,1-4+2,5) = q(i,1-4+2,5) - param_float(5)*q2(i,1-4+2,5)
   
   q(i,1-4+3,1) = q(i,1-4+3,1) - param_float(5)*q2(i,1-4+3,1)
   q(i,1-4+3,2) = q(i,1-4+3,2) - param_float(5)*q2(i,1-4+3,2)
   q(i,1-4+3,3) = q(i,1-4+3,3) - param_float(5)*q2(i,1-4+3,3)
   q(i,1-4+3,4) = q(i,1-4+3,4) - param_float(5)*q2(i,1-4+3,4)
   q(i,1-4+3,5) = q(i,1-4+3,5) - param_float(5)*q2(i,1-4+3,5)
   
      enddo

 CASE (8)

    
         do i=idloop(1),idloop(2) 
   
   q(i,ny+4+0,1) = q(i,ny+4+0,1) - param_float(5)*q2(i,ny+4+0,1)
   q(i,ny+4+0,2) = q(i,ny+4+0,2) - param_float(5)*q2(i,ny+4+0,2)
   q(i,ny+4+0,3) = q(i,ny+4+0,3) - param_float(5)*q2(i,ny+4+0,3)
   q(i,ny+4+0,4) = q(i,ny+4+0,4) - param_float(5)*q2(i,ny+4+0,4)
   q(i,ny+4+0,5) = q(i,ny+4+0,5) - param_float(5)*q2(i,ny+4+0,5)
   
   q(i,ny+4-1,1) = q(i,ny+4-1,1) - param_float(5)*q2(i,ny+4-1,1)
   q(i,ny+4-1,2) = q(i,ny+4-1,2) - param_float(5)*q2(i,ny+4-1,2)
   q(i,ny+4-1,3) = q(i,ny+4-1,3) - param_float(5)*q2(i,ny+4-1,3)
   q(i,ny+4-1,4) = q(i,ny+4-1,4) - param_float(5)*q2(i,ny+4-1,4)
   q(i,ny+4-1,5) = q(i,ny+4-1,5) - param_float(5)*q2(i,ny+4-1,5)
   
   q(i,ny+4-2,1) = q(i,ny+4-2,1) - param_float(5)*q2(i,ny+4-2,1)
   q(i,ny+4-2,2) = q(i,ny+4-2,2) - param_float(5)*q2(i,ny+4-2,2)
   q(i,ny+4-2,3) = q(i,ny+4-2,3) - param_float(5)*q2(i,ny+4-2,3)
   q(i,ny+4-2,4) = q(i,ny+4-2,4) - param_float(5)*q2(i,ny+4-2,4)
   q(i,ny+4-2,5) = q(i,ny+4-2,5) - param_float(5)*q2(i,ny+4-2,5)
   
   q(i,ny+4-3,1) = q(i,ny+4-3,1) - param_float(5)*q2(i,ny+4-3,1)
   q(i,ny+4-3,2) = q(i,ny+4-3,2) - param_float(5)*q2(i,ny+4-3,2)
   q(i,ny+4-3,3) = q(i,ny+4-3,3) - param_float(5)*q2(i,ny+4-3,3)
   q(i,ny+4-3,4) = q(i,ny+4-3,4) - param_float(5)*q2(i,ny+4-3,4)
   q(i,ny+4-3,5) = q(i,ny+4-3,5) - param_float(5)*q2(i,ny+4-3,5)
   
      enddo
