 
      do i=idloop(1),idloop(2) 

q(i,ny+2+0,1) = q(i,ny+2+0,1) - param_float(5)*q2(i,ny+2+0,1)
q(i,ny+2+0,2) = q(i,ny+2+0,2) - param_float(5)*q2(i,ny+2+0,2)
q(i,ny+2+0,3) = q(i,ny+2+0,3) - param_float(5)*q2(i,ny+2+0,3)
q(i,ny+2+0,4) = q(i,ny+2+0,4) - param_float(5)*q2(i,ny+2+0,4)
q(i,ny+2+0,5) = q(i,ny+2+0,5) - param_float(5)*q2(i,ny+2+0,5)

q(i,ny+2-1,1) = q(i,ny+2-1,1) - param_float(5)*q2(i,ny+2-1,1)
q(i,ny+2-1,2) = q(i,ny+2-1,2) - param_float(5)*q2(i,ny+2-1,2)
q(i,ny+2-1,3) = q(i,ny+2-1,3) - param_float(5)*q2(i,ny+2-1,3)
q(i,ny+2-1,4) = q(i,ny+2-1,4) - param_float(5)*q2(i,ny+2-1,4)
q(i,ny+2-1,5) = q(i,ny+2-1,5) - param_float(5)*q2(i,ny+2-1,5)

   enddo
