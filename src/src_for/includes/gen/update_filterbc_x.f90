     do j=idloop(3),idloop(4) 

q(nx+4+0,j,1) = q(nx+4+0,j,1) - param_float(5)*q2(nx+4+0,j,1)
q(nx+4+0,j,2) = q(nx+4+0,j,2) - param_float(5)*q2(nx+4+0,j,2)
q(nx+4+0,j,3) = q(nx+4+0,j,3) - param_float(5)*q2(nx+4+0,j,3)
q(nx+4+0,j,4) = q(nx+4+0,j,4) - param_float(5)*q2(nx+4+0,j,4)
q(nx+4+0,j,5) = q(nx+4+0,j,5) - param_float(5)*q2(nx+4+0,j,5)

q(nx+4-1,j,1) = q(nx+4-1,j,1) - param_float(5)*q2(nx+4-1,j,1)
q(nx+4-1,j,2) = q(nx+4-1,j,2) - param_float(5)*q2(nx+4-1,j,2)
q(nx+4-1,j,3) = q(nx+4-1,j,3) - param_float(5)*q2(nx+4-1,j,3)
q(nx+4-1,j,4) = q(nx+4-1,j,4) - param_float(5)*q2(nx+4-1,j,4)
q(nx+4-1,j,5) = q(nx+4-1,j,5) - param_float(5)*q2(nx+4-1,j,5)

q(nx+4-2,j,1) = q(nx+4-2,j,1) - param_float(5)*q2(nx+4-2,j,1)
q(nx+4-2,j,2) = q(nx+4-2,j,2) - param_float(5)*q2(nx+4-2,j,2)
q(nx+4-2,j,3) = q(nx+4-2,j,3) - param_float(5)*q2(nx+4-2,j,3)
q(nx+4-2,j,4) = q(nx+4-2,j,4) - param_float(5)*q2(nx+4-2,j,4)
q(nx+4-2,j,5) = q(nx+4-2,j,5) - param_float(5)*q2(nx+4-2,j,5)

q(nx+4-3,j,1) = q(nx+4-3,j,1) - param_float(5)*q2(nx+4-3,j,1)
q(nx+4-3,j,2) = q(nx+4-3,j,2) - param_float(5)*q2(nx+4-3,j,2)
q(nx+4-3,j,3) = q(nx+4-3,j,3) - param_float(5)*q2(nx+4-3,j,3)
q(nx+4-3,j,4) = q(nx+4-3,j,4) - param_float(5)*q2(nx+4-3,j,4)
q(nx+4-3,j,5) = q(nx+4-3,j,5) - param_float(5)*q2(nx+4-3,j,5)

     enddo
