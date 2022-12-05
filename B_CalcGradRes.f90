Subroutine B_CalcGradRes(NI,NJ,gradP,gradP_t,gradP_res)
implicit none
integer:: ni,nj
real,dimension(0:NI,0:NJ,2):: gradP, gradP_t, gradP_res

gradP_res(1:NI-1,1:NJ-1,:) = abs(gradP(1:NI-1,1:NJ-1,:)-gradP_t(1:NI-1,1:NJ-1,:))

write(*,*) 'Max error:', maxval(gradP_res)


End Subroutine 
