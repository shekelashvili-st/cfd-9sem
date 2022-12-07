Subroutine B_CalcDivRes(NI,NJ,divV,divV_t,divV_res)
implicit none
integer:: ni,nj
real,dimension(0:NI,0:NJ):: divV,divV_t,divV_res

divV_res=0
divV_res(1:NI-1,1:NJ-1) = abs(divV(1:NI-1,1:NJ-1)-divV_t(1:NI-1,1:NJ-1))/abs(divV_t(1:NI-1,1:NJ-1))

write(*,*) 'Max error, divergence:', maxval(divV_res)


End Subroutine 
