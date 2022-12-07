Subroutine B_CalcDivphiRes(NI,NJ,div,div_t,div_res)
implicit none
integer:: ni,nj
real,dimension(0:NI,0:NJ):: div,div_t,div_res

div_res = 0
div_res(1:NI-1,1:NJ-1) = abs(div(1:NI-1,1:NJ-1)-div_t(1:NI-1,1:NJ-1))/abs(div_t(1:NI-1,1:NJ-1))
write(*,*) 'Max error, divergence_phi:', maxval(div_res)


End Subroutine 
