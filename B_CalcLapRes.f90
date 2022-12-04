Subroutine B_CalcLapRes(NI,NJ,lapP,lapP_t,lapP_res)
implicit none
integer:: ni,nj
real,dimension(0:NI,0:NJ):: lapP,lapP_t,lapP_res

lapP_res = abs(lapP(1:NI-1,1:NJ-1)-lapP_t(1:NI-1,1:NJ-1))

write(*,*) 'Max error, laplacian:', maxval(lapP_res)


End Subroutine 
