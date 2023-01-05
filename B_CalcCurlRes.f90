Subroutine B_CalcCurlRes(NI,NJ,curlV,curlV_t,curlV_res)
implicit none
integer:: ni,nj
real,dimension(0:NI,0:NJ):: curlV,curlV_t,curlV_res

curlV_res = 0
curlV_res(1:NI-1,1:NJ-1) = abs(curlV(1:NI-1,1:NJ-1)-curlV_t(1:NI-1,1:NJ-1))/abs(curlV_t(1:NI-1,1:NJ-1))
write(*,*) 'Max error, curl:', maxval(curlV_res)


End Subroutine 
