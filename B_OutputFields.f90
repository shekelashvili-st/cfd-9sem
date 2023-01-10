Subroutine B_OutputFields(IO,NI,NJ,X,Y,P,V,GradP,GradP_res,divV,divV_res,lapP,lapP_res,divVP,divVP_res,curlV,curlV_res,T)
  Real,Dimension(NI,NJ):: X,Y
  Real,Dimension(0:NI,0:NJ)::P,divV,divV_res,lapP,lapP_res,divVP,divVP_res,curlV,curlV_res,T
  Real,Dimension(0:NI,0:NJ,2)::V,GradP,gradP_res
  character(:),allocatable::varss

  varss = 'VARIABLES = "X", "Y", "P", "U", "V", "GradPX", "GradPY", "GradPX_res", '
  varss = varss // '"GradPY_res", "DivV", "DivV_res", "LapP", "LapP_res", "DivVP", "DivVP_res", "curlV", "curlV_res", "T"'
  Write(IO,*) varss
  Write(IO,*) 'ZONE I=',NI,', J=',NJ,', DATAPACKING=BLOCK, VARLOCATION=([3-30]=CELLCENTERED)'
  Write(IO,'(100F14.7)') X(1:NI,1:NJ) 
  Write(IO,'(100F14.7)') Y(1:NI,1:NJ)
  Write(IO,'(100F14.7)') P(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') V(1:NI-1,1:NJ-1,1)
  Write(IO,'(100F14.7)') V(1:NI-1,1:NJ-1,2)
  Write(IO,'(100F14.7)') GradP(1:NI-1,1:NJ-1,1)
  Write(IO,'(100F14.7)') GradP(1:NI-1,1:NJ-1,2)
  Write(IO,'(100F14.7)') GradP_res(1:NI-1,1:NJ-1,1)
  Write(IO,'(100F14.7)') GradP_res(1:NI-1,1:NJ-1,2)
  Write(IO,'(100F14.7)') divV(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') divV_res(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') lapP(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') lapP_res(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') divVP(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') divVP_res(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') curlV(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') curlV_res(1:NI-1,1:NJ-1)
  Write(IO,'(100F14.7)') T(1:NI-1,1:NJ-1)
  
End Subroutine 
