Subroutine B_CalcLap(NI,NJ,X,Y,p,Lap,CellVolume,CellCenter,    &
&											IFaceVector,JFaceVector,  &
&											IFaceCenter,JFaceCenter)
implicit none
 integer::NI,NJ,I,J
 real,dimension(NI,NJ):: X,Y
 real,dimension(NI-1,NJ-1):: CellVolume ! scalar arrays
 real,dimension(0:NI,0:NJ):: p
 real:: CellCenter(0:NI,0:NJ,2), &
	   &IFaceCenter( NI,NJ-1,2),IFaceVector( NI,NJ-1,2), &
	   &JFaceCenter( NI-1,NJ,2),JFaceVector( NI-1,NJ,2),  &
&		gradP(0:NI,0:NJ,2)								 ! vector arrays
 real,dimension(0:NI,0:NJ):: Lap
 

real:: S,N,W,E,&			!Face values
&	   distN,distE			!Dist from center to face

Lap=0

 !В два обхода по всем ячейкам - по i и по j направлениям
 do j=1,NJ-1
	distE = norm2(CellCenter(1,j,:) - IFaceCenter(1,j,:))
	E = (p(1,j)-p(0,j))/distE
	do i=1,NI-1
		!Расчёт производной на грани
		distE = norm2(CellCenter(i+1,j,:) - CellCenter(i,j,:))
		
		W = E
	    E = (p(i+1,j)-p(i,j))/distE !+ &
!&					  (IFaceVector(i+1,j,:)-&
!&				      (CellCenter(i+1,j,:)-CellCenter(i,j,:)))
		
		
		Lap(i,j) = Lap(i,j) + 1/CellVolume(i,j) * &
&			(-W*norm2(IFaceVector(i,j,:)) + E*norm2(IFaceVector(i+1,j,:))) 	
	end do
 end do
 
  do i=1,NI-1
	distN = norm2(CellCenter(i,1,:) - JFaceCenter(i,1,:))
	N = (p(i,1)-p(i,0))/distN
	do j=1,NJ-1
		!Интерполяция на грани S->1, N->2, W->3, E->4
		distN = norm2(CellCenter(i,j+1,:) - CellCenter(i,j,:))
		
		S = N
	    N = (p(i,j+1)-p(i,j))/distN
		
		
		Lap(i,j) = Lap(i,j) + 1/CellVolume(i,j) * &
&			(-S*norm2(JFaceVector(i,j,:)) + N*norm2(JFaceVector(i,j+1,:))) 	
	end do
 end do
 

 
End Subroutine 
