Subroutine B_CalcDiv(NI,NJ,X,Y,V,DivV,CellVolume,CellCenter,    &
&											IFaceVector,JFaceVector,  &
&											IFaceCenter,JFaceCenter)
implicit none
 integer::NI,NJ,I,J
 real,dimension(NI,NJ):: X,Y
 real,dimension(NI-1,NJ-1):: CellVolume ! scalar arrays
 real,dimension(0:NI,0:NJ,2):: V
 real:: CellCenter(0:NI,0:NJ,2), &
	   &IFaceCenter( NI,NJ-1,2),IFaceVector( NI,NJ-1,2), &
	   &JFaceCenter( NI-1,NJ,2),JFaceVector( NI-1,NJ,2)			! vector arrays
 real,dimension(0:NI,0:NJ):: DivV
 

real:: S(2),N(2),W(2),E(2),&			!Face values
&	   distN(2),distE(2)    			!Dist from center to face

DivV=0

 !В два обхода по всем ячейкам - по i и по j направлениям
 do j=1,NJ-1
	E = V(0,j,:)
	do i=1,NI-1
		!Интерполяция на грани S->1, N->2, W->3, E->4
		distE(1) = norm2(IFaceCenter(i+1,j,:) - CellCenter(i,j,:)) 
		distE(2) = norm2(CellCenter(i+1,j,:) - IFaceCenter(i+1,j,:))
		
		W = E
		E = (V(i+1,j,:)*distE(1) + V(i,j,:)*distE(2))/sum(distE)
		
		
		DivV(i,j) = DivV(i,j) + 1/CellVolume(i,j) * &
&			(-dot_product(W,IFaceVector(i,j,:)) + dot_product(E,IFaceVector(i+1,j,:))) 	
	end do
 end do
 
  do i=1,NI-1
	N = V(i,0,:)
	do j=1,NJ-1
		!Интерполяция на грани S->1, N->2, W->3, E->4
		distN(1) = norm2(JFaceCenter(i,j+1,:) -  CellCenter(i,j,:))
		distN(2) = norm2(CellCenter(i,j+1,:) - JFaceCenter(i,j+1,:))
		
		S = N
		N = (V(i,j+1,:)*distN(1) + V(i,j,:)*distN(2))/sum(distN)
		
		
		DivV(i,j) = DivV(i,j) + 1/CellVolume(i,j) * &
&			(-dot_product(S,JFaceVector(i,j,:)) + dot_product(N,JFaceVector(i,j+1,:))) 	
	end do
 end do
 

 
End Subroutine 
