Subroutine B_CalcGradient(NI,NJ,X,Y,P,GradP,CellVolume,CellCenter,    &
&											IFaceVector,JFaceVector,  &
&											IFaceCenter,JFaceCenter)
implicit none
 integer::NI,NJ,I,J
 real,dimension(NI,NJ):: X,Y
 real,dimension(NI-1,NJ-1):: CellVolume ! scalar arrays
 real,dimension(0:NI,0:NJ):: P
 real:: CellCenter(0:NI,0:NJ,2), &
	   &IFaceCenter( NI,NJ-1,2),IFaceVector( NI,NJ-1,2), &
	   &JFaceCenter( NI-1,NJ,2),JFaceVector( NI-1,NJ,2)			! vector arrays
 real,dimension(0:NI,0:NJ,2):: GradP
 real,dimension(0:NI,0:NJ,2):: GradP_old
 

real:: S,N,W,E,&			!Face values
&	   distN(2),distE(2),&  !Dist from center to face
&	   EM(2),gradPE(2)		!Scew correction dist and gradient

 gradP_old = gradP
 gradP = 0
 !В два обхода по всем ячейкам - по i и по j направлениям
 do j=1,NJ-1
	E = P(0,j)
	do i=1,NI-1
		!Интерполяция на грани S->1, N->2, W->3, E->4
		distE(1) = norm2(IFaceCenter(i+1,j,:) - CellCenter(i,j,:)) 
		distE(2) = norm2(CellCenter(i+1,j,:) - IFaceCenter(i+1,j,:))
		
		W = E
		E = (P(i+1,j)*distE(1) + P(i,j)*distE(2))/sum(distE)
		!Коррекция для скошенных ячеек
		EM = IFaceCenter(i+1,j,:) - (CellCenter(i+1,j,:)*distE(1) + CellCenter(i,j,:)*distE(2))/sum(distE)
		gradPE = (gradP_old(i+1,j,:)*distE(1) + gradP_old(i,j,:)*distE(2))/sum(distE)
		E = E + dot_product(EM,gradPE)
		
		
		GradP(i,j,:) = GradP(i,j,:) + 1/CellVolume(i,j) * &
&			(-W*IFaceVector(i,j,:) + E*IFaceVector(i+1,j,:)) 	
	end do
 end do
 
  do i=1,NI-1
	N = P(i,0)
	do j=1,NJ-1
		!Интерполяция на грани S->1, N->2, W->3, E->4
		distN(1) = norm2(JFaceCenter(i,j+1,:) -  CellCenter(i,j,:))
		distN(2) = norm2(CellCenter(i,j+1,:) - JFaceCenter(i,j+1,:))
		
		S = N
		N = (P(i,j+1)*distN(1) + P(i,j)*distN(2))/sum(distN)
		!Коррекция для скошенных ячеек
		EM = JFaceCenter(i,j+1,:) - (CellCenter(i,j+1,:)*distN(1) + CellCenter(i,j,:)*distN(2))/sum(distN)
		gradPE = (gradP_old(i,j+1,:)*distN(1) + gradP_old(i,j,:)*distN(2))/sum(distN)
		N = N + dot_product(EM,gradPE)
		
		GradP(i,j,:) = GradP(i,j,:) + 1/CellVolume(i,j) * &
&			(-S*JFaceVector(i,j,:) + N*JFaceVector(i,j+1,:)) 	
	end do
 end do
 

 
End Subroutine 
