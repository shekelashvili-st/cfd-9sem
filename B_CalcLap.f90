Subroutine B_CalcLap(NI,NJ,X,Y,p,gradP,Lap,CellVolume,CellCenter,    &
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
 

real:: S,N,W,E,&					!Face values
&	   distN(2),distE(2),distC,&	!Dist from center to face
&	   r1(2),gradPE(2), gradPN(2)
Lap=0

 !В два обхода по всем ячейкам - по i и по j направлениям
 do j=1,NJ-1
	r1 = CellCenter(1,j,:) - IFaceCenter(1,j,:)
	distC = norm2(r1)
	r1 = r1/distC
	E = 5/3*(p(1,j)-p(0,j))/distC - 2/3*dot_product(gradP(1,j,:),IFaceVector(1,j,:)/norm2(IFaceVector(1,j,:))) + &
		dot_product((IFaceVector(1,j,:)/norm2(IFaceVector(1,j,:)) - r1),gradP(1,j,:))	
	
	do i=1,NI-2
		!Расчёт производной на грани
		distE(1) = norm2(IFaceCenter(i+1,j,:) - CellCenter(i,j,:)) 
		distE(2) = norm2(CellCenter(i+1,j,:) - IFaceCenter(i+1,j,:))
		r1 = CellCenter(i+1,j,:) - CellCenter(i,j,:)
		distC = norm2(r1)
		r1 = r1/distC
		
		W = E
		gradPE = (gradP(i+1,j,:)*distE(1) + gradP(i,j,:)*distE(2))/sum(distE)
	    E = (p(i+1,j)-p(i,j))/distC + &
&			dot_product((IFaceVector(i+1,j,:)/norm2(IFaceVector(i+1,j,:)) - r1),gradPE) 
		
		
		Lap(i,j) = Lap(i,j) + 1/CellVolume(i,j) * &
&			(-W*norm2(IFaceVector(i,j,:)) + E*norm2(IFaceVector(i+1,j,:)))
	end do
	!Обработка внешней границы
	W = E
	r1 = CellCenter(ni-1,j,:) - IFaceCenter(ni,j,:)
	distC = norm2(r1)
	r1 = r1/distC
	E = -(5/3*(p(ni-1,j)-p(ni,j))/distC - 2/3*dot_product(gradP(ni-1,j,:),-IFaceVector(ni,j,:)/norm2(IFaceVector(ni,j,:))) + &
	dot_product((-IFaceVector(ni,j,:)/norm2(IFaceVector(ni,j,:)) - r1),gradP(ni-1,j,:)))
	Lap(ni-1,j) = Lap(ni-1,j) + 1/CellVolume(ni-1,j) * &
&		(-W*norm2(IFaceVector(ni-1,j,:)) + E*norm2(IFaceVector(ni,j,:)))	
 end do
 
  do i=1,NI-1
	r1 = CellCenter(i,1,:) - JFaceCenter(i,1,:)
	distC = norm2(r1)
	r1 = r1/distC
	N = 5/3*(p(i,1)-p(i,0))/distC - 2/3*dot_product(gradP(i,1,:),JFaceVector(i,1,:)/norm2(JFaceVector(i,1,:))) + &
		dot_product((JFaceVector(i,1,:)/norm2(JFaceVector(i,1,:)) - r1),gradP(i,1,:))
	do j=1,NJ-2
		!Интерполяция на грани S->1, N->2, W->3, E->4
		distN(1) = norm2(JFaceCenter(i,j+1,:) -  CellCenter(i,j,:))
		distN(2) = norm2(CellCenter(i,j+1,:) - JFaceCenter(i,j+1,:))
		r1 = CellCenter(i,j+1,:) - CellCenter(i,j,:)
		distC = norm2(r1)
		r1 = r1/distC
		
		S = N
		gradPN = (gradP(i,j+1,:)*distN(1) + gradP(i,j,:)*distN(2))/sum(distN)
	    N = (p(i,j+1)-p(i,j))/distC + &
&			dot_product((JFaceVector(i,j+1,:)/norm2(JFaceVector(i,j+1,:)) - r1),gradPN) 
		
		
		Lap(i,j) = Lap(i,j) + 1/CellVolume(i,j) * &
&			(-S*norm2(JFaceVector(i,j,:)) + N*norm2(JFaceVector(i,j+1,:)))
	end do
	!Обработка внешней границы
	S = N
	r1 = CellCenter(i,nj-1,:) - JFaceCenter(i,nj,:)
	distC = norm2(r1)
	r1 = r1/distC
	N = -(5/3*(p(i,nj-1)-p(i,nj))/distC - 2/3*dot_product(gradP(i,nj-1,:),-JFaceVector(i,nj,:)/norm2(JFaceVector(i,nj,:))) + &
		dot_product((-JFaceVector(i,nj,:)/norm2(JFaceVector(i,nj,:)) - r1),gradP(i,nj-1,:)))
	
	Lap(i,nj-1) = Lap(i,nj-1) + 1/CellVolume(i,nj-1) * &
&		(-S*norm2(JFaceVector(i,nj-1,:)) + N*norm2(JFaceVector(i,nj,:)))
 end do
 

 
End Subroutine 
