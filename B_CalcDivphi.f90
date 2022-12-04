Subroutine B_CalcDivPhi(NI,NJ,X,Y,V,phi,gradPhi,Div,CellVolume,CellCenter,    &
&											IFaceVector,JFaceVector,  &
&											IFaceCenter,JFaceCenter)
implicit none
 integer::NI,NJ,I,J,sch=1
 real,dimension(NI,NJ):: X,Y
 real,dimension(NI-1,NJ-1):: CellVolume ! scalar arrays
 real,dimension(0:NI,0:NJ,2):: V,gradPhi
 real:: CellCenter(0:NI,0:NJ,2), &
	   &IFaceCenter( NI,NJ-1,2),IFaceVector( NI,NJ-1,2), &
	   &JFaceCenter( NI-1,NJ,2),JFaceVector( NI-1,NJ,2)			! vector arrays
 real,dimension(0:NI,0:NJ):: phi,Div
 

real:: S(2),N(2),W(2),E(2),&			!Face values
&	   r1(2),r2(2),        &			!Vectors from center to face
&	   distN(2),distE(2)    			!Dist from center to face

Div=0

 !В два обхода по всем ячейкам - по i и по j направлениям
 do j=1,NJ-1
	E = V(0,j,:)
	do i=1,NI-1
		!Интерполяция на грани S->1, N->2, W->3, E->4
		r1 = IFaceCenter(i+1,j,:) - CellCenter(i,j,:)
		r2 = CellCenter(i+1,j,:) - IFaceCenter(i+1,j,:)
		distE(1) = norm2(r1) 
		distE(2) = norm2(r2)
		
		W = E
		!Без скаляра
		E = (V(i+1,j,:)*distE(1) + V(i,j,:)*distE(2))/sum(distE)
		select case (sch)
		!Центральная схема
			case (1)
				E = E*(phi(i+1,j)*distE(1) + phi(i,j)*distE(2))/sum(distE)
			
			!FOU
			case (2)
				if (dot_product(E,IFaceVector(i+1,j,:)) >= 0) then
					E = E*phi(i,j)
				else
					E = E*phi(i+1,j)
				end if
			
			!SOU
			case (3)
				if (dot_product(E,IFaceVector(i+1,j,:)) >= 0) then
					E = E * (phi(i,j) + &
		&					dot_product(gradPhi(i,j,:),r1))
				else
					E = E * (phi(i+1,j) + &
		&					dot_product(gradPhi(i+1,j,:),r2))
				end if
		    case default
				print*, 'Error, select interpolation scheme'
			
		end select
		
		Div(i,j) = Div(i,j) + 1/CellVolume(i,j) * &
&			(-dot_product(W,IFaceVector(i,j,:)) + dot_product(E,IFaceVector(i+1,j,:))) 	
	end do
 end do
 
  do i=1,NI-1
	N = V(i,0,:)
	do j=1,NJ-1
		!Интерполяция на грани S->1, N->2, W->3, E->4
		r1 = JFaceCenter(i,j+1,:) -  CellCenter(i,j,:)
		r2 = CellCenter(i,j+1,:) - JFaceCenter(i,j+1,:)
		distN(1) = norm2(r1)
		distN(2) = norm2(r2)
		
		S = N
		!Без скаляра
		N = (V(i,j+1,:)*distN(1) + V(i,j,:)*distN(2))/sum(distN)
		select case (sch)
		!Центральная схема
			case (1)
				N = N*(phi(i,j+1)*distN(1) + phi(i,j)*distN(2))/sum(distN)
			
			!FOU
			case (2)
				if (dot_product(N,IFaceVector(i+1,j,:)) >= 0) then
					N = N*phi(i,j)
				else
					N = N*phi(i,j+1)
				end if
			
			!SOU
			case (3)
				if (dot_product(N,IFaceVector(i+1,j,:)) >= 0) then
					N = N * (phi(i,j) + &
		&					dot_product(gradPhi(i,j,:),r1))
				else
					N = N * (phi(i,j+1) + &
		&					dot_product(gradPhi(i,j+1,:),r2))
				end if
		    case default
				print*, 'Error, select interpolation scheme'
			
		end select

		
		Div(i,j) = Div(i,j) + 1/CellVolume(i,j) * &
&			(-dot_product(S,JFaceVector(i,j,:)) + dot_product(N,JFaceVector(i,j+1,:))) 	
	end do
 end do
 

 
End Subroutine 
