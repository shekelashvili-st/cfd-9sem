Subroutine B_CalcDivPhi(NI,NJ,X,Y,V,phi,gradPhi,Div,CellVolume,CellCenter,    &
&											IFaceVector,JFaceVector,  &
&											IFaceCenter,JFaceCenter, sch)
implicit none
 integer::NI,NJ,I,J,sch
 real,dimension(NI,NJ):: X,Y
 real,dimension(NI-1,NJ-1):: CellVolume ! scalar arrays
 real,dimension(0:NI,0:NJ,2):: V,gradPhi
 real:: CellCenter(0:NI,0:NJ,2), &
	   &IFaceCenter( NI,NJ-1,2),IFaceVector( NI,NJ-1,2), &
	   &JFaceCenter( NI-1,NJ,2),JFaceVector( NI-1,NJ,2)			! vector arrays
 real,dimension(0:NI,0:NJ):: phi,Div
 

real:: S(2),N(2),W(2),E(2),&			!Face values
&	   r1(2),r2(2),        &			!Vectors from center to face
&	   distN(2),distE(2),  &  			!Dist from center to face
&	   phi0, grad0(2)

Div=0

 !В два обхода по всем ячейкам - по i и по j направлениям
 do j=1,NJ-1
!Начальная грань 
	E = V(0,j,:)
	select case (sch)
	!Центральная схема
	case (1)
		E = E * phi(0,j)
	!FOU	
	case (2)
		if (dot_product(E,IFaceVector(1,j,:)) >= 0) then
			phi0 = 2*phi(0,j) - phi(1,j)
			E = E*phi0
		else
			E = E*phi(1,j)
		end if		
	!SOU
	case (3)
		!Расчёт и интерполяция градиента в заграничную ячейку // простой способ, 1 порядок точности на границе
		!Линейная экстраполяция из первых двух ячеек i=1 и i=2
		r2 = -CellCenter(1,j,:) + IFaceCenter(1,j,:)
		r1 = -r2
		grad0 = gradPhi(2,j,:) + &
&				norm2(-CellCenter(1,j,:)-CellCenter(2,j,:))/norm2(CellCenter(1,j,:)-CellCenter(2,j,:)) * (gradPhi(1,j,:) - gradPhi(2,j,:))
		if (dot_product(E,IFaceVector(1,j,:)) >= 0) then
			E = E * (phi(0,j) + &
&					dot_product(grad0,r1))
		else
			E = E * (phi(1,j) + &
&					dot_product(gradPhi(1,j,:),r2))
		end if		
		
	case default
		print*, 'Error, select interpolation scheme'						
	end select
	do i=1,NI-2
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
	!Последняя итерация на границе
		
	W = E
	!Без скаляра
	E = V(ni,j,:)
	select case (sch)
	!Центральная схема
	case (1)
		E = E * phi(ni,j)
	!FOU	
	case (2)
		if (dot_product(E,IFaceVector(ni,j,:)) >= 0) then
			E = E*phi(ni-1,j)
		else
			phi0 = 2*phi(ni,j) - phi(ni-1,j)
			E = E*phi0
		end if		
	!SOU
	case (3)
		!Расчёт и интерполяция градиента в заграничную ячейку // простой способ, 1 порядок точности на границе
		!Линейная экстраполяция из последних двух ячеек
		r1 = IFaceCenter(ni,j,:) - CellCenter(ni-1,j,:)
		r2 = -r1
		phi0 = 2*phi(ni,j) - phi(ni-1,j)
		grad0 = gradPhi(ni-2,j,:) + &
&				norm2(CellCenter(ni-1,j,:)+2*r1-CellCenter(ni-2,j,:))/norm2(CellCenter(ni-1,j,:)-CellCenter(ni-2,j,:)) &
&				* (gradPhi(ni-1,j,:) - gradPhi(ni-2,j,:))
		if (dot_product(E,IFaceVector(1,j,:)) >= 0) then
			E = E * (phi(ni-1,j) + &
&					dot_product(gradPhi(ni-1,j,:),r1))
		else
			E = E * (phi0 + &
&					dot_product(grad0,r2))
		end if		
		
	case default
		print*, 'Error, select interpolation scheme'
			
	end select
		
	Div(ni-1,j) = Div(ni-1,j) + 1/CellVolume(ni-1,j) * &
&		(-dot_product(W,IFaceVector(ni-1,j,:)) + dot_product(E,IFaceVector(ni,j,:))) 	
		
 end do
 
  do i=1,NI-1
!Начальная грань 
	N = V(i,0,:)
	select case (sch)
	!Центральная схема
	case (1)
		N = N * phi(i,0)
	!FOU	
	case (2)
		if (dot_product(N,JFaceVector(i,1,:)) >= 0) then
			phi0 = 2*phi(i,0) - phi(i,1)
			N = N*phi0
		else
			N = N*phi(i,1)
		end if		
	!SOU
	case (3)
		!Расчёт и интерполяция градиента в заграничную ячейку // простой способ, 1 порядок точности на границе
		!Линейная экстраполяция из первых двух ячеек i=1 и i=2
		r2 = -CellCenter(i,1,:) + JFaceCenter(i,1,:)
		r1 = -r2
		phi0 = 2*phi(i,0) - phi(i,1)
		grad0 = gradPhi(i,2,:) + &
&				norm2(CellCenter(i,1,:)-2*r1-CellCenter(i,2,:))/norm2(CellCenter(i,1,:)-CellCenter(i,2,:)) &
&			 * (gradPhi(i,1,:) - gradPhi(i,2,:))
		if (dot_product(N,JFaceVector(i,1,:)) >= 0) then
			N = N * (phi0 + &
&					dot_product(grad0,r1))
		else
			N = N * (phi(i,1) + &
&					dot_product(gradPhi(i,1,:),r2))
		end if		
		
	case default
		print*, 'Error, select interpolation scheme'						
	end select
	do j=1,NJ-2
		!Интерполяция на грани S->1, N->2, W->3, E->4
		r1 = JFaceCenter(i,j+1,:) -  CellCenter(i,j,:)
		r2 = -CellCenter(i,j+1,:) + JFaceCenter(i,j+1,:)
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
	!Последняя итерация на границе
		
	S = N
	!Без скаляра
	N = V(i,nj,:)
	select case (sch)
	!Центральная схема
	case (1)
		N = N * phi(i,nj)
	!FOU	
	case (2)
		if (dot_product(N,JFaceVector(i,nj,:)) >= 0) then
			N = N*phi(i,nj-1)
		else
			phi0 = 2*phi(i,nj) - phi(i,nj-1)
			N = N*phi0
		end if		
	!SOU
	case (3)
		!Расчёт и интерполяция градиента в заграничную ячейку // простой способ, 1 порядок точности на границе
		!Линейная экстраполяция из последних двух ячеек
		r1 = JFaceCenter(i,nj,:) - CellCenter(i,nj-1,:)
		r2 = -r1
		phi0 = 2*phi(i,nj) - phi(i,nj-1)
		grad0 = gradPhi(i,nj-2,:) + &
&				norm2(CellCenter(i,nj-1,:)+2*r1-CellCenter(i,nj-2,:))/norm2(CellCenter(i,nj-1,:)-CellCenter(i,nj-2,:)) &
&				* (gradPhi(i,nj-1,:) - gradPhi(i,nj-2,:))
		if (dot_product(N,JFaceVector(1,j,:)) >= 0) then
			N = N * (phi(i,nj-1) + &
&					dot_product(gradPhi(i,nj-1,:),r1))
		else
			N = N * (phi0 + &
&					dot_product(grad0,r2))
		end if		
		
	case default
		print*, 'Error, select interpolation scheme'
			
	end select
		
	Div(i,nj-1) = Div(i,nj-1) + 1/CellVolume(i,nj-1) * &
&		(-dot_product(S,JFaceVector(i,nj-1,:)) + dot_product(N,JFaceVector(i,nj,:))) 	

 end do
 

 
End Subroutine 
