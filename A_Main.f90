Program Main
  implicit none
  character(*), parameter:: InputFile='input.txt',OutputFile='data.plt',sol_file='fields.plt', &
							input_cavern='params.nml'							                ! names of input and output files
  character MeshFile*30       ! name of file with computational mesh
  integer::	i,j,ni,nj,m,nm,sch=3,sols = 0,cavern=0, iu
  integer, parameter:: IO = 12 ! input-output unit
  real,allocatable,dimension(:,:):: X,Y,P,CellVolume,DivV,DivV_t,DivV_res,lapP,lapP_t,lapP_res, &
&									DivVP,DivVP_t,DivVP_res,CurlV,CurlV_t,CurlV_res  	 ! scalar arrays
  real,allocatable,dimension(:,:,:):: CellCenter,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector,&
&									  GradP,GradP_t,GradP_res,V  ! vector arrays
!Решение уравнения переноса
  real:: Re, Pr, eps_T, eps, T1, T2, dt, rtmp, CFL
  real, allocatable:: T(:,:), T_new(:,:), gradT(:,:,:), conv(:,:), dif(:,:), res(:,:)

!===  READ INPUT FILE ===
  namelist /params_cavern/ Re, Pr, T1, T2, CFL, eps, nm
  open(newunit=iu, file=input_cavern)
  read(iu, nml=params_cavern)
  close(iu)
  
  WRITE(*,*) 'Read input file: ', InputFile
  OPEN(IO,FILE=InputFile)
  READ(IO,*) MeshFile  ! read name of file with computational mesh
  READ(IO,*) sch 	   ! 1 - linear, 2 - FOU, 3 - SOU
  READ(IO,*) sols 	   ! 0 - initialize fields from func, 1 - initialize from solution file
  READ(IO,*) cavern    ! 0 - без решения ур. переноса, 1 - с решением ур. переноса
  CLOSE(IO)

!===   READ NODES NUMBER (NI,NJ) FROM FILE WITH MESH ===
  WRITE(*,*) 'Read nodes number from file: ', MeshFile
  OPEN(IO,FILE = MeshFile)
  READ(IO,*) NI,NJ
  WRITE(*,*) 'NI, NJ = ',NI,NJ

!=== ALLOCATE ALL ARRAYS ===
  WRITE(*,*) 'Allocate arrays'       
  allocate(X(NI,NJ)) ! mesh nodes X-coordinates
  allocate(Y(NI,NJ)) ! mesh nodes Y-coordinates
  allocate(P(0:NI,0:NJ))   ! Pressure
  allocate(CellVolume(NI-1,NJ-1))   ! Cell Volumes    
  allocate(CellCenter(0:NI,0:NJ,2)) ! Cell Centers
  allocate(IFaceCenter( NI,NJ-1,2)) ! Face Centers for I-faces
  allocate(IFaceVector( NI,NJ-1,2)) ! Face Vectors for I-faces
  allocate(JFaceCenter( NI-1,NJ,2)) ! Face Centers for J-faces
  allocate(JFaceVector( NI-1,NJ,2)) ! Face Vectors for J-faces
  allocate(GradP(0:NI,0:NJ,2),GradP_t(0:NI,0:NJ,2),GradP_res(0:NI,0:NJ,2))  		! Pressure gradients array
  allocate(V(0:NI,0:NJ,2),DivV(0:NI,0:NJ),DivV_t(0:NI,0:NJ),DivV_res(0:NI,0:NJ), &
&			             DivVP(0:NI,0:NJ),DivVP_t(0:NI,0:NJ),DivVP_res(0:NI,0:NJ))	! Velocity vector and divergence arrays
  allocate(lapP(0:NI,0:NJ),lapP_t(0:NI,0:NJ),lapP_res(0:NI,0:NJ))
  allocate(CurlV(0:NI,0:NJ),CurlV_t(0:NI,0:NJ),CurlV_res(0:NI,0:NJ))
  allocate(T(0:NI,0:NJ),T_new(0:NI,0:NJ),gradT(0:NI,0:NJ,2),conv(0:NI,0:NJ),dif(0:NI,0:NJ),res(0:NI,0:NJ))
  gradP = 0
  divV = 0
  divVP = 0
  lapP = 0
  CurlV = 0
  T = 0
  gradT = 0
  conv = 0
  dif = 0
  

!===  READ GRID ===
  WRITE(*,*) 'Read mesh from file: ', MeshFile
  READ(IO,*) ((X(I,J),Y(I,J),I=1,NI),J=1,NJ)
  CLOSE(IO)

!=== CALCULATE METRIC ===
  WRITE(*,*) 'Calculate metric'       
  Call B_CalcMetric(NI,NJ,X,Y,CellCenter,CellVolume,IFaceCenter,IFaceVector,JFaceCenter,JFaceVector) 
  
!=== INITIATE FIELDS ===
  WRITE(*,*) 'Initiate fields'       
  DO  J = 0,NJ
    DO  I = 0,NI
      P(I,J) = Pressure(CellCenter(I,J,1),CellCenter(i,j,2))
	  V(I,J,:) = rVelocity(CellCenter(I,J,1),CellCenter(i,j,2))
	  GradP_t(I,J,:) = rGradP_ter(CellCenter(I,J,1),CellCenter(i,j,2))
	  divV_t(i,j) = rdivV_ter(CellCenter(I,J,1),CellCenter(i,j,2))
	  divVP_t(i,j) = rdivVP_ter(CellCenter(I,J,1),CellCenter(i,j,2))
	  lapP_t(i,j) = rlapP_ter(CellCenter(I,J,1),CellCenter(i,j,2))
	  curlV_t(i,j) = rCurlV_ter(CellCenter(I,J,1),CellCenter(i,j,2))
    ENDDO
  ENDDO

!=== INITIATE FIELDS ===
  if (sols == 1) then
	P = 0; V = 0 
	open(io,file=sol_file)
	read(io,*)
	read(io,*)
	read(io,*) ((rtmp,rtmp,V(i,j,1),V(i,j,2),rtmp,P(i,j),rtmp,rtmp, i=0,NI), J=0,NJ)
	close(io)
  end if

!=== CALCULATE GRADIENT ===
  WRITE(*,*) 'Calculate derivatives'
  do i=1,20
  Call B_CalcGradient(NI,NJ,X,Y,P,GradP,CellVolume,CellCenter,    	  &
&											IFaceVector,JFaceVector,  &
&											IFaceCenter,JFaceCenter)
  Call B_CalcGradRes(ni,nj,GradP,GradP_t,GradP_res)
  enddo
!===CALCULATE DIVERGENCE ===
  call B_CalcDiv(NI,NJ,X,Y,V,DivV,CellVolume,CellCenter,    &
&											IFaceVector,JFaceVector,  &
&											IFaceCenter,JFaceCenter)
  call B_CalcDivRes(NI,NJ,divV,divV_t,divV_res)

  call B_CalcDivphi(NI,NJ,X,Y,V,P,gradP,DivVP,CellVolume,CellCenter,    &
&											IFaceVector,JFaceVector,  &
&											IFaceCenter,JFaceCenter,sch)
  call B_CalcDivphiRes(NI,NJ,divVP,divVP_t,divVP_res)

!===CALCULATE LAPLACIAN ===
  call B_CalcLap(NI,NJ,X,Y,p,gradP,lapP,CellVolume,CellCenter,    &
&											IFaceVector,JFaceVector,  &
&											IFaceCenter,JFaceCenter)
  call B_CalcLapRes(NI,NJ,lapP,lapP_t,lapP_res)

!===CALCULATE CURL ===
  call B_CalcCurl(NI,NJ,X,Y,V,CurlV,CellVolume,CellCenter,    &
&											IFaceVector,JFaceVector,  &
&											IFaceCenter,JFaceCenter)
  call B_CalcCurlRes(NI,NJ,curlV,curlV_t,curlV_res)

! Решение уравнения конвективно-диффузионного переноса
dt = CFL * 0.02 / maxval(V)
!Граничные условия
T = 100
T(0,:) = T1
T(NI,:) = T2
T(:,0) = T(:,1)
T(:,NJ) = T(:,NJ-1)


  do m=1,nm
		!Расчёт невязки
		Call B_CalcGradient(NI,NJ,X,Y,T,GradT,CellVolume,CellCenter,    	  &
		&											IFaceVector,JFaceVector,  &
		&											IFaceCenter,JFaceCenter)
		call B_CalcDivphi(NI,NJ,X,Y,V,T,gradT,conv,CellVolume,CellCenter,    &
		&											IFaceVector,JFaceVector,   &
		&											IFaceCenter,JFaceCenter,sch)
		call B_CalcLap(NI,NJ,X,Y,T,gradT,dif,CellVolume,CellCenter,    &
		&											IFaceVector,JFaceVector,  &
		&											IFaceCenter,JFaceCenter)
		res = dif/(Re*Pr) - conv
		
		!Расчёт температуры на следующем слое
		T_new = T + res*dt
		
		T(:,0) = T(:,1)
		T(:,NJ) = T(:,NJ-1)
		!Проверка сходимости
		eps_T = maxval(abs(T_new(1:ni-1,1:nj-1) - T(1:ni-1,1:nj-1))) / maxval(abs(T_new(1:ni-1,1:nj-1)))
		print*, 'Iteration:', m, 'eps_T:', eps_T
		
		if (eps_T < eps) exit
		T = T_new
		
  end do

!=== OUTPUT FIELDS ===
  WRITE(*,*) 'Output fields to file: ', OutputFile       
  Open(IO,FILE=OutputFile)
  Call B_OutputFields(IO,NI,NJ,X,Y,P,V,GradP,GradP_res,divV,divV_res,lapP,lapP_res,divVP,divVP_res,curlV,curlV_res,T)
  Close(IO)


	contains
	
Function Pressure(X,Y)
  REAL:: PRESSURE, X, Y
  Pressure = x**4+y**4
End Function

Function rVelocity(X,Y)
  REAL:: rVelocity(2), X, Y
  rVelocity(1) = x*y**3
  rVelocity(2) = 2*x**3*y
End Function

Function rGradP_ter(x,y)
  real:: rGradP_ter(2), x, y
  rGradP_ter(1) = 3*x**2
  rGradP_ter(2) = 3*y**2
End Function

Function rdivV_ter(X,Y)
  REAL:: rdivV_ter, X, Y
  rdivV_ter = 2
End Function

Function rdivVP_ter(X,Y)
  REAL:: rdivVP_ter, X, Y
  rdivVP_ter = 2
End Function

Function rlapP_ter(X,Y)
  REAL:: rlapP_ter, X, Y
  rlapP_ter = 12*x**2+12*y**2
End Function

Function rCurlV_ter(X,Y)
  REAL:: rCurlV_ter, X, Y
  rCurlV_ter = 6*x**2*y-3*y**2*x
End Function

END PROGRAM Main  
