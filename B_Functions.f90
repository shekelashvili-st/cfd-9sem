Function Pressure(X,Y)
  REAL:: PRESSURE, X, Y
  Pressure = x**2+y**2
End Function

Function Velocity(X,Y)
  REAL:: Velocity(2), X, Y
  Velocity(1) = x
  Velocity(2) = y
End Function

Function GradP_ter(x,y)
  real:: GradP_ter, x, y
  GradP_ter = 1 
End Function

Function divV_ter(X,Y)
  REAL:: divV_ter, X, Y
  divV_ter = 2
End Function

Function lapP_ter(X,Y)
  REAL:: lapP_ter, X, Y
  lapP_ter = 4
End Function