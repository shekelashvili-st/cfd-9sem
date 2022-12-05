Function Pressure(X,Y)
  REAL:: PRESSURE, X, Y
  Pressure = x**2+y**2
End Function

Function rVelocity(X,Y)
  REAL:: rVelocity(2), X, Y
  rVelocity(1) = x
  rVelocity(2) = y
End Function

Function rGradP_ter(x,y)
  real:: rGradP_ter, x, y
  rGradP_ter = 1 
End Function

Function rdivV_ter(X,Y)
  REAL:: rdivV_ter, X, Y
  rdivV_ter = 2
End Function

Function rdivVP_ter(X,Y)
  REAL:: rdivVP_ter, X, Y
  rdivVP_ter = 4*(x*x+y*y)
End Function

Function rlapP_ter(X,Y)
  REAL:: rlapP_ter, X, Y
  rlapP_ter = 4
End Function