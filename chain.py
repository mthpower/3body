from physics import *
import pylab

def runkut(n, x, y, h):
	"Advances the solution of diff eqn defined by derivs from x to x+h"
	y0=y[:]
	k1=derivs(n, x, y)
	for i in range(1,n+1): y[i]=y0[i]+0.5*h*k1[i]
	k2=derivs(n, x+0.5*h, y)
	for i in range(1,n+1): y[i]=y0[i]+h*(0.2071067811*k1[i]+0.2928932188*k2[i])
	k3=derivs(n, x+0.5*h, y)
	for i in range(1,n+1): y[i]=y0[i]-h*(0.7071067811*k2[i]-1.7071067811*k3[i])
	k4=derivs(n, x+h, y)
	for i in range(1,n+1):
		a=k1[i]+0.5857864376*k2[i]+3.4142135623*k3[i]+k4[i]
		y[i]=y0[i]+0.16666666667*h*a
	
	x+=h
	return (x,y)

#----------------------------------------------------------------------------

def derivs(n, x, y):
	"The function DERIVS calculates y' from x and y"
	dy=[0 for i in range(0,n+1)]
	if x==0.0: dy[1]=-omsq
	else:      dy[1]=y[2]/x

	dy[2]=-omsq*y[1]
	return dy
#----------------------------------------------------------------------------

N=40
omsq=1.

x=0.0; y=[0, 1.0, 0.0]             # Set Boundary Conditions



# Calculate and print the solution for x= 0 to 1
# h=0.025 is used this may not be an appropriate value

z=[0 for j in range(0,N)]
zeroline = [0 for j in range(0,N)]

for j in range(0,N):
	(x,y) = runkut(2, x, y, 1.0/N)
	print x, y[1], y[2]
	z[j]=y[1]

pylab.plot(zeroline)
pylab.plot(z)
pylab.show()
