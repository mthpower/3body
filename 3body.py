#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  3body.py
#  
#  Copyright 2012 Matthew <matthew@matthew-VirtualBox>
#-----------------------------------------------------------------------------
#imports

from physics import *
import pylab
import math


#-----------------------------------------------------------------------------
#global vars

N=5000
omega=0.0033
omsq=omega*omega
time=1000000

#-----------------------------------------------------------------------------
#coords and eqns

def X1(T):
    X1=math.cos(T)
    return X1

def Y1(T):
    Y1=math.sin(T)
    return Y1

def X2(T):
    X2=-1*math.cos(T)
    return X2

def Y2(T):
    Y2=-1*math.sin(T)
    return Y2

def R1(x,y,T):
    r1=math.sqrt(((x-X1(T))*(x-X1(T)))+((y-Y1(T))*(y-Y1(T))))
    return r1

def R2(x,y,T):
    r2=math.sqrt(((x-X2(T))*(x-X2(T)))+((y-Y2(T))*(y-Y2(T))))
    return r2

def I(x,y,T):
    I=x*math.cos(T)+y*math.sin(T)
    return I

def J(x,y,T):
    J=-1*x*math.sin(T)+y*math.cos(T)
    return J

#----------------------------------------------------------------------------

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
    dy[1]=-4*omsq*(((y[2]-X1(x))/((R1(y[2],y[4],x))**3))+((y[2]-X2(x))/((R2(y[2],y[4],x))**3)))
    dy[2]=y[1]
    dy[3]=-4*omsq*(((y[4]-Y1(x))/((R1(y[2],y[4],x))**3))+((y[4]-Y2(x))/((R2(y[2],y[4],x))**3)))
    dy[4]=y[3]
    return dy
    
#----------------------------------------------------------------------------




x=0.0; y=[0, 0.002, 10.0, -0.002, 10.0]

zx=[0 for j in range(0,N)]
zy=[0 for j in range(0,N)]
zxp=[0 for j in range(0,N)]
zyp=[0 for j in range(0,N)]
zX1=[0 for j in range(0,N)]
zY1=[0 for j in range(0,N)]
zX2=[0 for j in range(0,N)]
zY2=[0 for j in range(0,N)]


for j in range(0,N):
    (x,y) = runkut(4, x, y, time/N)
    #print x, y[2], y[4] 
    zxp[j]=y[2]
    zyp[j]=y[4] 
    zx[j]=I(y[2], y[4], x)
    zy[j]=J(y[2], y[4], x)
    zX1[j]=X1(x)
    zY1[j]=Y1(x)    
    zX2[j]=X2(x)
    zY2[j]=Y2(x)

pylab.plot(zx,zy)
pylab.plot(zxp,zyp)
pylab.plot(zX1,zY1)
pylab.plot(zX2,zY2)
pylab.show()
