#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Tkinter import *
from Tkinter import _tkinter
from numpy import *
from numpy.linalg import eigh as _eigh, solve as _solve
from numpy.fft import fft as _fft, ifft as _ifft
from FileDialog import *
import random

def _maxmin(y,n):
    y1=1.0e30;  y2=-y1
    for i in range(0,n):
        if type(y[i]) is not NoneType:
            if y[i] < y1: y1=y[i]
            if y[i] > y2: y2=y[i]
    
    if y1 > 0.0: y1=0.0
    if y2 < 0.0: y2=0.0
    y1,y2 = (1.05*y1-0.05*y2), (1.05*y2-0.05*y1)
    return (y1,y2)

#---------------------------------------------------------------------------

def _layout(z1, z2):
    j=floor(log10(z2-z1))
    y=pow(10, j)
    if fabs(z2) > fabs(z1): m2=floor(log10(fabs(z2)))
    else:                   m2=floor(log10(fabs(z1)))
    if (m2<0): m2=0
    a=(z2-z1)/y
    dz=.2
    if a>1.5: dz=0.5
    if a>3.5: dz=1.0
    if a>7.0: dz=2.0
    m1 = int((dz < 1.0))    
    m1=m1+fabs(j)
    dz=dz*y
    if dz<1.0:
        d=8*(m1+m2+3)
        s = "%" + str(int(m1+m2+3)) + "." + str(int(m1)) + "f"
    else:
        d=8*(m2+1)
        s = "%" + str(int(m2+1)) + "g"
    
    n1=-int(floor(-(z1/dz-1.0E-8)))
    n2= int(floor(  z2/dz+1.0E-8))
    return (dz, n1, n2, d, s)

#---------------------------------------------------------------------------

def _axes(gx, gy, gh, gw, x1, x2, y1, y2):
    if y1==y2:
        y1 = 0.5*y1
        y2 = 2.0*y2
    if y2==0.0: y2 = 1.0
    
    gp = float(gw)/(x2-x1)
    gq = -gp*x1+gx
    gr = -float(gh)/(y2-y1)
    gs = -gr*y1+gy
    
    _line(gp*x1+gq-16, gs, gp*x1+gq+gw+16, gs, "black")
    _line(gq, gr*y1+gs+16, gq, gr*y1+gs-gh-16, "black")
    
    (dz, n1, n2, d, s) = _layout(x1,x2)
    for i in range(n1,n2+1):
        j = gp*i*dz + gq
        if i != 0:
            _line(j, gs, j, gs+8, "black")
            _graphtext(j-d/2, gs+12, s % (i*dz), "black")
    
    (dz, n1, n2, d, s) = _layout(y1,y2)
    for i in range(n1,n2+1):
        j = gr*i*dz + gs
        if i != 0:
            _line(gq, j, gq-8, j, "black")
            _graphtext(gq-d-16, j-4, s % (i*dz), "black")
    
    _graphtext(gq-16, gs+8, "0", "black")
    return (gp, gq, gr, gs)

#-----------------------------------------------------------------------------

def _line(x1, y1, x2, y2, col):
    _c.create_line(int(x1), int(y1), int(x2), int(y2), fill=col, tag="graphical")

#-----------------------------------------------------------------------------

def _bar(x1, y1, x2, y2, col):
    _c.create_rectangle(int(x1), int(y1), int(x2), int(y2), width=0, fill=col, tag="graphical")

#-----------------------------------------------------------------------------

def _graphtext(x, y, t, col):
    _c.create_text(int(x), int(y), text=t, fill=col, anchor="nw", tag="graphical", font="helvetica 10")

#-----------------------------------------------------------------------------

def _testgraph():
    global _graph_number, _c
    try:
        _ccopy=_c
    except NameError:
# _c does not exist so create it
        _c=Canvas(width=640, height=480, background="grey60")
        _c.pack()
        _c.master.title("Python Graph")
        _graph_number=0
    else:
# _c does exist so test if widget exists    
        try:
            test=Misc.winfo_exists(_c)
        except _tkinter.TclError:
# widget has been destroyed so delete _c and start again
            del(_c)
            _c=Canvas(width=640, height=480, background="grey60")
            _c.pack()
            _c.master.title("Python Graph")
            _graph_number=0
        else:
# widget exists so superimpose on existing graph
            _graph_number+=1

#------------------------------------------------------------------------

def graph(y):
# Plots list or array y.  Some list elements may have value None. The 
# corresponding points are not plotted, but no error is incurred.
    global _axes_params
    _testgraph()
    n=len(y)
    y1,y2=_maxmin(y, n)
    if _graph_number==0:
        _axes_params = _axes(80, 440, 400, 500, 0, n, y1, y2)
    graph_colour=["red", "yellow", "green", "blue", "violet"]
    if _graph_number<5: colour=graph_colour[_graph_number]
    else: colour="white"
    
    gp, gq, gr, gs = _axes_params
    for i in range(0,n-1):
        if (type(y[i]) is not NoneType) and (type(y[i+1]) is not NoneType):
            _line(gp*i+gq, gr*y[i]+gs, gp*(i+1)+gq, gr*y[i+1]+gs, colour)

#------------------------------------------------------------------------

def histogram(y):
    global _axes_params
    _testgraph()
    graph_colour=["red", "yellow", "green", "blue"]
    if _graph_number<4: colour=graph_colour[_graph_number]
    else: colour="white"
    n = len(y)
    y1, y2 = _maxmin(y, n)
    
    if _graph_number==0: _axes_params = _axes(80, 200, 180,200,0,n,y1,y2)
    if _graph_number==1: _axes_params = _axes(370,200, 180,200,0,n,y1,y2)
    if _graph_number==2: _axes_params = _axes(80, 440, 180,200,0,n,y1,y2)
    if _graph_number>=3: _axes_params = _axes(370,440, 180,200,0,n,y1,y2)
    gp, gq, gr, gs = _axes_params
    for i in range(0,n-1):
        if type(y[i]) is not NoneType:
            _bar(gp*i+gq, gr*y[i]+gs, gp*i+gq+gp, gs,colour)

#-----------------------------------------------------------------------------

def clear_graph():
    global _graph_number
    try: _c.delete("graphical")
    except _tkinter.TclError: pass
    except NameError: pass  
    
    _graph_number=-1

#-----------------------------------------------------------------------------

def axes(x1, x2, y1, y2):
    global _axes_params
    _testgraph()
    clear_graph()
    _axes_params=_axes(80, 440, 400, 500, x1, x2, y1, y2)

#-----------------------------------------------------------------------------

def point(x,y,colour_number):
    gp, gq, gr, gs = _axes_params
    graph_colour=["red", "yellow", "green", "blue", "violet"]
    if colour_number<5: colour=graph_colour[colour_number-1]
    else: colour="white"
    
    _graphtext(gp*x+gq-3, gr*y+gs-6, "+",colour)

#-----------------------------------------------------------------------------

def line(x1,y1,x2,y2,colour_number):
    gp, gq, gr, gs = _axes_params
    graph_colour=["red", "yellow", "green", "blue", "violet"]
    if colour_number<5: colour=graph_colour[colour_number-1]
    else: colour="white"
    
    _line(gp*x1+gq, gr*y1+gs, gp*x2+gq, gr*y2+gs, colour)

#-----------------------------------------------------------------------------

def contour(a):
    global _graph_number, _axes_params
    colours= ["violet", "blue2", "DeepSkyBlue1", "cyan2", "SeaGreen1",\
              "green2", "GreenYellow", "yellow2", "orange1", "red"]
    n=len(a)
    gx, gy, gh, gw = 80, 440, 400, 400
    _graph_number=0
    _testgraph()
    x1, x2, y1, y2 = 0.0, n-1, 0.0, n-1
    gp=gw/(x2-x1)
    gq=-gp*x1+gx
    gr=-gh/(y2-y1)
    gs=-gr*y1+gy
    _axes_params = (gp, gq, gr, gs)
    _line(gx,    gy,    gx+gw, gy,    "white")
    _line(gx+gw, gy,    gx+gw, gy-gh, "white")
    _line(gx+gw, gy-gh, gx,    gy-gh, "white")
    _line(gx,    gy-gh, gx,    gy,    "white")
    max=-1.0e30
    min= 1.0e30
    for i in range(0,n):
        for j in range(0,n):
            if a[i,j]>max: max=a[i,j]
            if a[i,j]<min: min=a[i,j]
    
    for k in range(0,10):
        l=min+(k+1)*(max-min)/11.0
        for i in range(0,n-1):
            for j in range(0,n-1):
                a0=0.25*(a[i,j]+a[i+1,j]+a[i,j+1]+a[i+1,j+1])
                _triangle(i, j, i+1.0, j, i+0.5, j+0.5,\
                          a[i,j], a[i+1,j], a0, l, colours[k])
                _triangle(i+1.0, j, i+1.0, j+1.0, i+0.5, j+0.5,\
                          a[i+1,j], a[i+1,j+1], a0, l, colours[k])
                _triangle(i+1.0, j+1.0, i, j+1.0, i+0.5, j+0.5,\
                          a[i+1,j+1], a[i,j+1], a0, l, colours[k])
                _triangle(i,j+1.0 , i, j, i+0.5, j+0.5,
                          a[i,j+1], a[i,j], a0, l, colours[k])

def _triangle(x1, y1, x2, y2, x3, y3, b1, b2, b3, l, col):
    x=[0,x1,x2,x3]
    y=[0,y1,y2,y3]
    b=[0,b1,b2,b3]
    s1=_side(1,2,l,b)
    s2=_side(2,3,l,b)
    s3=_side(3,1,l,b)
    if (s1 and s2): _join(1, 2, 2 ,3, x, y, l, b, col)
    if (s2 and s3): _join(2, 3, 3, 1, x, y, l, b, col)
    if (s3 and s1): _join(3, 1, 1, 2, x, y, l, b, col)

def _side(m, n, l, b):
    return not ( (l>=b[m]) ^ (l<b[n]) )

def _join(n1, n2, m1, m2, x, y, l, b, col):
    gp, gq, gr, gs = _axes_params
    x1=_cross( x[n1], x[n2], n1, n2, l, b)
    y1=_cross( y[n1], y[n2], n1, n2, l, b) 
    x2=_cross( x[m1], x[m2], m1, m2, l, b)
    y2=_cross( y[m1], y[m2], m1, m2, l, b)
    _line(gp*x1+gq, gr*y1+gs,gp*x2+gq, gr*y2+gs, col)

def _cross(u, v, n, m, l, b):
    return u+(l- b[n]) * (v-u)/(b[m] - b[n])

#----------------------------------------------------------------------------

def postscript():
    f=Frame()
    dialog=SaveFileDialog(f,)
    filename=dialog.go(pattern="*.ps")
    _c.postscript(file=filename, pagewidth="10.0c", pageheight="10.0c")

#----------------------------------------------------------------------------

def uniform(n):
    global _first_random
    if n==0: random.seed(1234); _first_random=0
    if _first_random==1: random.seed(None); _first_random=0
    if n==1: return random.random()
    else: return floor(n*random.random()+1)

#-----------------------------------------------------------------------------    
_first_random=1
#-----------------------------------------------------------------------------

def gaussian(x):
    return random.gauss(0.0,x)

#--------------------------------------------------------------------------    

def poisson(x):
    g=exp(-x)
    em=0
    t=random.random()
    while t>g: t*=random.random(); em+=1
    
    return em

#-----------------------------------------------------------------------------

def capitalj(n,x):
    even=True
    if (x==0.0) and (n==0): return 1.0
    if (x==0.0) and (n> 0): return 0.0
    nn=2*floor(.5*(abs(x)+n+20.0))
    a=0.0; b=1.0E-30; y=0.0; z=1.0
    i=nn
    while i>=0:
        if i==n: z=b
        c=2.0*i*b/x-a; a=b; b=c;
        if abs(c)>1.0e10: a*=1.0e-20; b*=1.0e-20; y*=1.0e-20; z*=1.0e-20
        if even>0: y=y+a
        even= not even
        i-=1
    
    y=y+y-a
    return z/y

#----------------------------------------------------------------------------

def littlej(n,x):
    if (x==0.0) and (n==0): return 1.0
    if (x==0.0) and (n> 0): return 0.0
    nn=2*floor(.5*(abs(x)+n+20.0))
    a=z=0.0;  b=1.0E-30; z=1.0
    i=nn
    while i>0:
        c=(2*i+1)*b/x-a; a=b; b=c
        if abs(c)>1.0e10: a*=1.0e-20; b*=1.0e-20; z*=1.0e-20
        if i==(n+1): z=b
        i-=1
    
    return z*sin(x+1.0e-10)/(c*(x+1.0e-10))

#-----------------------------------------------------------------------------

def erfc(x):
    y=1.0+x*(.0705230784+x*(.0422820123+x*(.0092705272+x*(1.520143e-4+x*\
    (2.765672e-4+x*4.30638e-5)))))
    y=y*y; y=y*y; y=y*y; y=y*y
    return 1.0/y

#-----------------------------------------------------------------------------


def fftmy(n, x, y, sn):
    z=[x[i]+y[i]*1j for i in range(0,n)]
    z = array(z)
    if sn>0.0: z=_ifft(z)*n
    else:      z=_fft(z)
    xx = zeros(n+1,float); yy = zeros(n+1,float)
    for i in range(0,n):
        xx[i]=z[i].real; yy[i]=z[i].imag
    return (xx,yy)


def fft(n, x, y, sn):
    z=[x[i]+y[i]*1j for i in range(0,n)]
    z = array(z)
    if sn>0.0: z=_ifft(z)*n
    else:      z=_fft(z)
    xx = zeros(n+1,float); yy = zeros(n+1,float)
    for i in range(0,n):
        xx[i]=z[i].real; yy[i]=z[i].imag
    return (xx,yy)

def fft_BUG(n, x, y, sn):
    z=[x[i]+y[i]*1j for i in range(1,n+1)]
    z = array(z)
    if sn>0.0: z=_ifft(z)*n
    else:      z=_fft(z)
    for i in range(1,n+1):
        x[i]=z[i-1].real; y[i]=z[i-1].imag
    return (x,y)


#----------------------------------------------------------------------------

def mateig(A):
    (d,A)=_eigh(A)
    return (A,d)

#-----------------------------------------------------------------------------

def matdiv(n,A,b):
    return _solve(A,b)

def pause():
    junk=raw_input("Press return")

# Make sure a new import clears any old graphs
clear_graph()

#Wm.wm_protocol(_c.master,"WM_DELETE_WINDOW",_c.destroy)
#obey Wm.wm_protocol(_c.master,"WM_DELETE_WINDOW")

#import types, defines NoneType ListType, Numeric defines ArrayType.
# use "is" comparison to test type().

# can use repeat to copy arrays b=repeat(a,1), better use b=array(a)!!!!!
# better still b=a[:]

