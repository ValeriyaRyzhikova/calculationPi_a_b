import numpy
from numpy import linalg as LA

def index(n1, n2):
    return int((n1+n2)*(n1+n2+1)/2+n1)

N = 2
sN=round((N+1)*(N+2)/2)
m1 = 0.6
m2 = 1.5
r1 = 0.2 
r2 = 0.1
r0 = 1 - r1 - r2
q = 0.333
l = 0.8
sig=0.2
C=1
def factorial(a):
    result =1
    for i in range(a):
        result*=(i+1);
    return result

def L(k1, k2, x):
    c = factorial(N)/(factorial(k1)*factorial(k2))
    return c*(m1*m2*(1-r1))**(N-k1-k2)*(m1*(1-q))**k2*(m2*q)**k1*(l+x)**(k1+k2)


def returnR(x):
    R=[[0 for col in range(N+1)] for row in range(N+1)]
    sum=0

    for n1 in range(N+1):
        for n2 in range(N+1-n1):
            R[n1][n2]=L(n1,n2,x)
            sum+=R[n1][n2]
    
    for n1 in range(N+1):
        for n2 in range(N+1-n1):
            R[n1][n2]=R[n1][n2]/sum

    return R


def returnA(x):
    sum1 = 0
    sum2 = 0
    sum3 = 0
    sum4 = 0
    R=returnR(x)
    for n1 in range(N+1):
        sum1 =sum1+ R[n1][N-n1]

    for n1 in range(N):
        for n2 in range(N-n1):
            sum2 =sum2+ R[n1][n2]

    for n1 in range(N+1):
        for n2 in range(N+1-n1):
            sum3 =sum3 +n1*R[n1][n2]

    for n1 in range(N+1):
        for n2 in range(N+1-n1):
            sum4 =sum4 + n2 * R[n1][n2]

    return l * sum1 - x * sum2 + m1 * r2 * sum3 + m2 * r2 * sum4

def E(a,b):
    if a==b:
        return 1
    else:
        return 0

def notE(a,b):
    if a==b:
        return 0
    else:
        return 1

def returnLeft(x):
    a=[[0 for col in range(sN)] for row in range(sN)]

    for n1 in range(N+1):
        for n2 in range(N+1-n1):
            a[index(n1,n2)][index(n1,n2)]=-(l+m1*n1+m2*n2+x)+m1*r1*q*n1+m2*r1*(1-q)*n2+l*E(N, n1+n2)-x*notE(N, n1+n1)
            if notE(0, n1)==1:
                a[index(n1,n2)][index(n1-1,n2)]=l*q+x*q
                a[index(n1,n2)][index(n1-1,n2+1)]=m2*r1*q*(n2+1)
            if notE(0, n2)==1:
                a[index(n1,n2)][index(n1,n2-1)]=l*(1-q)+x*(1-q)
                a[index(n1,n2)][index(n1+1,n2-1)]=m1*r1*(1-q)*(n1+1)
            if notE(N, n1+n2)==1:
                a[index(n1,n2)][index(n1+1,n2)]=m1*r0*(n1+1)+m1*r2*(n1+1)
                a[index(n1,n2)][index(n1,n2+1)]=m2*r0*(n2+1)+m2*r2*(n2+1)

    for i in range(sN):
        a[sN-1][i]=1
    return a

def returnRight(x):
    R=returnR(x)
    ax=returnA(x)
    b=[0 for row in range(sN)]

    for n1 in range(N+1):
        for n2 in range(N+1-n1):
            b[index(n1,n2)]=R[n1][n2]*ax
            if notE(0, n1)==1:
                 b[index(n1,n2)]+=x*q*R[n1-1][n2]
            if notE(0, n2)==1:
                b[index(n1,n2)]+=x*(1-q)*R[n1][n2-1]
            if notE(N, n1+n2)==1:
                 b[index(n1,n2)]-=m1*r2*(n1+1)*R[n1+1][n2]
                 b[index(n1,n2)]-=m2*r2*(n2+1)*R[n1][n2+1]
            if E(N, n1+n2)==1:
                b[index(n1,n2)]-=l*R[n1][n2]

    b[sN-1]=0
    return b



def returnG(x):
    a=returnLeft(x)
    b=returnRight(x)
    g=LA.solve(a,b)
    return g


def returnB(x):
    R=returnR(x)
    g=returnG(x)
    result=0
    for i in range(N+1):
        j=N-i
        result+=2*l*g[round((i+j)*(i+j+1)/2)+i]

    for i in range(1, N+1):
        for j in range(N+1-i):
            result+=2*m1*r2*i*g[round((i+j)*(i+j+1)/2)+i]

    for i in range(N):
        for j in range(1, N+1-i):
            result+=2*m2*r2*j*g[round((i+j)*(i+j+1)/2)+i]

    for i in range(N):
        for j in range(N-i):
            result-=2*x*g[round((i+j)*(i+j+1)/2)+i]

    for i in range(N):
        for j in range(N-i):
            result+=2*x*R[i][j]
    
    ax=returnA(x)
    result=result+ax
  
    return result


from sympy import *
x=symbols("x", integer=True, positive=True)
a=returnA(x)
b=returnB(x)

from matplotlib import pyplot as plt
Nx=30
xCh=numpy.linspace(0, 10, Nx)
yA= [returnA(i) for i in xCh]
yB= [returnB(i) for i in xCh]

z=numpy.linspace(0, Nx, Nx+1)

newPi=[i.evalf(3) for i in Pi]
yP=returnP(z, newPi)
plt.plot(z,yP)


plt.plot(xCh,yA,label='a(x)')
plt.plot(xCh,yB,label='b(x)')
plt.xlable("x")
plt.ylable("a, b")
plt.legend()
plt.grid(True)
plt.show()


f = open('Pi2.txt', 'w')
for index in Pi:
    f.write(str(index) + '\n')
f.close()
vR= [0 for i in range(sN)]
for i in range(N+1):
    for j in range(N+1-i):
        vR[int((i + j) * (i + j + 1) / 2 + i)] = R[i][j]