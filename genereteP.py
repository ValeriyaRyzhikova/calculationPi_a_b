from sympy import *
x=symbols("x", integer=True, positive=True)
import numpy
from numpy import linalg as LA
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
        result=result*(i+1);
    return result

def L(k1, k2, x):
    c = factorial(N)/(factorial(k1)*factorial(k2))
    return c*(m1*m2*(1-r1))**(N-k1-k2)*(m1*(1-q))**k2*(m2*q)**k1*(l+x)**(k1+k2)


def returnR(x):
    R=[[0 for col in range(N+1)] for row in range(N+1)]
    sum=0
    n1=0
    while n1<=N:
        n2=0
        while n2<=N-n1:
            R[n1][n2]=L(n1,n2,x)
            sum=sum+R[n1][n2]
            n2=n2+1
        n1=n1+1
    
    n1=0
    while n1<=N:
        n2=0
        while n2<=N-n1:
            R[n1][n2]=R[n1][n2]/sum
            n2=n2+1
        n1=n1+1
    return R


def returnA(x):
    sum1 = 0
    sum2 = 0
    sum3 = 0
    sum4 = 0
    n1=0
    R=returnR(x)
    while n1<=N:
        sum1 =sum1+ R[n1][N-n1]
        n1=n1+1
    n1=0
    while n1<=N:
        n2=0
        while n2<N-n1:
            sum2 =sum2+ R[n1][n2]
            n2=n2+1
        n1=n1+1
    n1=1
    while n1<=N:
        n2=0
        while n2<=N-n1:
            sum3 =sum3 +n1*R[n1][n2]
            n2=n2+1
        n1=n1+1
    n1=0
    while n1<=N:
        n2=1
        while n2<=N-n1:
            sum4 += n2 * R[n1][n2]
            n2=n2+1
        n1=n1+1
    return l * sum1 - x * sum2 + m1 * r2 * sum3 + m2 * r2 * sum4


def returnLeft(x):
    a=[[0 for col in range(sN)] for row in range(sN)]

    n1=0
    n2=0
    a[round((n1+n2)*(n1+n2+1)/2)+n1][0]=-(l+x)
    a[round((n1+n2)*(n1+n2+1)/2)+n1][2]=m1*r0+m1*r2
    a[round((n1+n2)*(n1+n2+1)/2)+n1][1]=m2*r0+m2*r2

    n1=0
    n2=1
    while n2<N:
        a[round((n1+n2)*(n1+n2+1)/2)+n1][round((n2)*(n2+1)/2)]=-(l+m2*n2+x)+m2*r1*(1-q)*n2
        a[round((n1+n2)*(n1+n2+1)/2)+n1][round((n2)*(n2-1)/2)]=l*(1-q)+x*(1-q)
        a[round((n1+n2)*(n1+n2+1)/2)+n1][round((n2+1)*(n2+2)/2)+1]=m1*r0+m1*r2
        a[round((n1+n2)*(n1+n2+1)/2)+n1][round((n2+2)*(n2+1)/2)]=m2*r0*(n2+1)+m2*r2*(n2+1)
        a[round((n1+n2)*(n1+n2+1)/2)+n1][round((n2)*(n2+1)/2+1)]=m1*r1*(1-q)
        n2=n2+1

    n1=1
    n2=0
    while n1<N:
        a[round((n1+n2)*(n1+n2+1)/2)+n1][round((n1)*(n1+1)/2)+n1]=-(l+m1*n1+x)+m1*r1*q*n1
        a[round((n1+n2)*(n1+n2+1)/2)+n1][round((n2)*(n1-1)/2)+n1-1]=l*q+x*q
        a[round((n1+n2)*(n1+n2+1)/2)+n1][round((n1+1)*(n1+2)/2)+n1+1]=m1*r0*(n1+1)+m1*r2*(n1+1)
        a[round((n1+n2)*(n1+n2+1)/2)+n1][round((n1+2)*(n1+1)/2)+n1]=m2*r0+m2*r2
        a[round((n1+n2)*(n1+n2+1)/2)+n1][round((n1)*(n1+1)/2)+n1-1]=m2*r1*q
        n1=n1+1
    
    n1=1
    n2=1
    while n1<N:
        n2=1
        while n2<N-n1:
            a[round((n1+n2)*(n1+n2+1)/2)+n1][round((n1+n2)*(n1+n2+1)/2)+n1]=-(l+m1*n1+m2*n2+x)+m1*r1*q*n1+m2*r1*(1-q)*n2
            a[round((n1+n2)*(n1+n2+1)/2)+n1][round((n1+n2)*(n1+n2-1)/2)+n1-1]=l*q+x*q
            a[round((n1+n2)*(n1+n2+1)/2)+n1][round((n1+n2)*(n1+n2-1)/2)+n1]=l*(1-q)+x*(1-q)
            a[round((n1+n2)*(n1+n2+1)/2)+n1][round((n1+n2+2)*(n1+n2+1)/2)+n1+1]=m1*r0*(n1+1)+m1*r2*(n1+1)
            a[round((n1+n2)*(n1+n2+1)/2)+n1][round((n1+n2+2)*(n1+n2+1)/2)+n1]=m2*r0*(n2+1)+m2*r2*(n2+1)
            a[round((n1+n2)*(n1+n2+1)/2)+n1][round((n1+n2)*(n1+n2+1)/2)+n1+1]=m1*r1*(1-q)*(n1+1)
            a[round((n1+n2)*(n1+n2+1)/2)+n1][round((n1+n2)*(n1+n2+1)/2)+n1-1]=m2*r1*q*(n2+1)
            n2=n2+1
        n1=n1+1

    n1=0
    n2=N
    a[round((n1+n2)*(n1+n2+1)/2)+n1][round((N+1)*N/2)]= -N*m2+N*m2*r1*(1-q)
    a[round((n1+n2)*(n1+n2+1)/2)+n1][round((N-1)*N/2)]=l*(1-q)+x*(1-q)
    a[round((n1+n2)*(n1+n2+1)/2)+n1][round((N+1)*N/2)+1]=m1*r1*(1-q)

    n1=N
    n2=0
    a[round((n1+n2)*(n1+n2+1)/2)+n1][round((N+1)*N/2)+N]= -N*m1+N*m1*r1*q
    a[round((n1+n2)*(n1+n2+1)/2)+n1][round((N-1)*N/2)+N-1]=l*q+x*q
    a[round((n1+n2)*(n1+n2+1)/2)+n1][round((N+1)*N/2)+N-1]=m2*r1*q

    n1=1
    while n1<N:
        n2=N-n1
        a[round((n1+n2)*(n1+n2+1)/2)+n1][round((N)*(N+1)/2)+n1]=-(m1*n1+m2*n2)+m1*r1*q*n1+m2*r1*(1-q)*n2
        a[round((n1+n2)*(n1+n2+1)/2)+n1][round((N)*(N-1)/2)+n1-1]=l*q+x*q
        a[round((n1+n2)*(n1+n2+1)/2)+n1][round((N)*(N-1)/2)+n1]=l*(1-q)+x*(1-q)
        a[round((n1+n2)*(n1+n2+1)/2)+n1][round((N)*(N+1)/2)+n1+1]=m1*r1*(1-q)*(n1+1)
        a[round((n1+n2)*(n1+n2+1)/2)+n1][round((N)*(N+1)/2)+n1-1]=m2*r1*q*(n2+1)
        n1=n1+1
    
    for i in range(sN):
        a[sN-1][i]=1
    return a

def returnRight(x):
    R=returnR(x)
    ax=returnA(x)
    b=[0 for row in range(sN)]
    b[0]=R[0][0]*ax-m2*r2*R[0][1]-m1*r2*R[1][0]
    n1=0
    n2=1
    while n2<N:
        b[round((n1+n2)*(n1+n2+1)/2)+n1]=R[0][n2]*ax+x*(1-q)*R[0][n2-1]-m1*r2*R[1][n2]-m2*r2*(n2+1)*R[0][n2+1]
        n2=n2+1
    n1=1
    n2=0
    while n1<N:
        b[round((n1+n2)*(n1+n2+1)/2)+n1]=R[n1][0]*ax+x*q*R[n1-1][0]-m1*r2*(n1+1)*R[n1+1][0]-m2*r2*R[n1][1]
        n1=n1+1

    n1=1
    n2=1
    while n1<N:
        n2=1
        while n2<N-n1:
            b[round((n1+n2)*(n1+n2+1)/2)+n1]=R[n1][n2]*ax+x*q*R[n1-1][n2]+x*(1-q)*R[n1][n2-1]-m1*r2*(n1+1)*R[n1+1][n2]-m2*r2*(n2+1)*R[n1][n2+1]
            n2=n2+1
        n1=n1+1

    n1=0
    n2=N
    b[round((n1+n2)*(n1+n2+1)/2)+n1]=R[n1][n2]*ax+x*(1-q)*R[n1][n2-1]-l*R[n1][n2]

    n1=N
    n2=0
    b[round((n1+n2)*(n1+n2+1)/2)+n1]=R[n1][n2]*(ax-l)+x*q*R[n1-1][n2]

    n1=1
    while n1<N:
        n2=N-n1
        b[round((n1+n2)*(n1+n2+1)/2)+n1]=R[n1][n2]*ax+x*q*R[n1-1][n2]+x*(1-q)*R[n1][n2-1]-l*R[n1][n2]
        n1=n1+1
    b[sN-1]=0
    return b



def returnG(x):
    a=returnLeft(x)
    b=returnRight(x)
    g=LA.solve(a,b)
    return g

def returnG(x):
    a=returnLeft(x)
    b=returnRight(x)
    g = MatrixSymbol('g', sN, 1)
    result= solve(Matrix(a)*Matrix(g)-Matrix(b),Matrix(g),minimal=True)
    sResult=[0 for row in range(sN)]
    for i in range(sN):
        sResult[i]=result[g[i]]
    return sResult

def returnB(x):
    R=returnR(x)
    g=returnG(x)
    result=0
    for i in range(N+1):
       j=N-i
       result=result+2*l*g[round((i+j)*(i+j+1)/2)+i]

    for i in range(N+1):
        for j in range(N+1-i):
            result=result+2*m1*r2*i*g[round((i+j)*(i+j+1)/2)+i]

    for i in range(N+1):
        for j in range(N+1-i):
            result=result+2*m2*r2*j*g[round((i+j)*(i+j+1)/2)+i]

    for i in range(N):
        for j in range(N-i):
            result=result-2*x*g[round((i+j)*(i+j+1)/2)+i]

    for i in range(N):
        for j in range(N-i):
            result=result+2*x*R[i][j]
    
    ax=returnA(x)
    result=result+ax
    
    return result

def returnIntegralSimvol(x):
    a=returnA(x)
    b=returnB(x)
    integ=(a/b).evalf(3)
    integ=simplify(integ)
    integ=nsimplify(integ)
    result = integrate(integ, x)
    return result

def returnPi(z):
    integ=returnIntegralSimvol(x)
    integ=integ
    b=returnB(x)
    result=[(C/b.subs(x,i*sig)*exp(2/sig*integ.subs(x,i*sig))) for i in z]
    return result

def returnP(z, Pi):
    sum=0
    for i in z:
        sum = sum + Pi[round(i)]
    result=[Pi[round(i)]/sum for i in z]
    return result

Nx=30
z=numpy.linspace(0, Nx, Nx+1)
Pi=returnPi(z)

f = open('Pi10.txt', 'w')
for index in Pi:
    f.write(str(index) + '\n')
f.close()

newPi=[i.evalf(3) for i in Pi]
f = open('newPi10.txt', 'w')
for index in Pi:
    f.write(str(index) + '\n')
f.close()

yP=returnP(z, newPi)
from matplotlib import pyplot as plt
plt.plot(z,yP)
plt.show()

f = open('newPi2.txt')
for line in f:int(line)