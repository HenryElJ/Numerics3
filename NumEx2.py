#Numerics 3 sheet 4 exercise 2

import numpy as np
import scipy.sparse as sp
import scipy.linalg as sl

inf_norm = []
l2_norm_list = []

def g(x,y):
    return (x**4)*(y**5)-(17*np.sin(x*y))
def f(x,y):
    return (-40)*(x**3)*(y**4)-(17*x*y*np.sin(x*y))
def u(x, y):
    return (x ** 4) * (y ** 5) - (17 * np.sin(x * y))

for n in [2, 3, 4, 5, 6, 7, 8]:
    h = 2**(-n)
    N = 2**n
    M = [] # matrix representing the equations of the stencil
    #range over euqations (rows)
    #We'll solve M*u=v(f,g), where u is ordered as (1,1), (1,2)...,(1,n),(2,1),...
    for i in range(1,N+2):
            for j in range(1,N+2):
                #boundary
                if i==1 or j==1 or i==(N+1) or j==(N+1):
                    row = [0 for t in range((N+1)**2)]
                    row[(i-1)*(N+1)+(j-1)]= 1
                    M.append(row)
                else:
                    #internal
                    row = [0 for t in range((N+1)**2)]
                    row[(i-1)*(N+1)+(j-1)]= 10/(3*h)
                    row[(i-2)*(N+1)+(j-1)] = -2/(3*h)
                    row[(i-1)*(N+1)+(j-2)] = -2/(3*h)
                    row[(i)*(N+1)+(j-1)] = -2/(3*h)
                    row[(i-1)*(N+1)+(j)] = -2/(3*h)
                    row[(i-2)*(N+1)+(j-2)] = -1/(6*h)
                    row[(i-2)*(N+1)+(j)] = -1/(6*h)
                    row[(i)*(N+1)+(j-2)] = -1/(6*h)
                    row[(i)*(N+1)+(j)] = -1/(6*h)
                    M.append(row)
    #Mat = sp.csr_array(M) #sparse storage
    
    v = []
    for i in range(1,N+2):
            for j in range(1,N+2):
                if i==1 or j==1 or i==(N+1) or j==(N+1):
                    v.append([g(i*h,j*h)])
                else:
                    #internal
                    v.append([f(i*h,j*h)])
    #Vec = sp.csr_array(v) #sparse storage
    Sol = sl.solve(M,v) #I don't know why but can't access scipy.sparse.linalg.spsolve ...
    
    U = []
    for i in range(1,N+2):
            for j in range(1,N+2):
                if i==1 or j==1 or i==(N+1) or j==(N+1):
                    U.append([u(i*h,j*h)])
                else:
                    #internal
                    U.append([u(i*h,j*h)])

    max_norm = 0
    for i in range(2,N+1):
            for j in range(2,N+1):
                tmp = np.abs(U[i*(N+1)+j][0] - Sol[i*(N+1)+j][0])
                if tmp > max_norm:
                    max_norm = tmp
        
    print(f"Infinite Norm Error for {n} = {max_norm}")
    inf_norm.append(max_norm)

    # L2 Norm normalised by number of gridpoints
    l2_norm = 0
    for i in range(2,N+1):
            for j in range(2,N+1):
                tmp= (np.abs(U[i*(N+1)+j][0] - Sol[i*(N+1)+j][0])) ** 2
                l2_norm = l2_norm + tmp
    
    l2_norm = np.sqrt(l2_norm)/((N - 1)**2)
    print(f"L2 Norm Error for {n} = {l2_norm}")
    l2_norm_list.append(l2_norm)
