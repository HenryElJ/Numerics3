import numpy as np
import scipy.sparse as sp
import scipy.linalg as sl

n = 2
h = 2**(-n)
M = [] # matrix
#range over euqations (rows)
#We'll solve M*u=v(f,g), where u is ordered as (1,1), (1,2)...,(1,n),(2,1),...
for i in range(1,n+1):
        for j in range(1,n+1):
            #boundary
            if i==1 or j==1 or i==n or j==n:
                row = [0 for t in range(n**2)]
                row[(i-1)*n+(j-1)]= 1
                M.append(row)
            else:
                #internal
                row = [0 for t in range(n**2)]
                row[(i-1)*n+(j-1)]= 4/h
                row[(i-2)*n+(j-1)] = -1/h
                row[(i-1)*n+(j-2)] = -1/h
                row[(i)*n+(j-1)] = -1/h
                row[(i-1)*n+(j)] = -1/h
                M.append(row)
Mat = sp.csr_array(M) #sparse storage
#v(f,g)
v = []
def g(x,y):
     return (x**4)*(y**5)-(17*np.sin(x*y))
def f(x,y):
     return (-40)*(x**3)*(y**4)-(17*x*y*np.sin(x*y))
for i in range(1,n+1):
        for j in range(1,n+1):
            if i==1 or j==1 or i==n or j==n:
                v.append([g(i*h,j*h)])
            else:
                #internal
                v.append([f(i*h,j*h)])
Vec = sp.csr_array(v) #sparse storage
Sol = sl.solve(M,v) #fingure out solve with sparse format...
print(Sol)






    

        