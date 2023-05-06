import numpy as np
import scipy.sparse as sp
import scipy.linalg as sl

inf_norm = []
l2_norm = []

for n in [2, 3, 4, 5, 6, 7, 8]:
    h = 2**(-n)
    M = [] # matrix
    #range over euqations (rows)
    #We'll solve M*u=v(f,g), where u is ordered as (1,1), (1,2)...,(1,n),(2,1),...
    for i in range(1,n+2):
            for j in range(1,n+2):
                #boundary
                if i==1 or j==1 or i==(n+1) or j==(n+1):
                    row = [0 for t in range((n+1)**2)]
                    row[(i-1)*(n+1)+(j-1)]= 1
                    M.append(row)
                else:
                    #internal
                    row = [0 for t in range((n+1)**2)]
                    row[(i-1)*(n+1)+(j-1)]= 4/h
                    row[(i-2)*(n+1)+(j-1)] = -1/h
                    row[(i-1)*(n+1)+(j-2)] = -1/h
                    row[(i)*(n+1)+(j-1)] = -1/h
                    row[(i-1)*(n+1)+(j)] = -1/h
                    M.append(row)
    Mat = sp.csr_array(M) #sparse storage
    #v(f,g)
    v = []
    def g(x,y):
        return (x**4)*(y**5)-(17*np.sin(x*y))
    def f(x,y):
        return (-40)*(x**3)*(y**4)-(17*x*y*np.sin(x*y))
    for i in range(1,n+2):
            for j in range(1,n+2):
                if i==1 or j==1 or i==(n+1) or j==(n+1):
                    v.append([g(i*h,j*h)])
                else:
                    #internal
                    v.append([f(i*h,j*h)])
    Vec = sp.csr_array(v) #sparse storage
    Sol = sl.solve(M,v) #figure out solve with sparse format...
    # print(Sol)

    def u(x, y):
        return (x ** 4) * (y ** 5) - (17 * np.sin(x * y))

    U = []
    for i in range(1,n+2):
            for j in range(1,n+2):
                if i==1 or j==1 or i==(n+1) or j==(n+1):
                    U.append([u(i*h,j*h)])
                else:
                    #internal
                    U.append([u(i*h,j*h)])

    max_u = 0
    for i in range(2,n+1):
            for j in range(2,n+1):
                tmp = np.abs(U[i*(n+1)+j][0] - Sol[i*(n+1)+j][0])
                if tmp > max_u:
                    max_u = tmp
        
    print(f"Infinite Norm Error for {n} = {max_u}")
    inf_norm.append(max_u)

# L2 Norm normalised by number of gridpoints

    max_u_sum = 0
    for i in range(2,n+1):
            for j in range(2,n+1):
                tmp_sum = (np.abs(U[i*(n+1)+j][0] - Sol[i*(n+1)+j][0])) ** 2

                max_u_sum = max_u_sum + tmp_sum
        
    print(f"L2 Norm Error for {n} = {np.sqrt(max_u_sum)/((n - 1)**2)}")
    l2_norm.append(np.sqrt(max_u_sum)/((n - 1)**2))

inf_norm
l2_norm