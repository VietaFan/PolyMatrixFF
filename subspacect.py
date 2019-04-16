from itertools import product
def getTables(p, r):
    if r != 1:
        raise Exception('not implemented yet')
    add = []
    mult = []
    for i in range(p):
        add.append([(i+j)%p for j in range(p)])
        mult.append([(i*j)%p for j in range(p)])
    return add, mult

def matadd(M1, M2, tables):
    add, mul = tables
    return [[add[M1[i][j]][M2[i][j]] for j in range(len(M1[0]))] for i in range(len(M1))]

def matmul(M1, M2, tables):
    add, mul = tables
    ans = []
    N = len(M1)
    for i in range(N):
        row = []
        for j in range(N):
            elem = 0
            for k in range(N):
                elem = add[elem][mul[M1[i][k]][M2[k][j]]]
            
            row.append(elem)
        ans.append(row)
    return ans

def eye(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]

def zeros(n):
    return [[0 for j in range(n)] for i in range(n)]

def commutes(M1, M2, tables):
    return matmul(M1, M2, tables) == matmul(M2, M1, tables)

def anticommutes(M1, M2, tables):
    return matadd(matmul(M1,M2,tables), matmul(M2,M1,tables), tables) == zeros(len(M1))

def comm(X, mats, tables):
    S = []
    for Y in mats:
        if commutes(X,Y,tables):
            S.append(Y)
    return S

def acomm(X, mats, tables):
    S = []
    for Y in mats:
        if anticommutes(X,Y,tables):
            S.append(Y)
    return S

# matrix kernel - kernel of X as an operator from GL(V) to GL(V) 
def matker(X, mats, tables):
    ker = []
    for Y in mats:
        if matmul(X,Y,tables) == matmul(Y,X,tables) == zeros(len(X)):
            ker.append(Y)
    return ker

def poly(X):

tables = getTables(3,1)
M_23 = [[[a,b], [c,d]] for (a,b,c,d) in product([0,1,2],[0,1,2],[0,1,2],[0,1,2])]
counts = []
for i in range(81):
    x = 0
    for j in range(81):
        if commutes(M_23[i],M_23[j],tables):
            x += 1
    counts.append(x)

acounts = []
for i in range(81):
    x = 0
    for j in range(81):
        if anticommutes(M_23[i],M_23[j],tables):
            x += 1
    acounts.append(x)

kernels = []
for i in range(81):
    kernels.append(len(matker(M_23[i],M_23,tables)))
print(counts)
print(acounts)
print(kernels)
print(acomm(eye(2), M_23, tables))
