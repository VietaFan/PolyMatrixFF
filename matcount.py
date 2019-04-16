from finfield import *
# count solutions in 2x2 matrices to [[X,Y],Z] = 0
# first count solutions to [X,Y] = M, for each M

#returns a dictionary mapping each M in F_q to the number of
#pairs (X,Y) of 2x2 matrices over F_q with [X,Y] = M.
def countLieSols2(q):
    counts = {}
    reprs = {}
    nz_mats = []
    for b in GF(q):
        for c in GF(q):
            for d in GF(q):
                if b != 0 or c != 0 or d != 0:
                    nz_mats.append(Matrix([[gf0(q),b],[c,d]]))
                nz_mats.append(Matrix([[gf1(q),b],[c,d]]))
    for X in nz_mats:
        for Y in nz_mats:
            comm = X*Y-Y*X
            if X.mat[0][0] == 0 and Y.mat[0][0] == 0:
                h = comm.tuplify()
                if h not in counts:
                    counts[h] = 0
                    reprs[h] = comm
                counts[h] += 1
                continue
            for x in GF(q):
                if x == 0:
                    continue
                h = (x*comm).tuplify()
                if h not in counts:
                    counts[h] = 0
                    reprs[h] = x*comm
                if X.mat[0][0] == 0 or Y.mat[0][0] == 0:
                    counts[h] += 1
                else:
                    counts[h] += q-1
    counts[zeros(2,q).tuplify()] += 2*q**4-1 # add in pairs where one term is 0
    return counts, reprs

def countMSols2(q):
    M_mats = []
    for a in GF(q):
        for b in GF(q):
            for c in GF(q):
                M_mats.append(Matrix([[a,b],[c,-a]]))
    Z_mats = []
    for a in GF(q):
        for b in GF(q):
            for c in GF(q):
                for d in GF(q):
                    Z_mats.append(Matrix([[a,b],[c,-a]]))
    solcts = {}
    reprs = {}
    zero = zeros(2,q)
    for M in M_mats:
        h = M.tuplify()
        reprs[h] = M
        solcts[h] = 0
        for Z in Z_mats:
            if M*Z-Z*M == zero:
                solcts[h] += 1
    return solcts, reprs

pairify = lambda solcts: sorted(map(lambda x: (x,solcts[x]),solcts),
                                key=lambda x: x[1], reverse=True)



'''
def polyInterp(x,y):
    if max(map(abs, y)) == 0:
        return []
    q = list(map(lambda a,b: b%a, x, y))
    d = list(map(lambda a,b: b//a, x, y))
    print(q,d)
    if len(set(q)) > 1:
        return None
    return [q[0]] + polyInterp(x,d)
'''

def C(X, q):
    ctr = 0
    for a in GF(q):
        for b in GF(q):
            for c in GF(q):
                for d in GF(q):
                    Y = Matrix([[a,b],[c,d]])
                    if X*Y == Y*X:
                        ctr += 1
    return ctr

