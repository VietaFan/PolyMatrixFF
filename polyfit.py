def polyfit(X,Y,d,maxcoeff=10,avoiding=set()):
    if d == 0:
        if len(set(Y)) == 1 and abs(Y[0]) < maxcoeff:
            return [Y[0]]
        else:
            return None
    for i in range(len(X)):
        if maxcoeff*sum([abs(X[i])**k for k in range(d+1)]) < abs(Y[i]):
            return None
    P = polyfit(X, Y, d-1, maxcoeff, avoiding)
    if P != None and tuple([0]+P) not in avoiding:
        return [0]+P
    for k in range(1, maxcoeff+1):
        P = polyfit(X, [y-k*x**d for x,y in zip(X,Y)], d-1, maxcoeff, avoiding)
        if P != None and tuple([k]+P) not in avoiding:
            return [k]+P
        P = polyfit(X, [y+k*x**d for x,y in zip(X,Y)], d-1, maxcoeff, avoiding)
        if P != None and tuple([-k]+P) not in avoiding:
            return [-k]+P

def getpoly(coeffs):
    return lambda q: sum([q**k*coeffs[-1-k] for k in range(len(coeffs))])
