# code for implementing GF(p^n) courtesy of
# https://stackoverflow.com/questions/48065360/interpolate-polynomial-over-a-finite-field/48067397#48067397


import itertools
from functools import reduce
from sympy import symbols, Dummy
from sympy.polys.domains import ZZ
from sympy.polys.galoistools import (gf_irreducible_p, gf_add, \
                                     gf_sub, gf_mul, gf_rem, gf_gcdex)
from sympy.ntheory.primetest import isprime

class GF():
    def __init__(self, p, n=None):
        if n is None:
            L = primeFact(p)
            if len(set(L)) > 1:
                raise ValueError('order must be a prime power')
            p = L[0]
            n = len(L)
        p, n = int(p), int(n)
        if not isprime(p):
            raise ValueError("p must be a prime number, not %s" % p)
        if n <= 0:
            raise ValueError("n must be a positive integer, not %s" % n)
        self.p = p
        self.n = n
        if n == 1:
            self.reducing = [1, 0]
        else:
            for c in itertools.product(range(p), repeat=n):
              poly = (1, *c)
              if gf_irreducible_p(poly, p, ZZ):
                  self.reducing = poly
                  break
        self.iterpos = 0

    def add(self, x, y):
        return gf_add(x, y, self.p, ZZ)

    def sub(self, x, y):
        return gf_sub(x, y, self.p, ZZ)

    def mul(self, x, y):
        return gf_rem(gf_mul(x, y, self.p, ZZ), self.reducing, self.p, ZZ)

    def inv(self, x):
        s, t, h = gf_gcdex(x, self.reducing, self.p, ZZ)
        return s

    def eval_poly(self, poly, point):
        val = []
        for c in poly:
            val = self.mul(val, point)
            val = self.add(val, c)
        return val

    def __iter__(self):
        return self

    def __next__(self):
        if self.iterpos == self.p**self.n:
            raise StopIteration
        position = []
        x = self.iterpos
        for i in range(self.n):
            position.append(x%self.p)
            x //= self.p
        self.iterpos += 1
        return GF_elt(position[::-1], self)


def polyify(L):
    def get_poly_term(c,d):
        if c == 0:
            return ''
        elif d == 0:
            return str(c)
        elif d == 1:
            if c == 1:
                return 'x'
            elif c == -1:
                return '-x'
            return str(c)+'x'
        else:
            if c == 1:
                return 'x^'+str(d)
            elif c == -1:
                return '-x^'+str(d)
            return str(c)+'x^'+str(d)
    n = len(L)-1
    pairs = map(lambda x: (L[x],n-x), range(n+1))
    reps = map(lambda x: get_poly_term(*x), pairs)
    nzreps = list(filter(len, reps))
    terms = [nzreps[0]]+list(map(lambda x: x if x[0] == '-' else '+'+x, nzreps[1:]))
    return ''.join(terms)

# Scott Aaronson
# Quantum Computing since Democritus

class GF_elt():
    def __init__(self, coeffs, field, rep='poly', repvar='x'):
        if type(coeffs) == int:
            coeffs = [coeffs]
        else:
            n = 0
            while n < len(coeffs) and coeffs[n] == 0:
                n += 1
            coeffs = coeffs[n:]
        self.field = field
        self.rep = rep
        self.repvar = repvar
        self.coeffs = [((x % field.p)+field.p)%field.p for x in coeffs] # the seemingly unnecessary extra step ensures nonnegative coefficients
    def __add__(self, other):
        if type(other) == type(self):
            other = other.coeffs
        elif type(other) == int:
            other = [other]
        return GF_elt(self.field.add(self.coeffs, other), self.field)
    def __sub__(self, other):
        return self + (other*(-1))
    def __radd__(self, other):
        return self + other
    def __rmul__(self, other):
        return self * other
    def __rsub__(self, other):
        return (-1)*(self-other)
    def __mul__(self, other):
        if type(other) == int:
            return GF_elt([x*other for x in self.coeffs], self.field)
        elif type(other) == type(self):
            other = other.coeffs
        elif type(other) == Matrix:
            return other*self
        return GF_elt(self.field.mul(self.coeffs, other), self.field)
    def __truediv__(self, other):
        if sum(self.coeffs) == 0:
            raise ZeroDivisionError('attempted division by 0 in GF(%s^%s)' % (self.field.p,self.field.n))
        elif type(other) == type(self):
            other = other.coeffs
        return self * GF_elt(self.field.inv(other), self.field)
    def __repr__(self):
        if self.rep == 'tuple':
            if len(self.coeffs) == 0:
                return '0'
            return ','.join(map(str, self.coeffs))
        elif self.rep == 'poly':
            if len(self.coeffs) == 0:
                return '0'
            return polyify(self.coeffs)
    def __pow__(self, n):
        if n < 0:
            return GF_elt(self.field.inv(self.coeffs), self.field) ** (-n)
        if n == 0:
            return GF_elt([1], self.field)
        if n > 0:
            return self*(self**(n-1))
    def tuplify(self):
        return tuple(self.coeffs)
    def __eq__(self,x):
        if type(x) == GF_elt:
            return self.tuplify() == x.tuplify()
        if type(x) == int:
            if len(self.coeffs) > 1:
                return False
            if len(self.coeffs) == 0:
                return x == 0
            return x == self.coeffs[0]
    def __neg__(self):
        return 0-self
def gf0(p,r=None):
    L = primeFact(p)
    p, r = L[0], len(L)
    return GF_elt([], GF(p,r))
def gf1(p,r=None):
    L = primeFact(p)
    p, r = L[0], len(L)
    return GF_elt([1], GF(p,r))
def gfx(p,r=None):
    L = primeFact(p)
    p, r = L[0], len(L)
    return GF_elt([1,0], GF(p,r))

def primeFact(n):
    primeList = []
    while n > 1:
        going = 0
        for x in range(2, int(n**0.5)+1):
            while not n%x:
                n //= x
                primeList.append(x)
                going = 1
        if not going:
            primeList.append(n)
            return primeList
    return primeList

class Matrix():
    def __init__(self, A):
        self.mat = A
        self.d = len(A)
        self.p = A[0][0].field.p
        self.r = A[0][0].field.n
    def __rmul__(self, B):
        return self*B
    def __mul__(self, B):
        if type(B) in [int, GF_elt]:
            return Matrix([[B*x for x in self.mat[i]] for i in range(self.d)])
        new_mat = []
        if B.d != self.d:
            raise Exception('Error: matrices of different dimensions')
        for i in range(self.d):
            row = []
            for j in range(self.d):
                x = gf0(self.p,self.r)
                for k in range(self.d):
                    x += self.mat[i][k]*B.mat[k][j]
                row.append(x)
            new_mat.append(row)
        return Matrix(new_mat)
    def __add__(self, B):
        new_mat = []
        if B.d != self.d:
            raise Exception('Error: matrices of different dimensions')
        for i in range(self.d):
            row = []
            for j in range(self.d):
                row.append(self.mat[i][j]+B.mat[i][j])
            new_mat.append(row)
        return Matrix(new_mat)
    def __repr__(self):
        return '\n'.join(['['+'\t'.join([str(x) for x in self.mat[i]])+']' for i in range(self.d)])
    def __sub__(self, B):
        new_mat = []
        if B.d != self.d:
            raise Exception('Error: matrices of different dimensions')
        for i in range(self.d):
            row = []
            for j in range(self.d):
                row.append(self.mat[i][j]-B.mat[i][j])
            new_mat.append(row)
        return Matrix(new_mat)
    def tuplify(self):
        return tuple(self.mat[i][j].tuplify() for (i,j) in itertools.product(range(self.d), repeat=2))
    def __eq__(self, other):
        return self.tuplify() == other.tuplify()
def zeros(d,q):
    zero = gf0(q)
    return Matrix([[zero for i in range(d)] for j in range(d)])

def ones(d,q):
    one = gf1(q)
    return Matrix([[one for i in range(d)] for j in range(d)])

def eye(d,q):
    Z = zeros(d,q)
    for i in range(d):
        Z.mat[i][i] += 1
    return Z
def getGenerator(F):
    q = (F.p)**(F.n)
    for x in F:
        gen = True
        for i in range(1,q-1):
            if x**i in [gf0(q), gf1(q)]:
                gen = False
                break
        if gen:
            return x

def writeFieldTable(q, filename):
    file = open(filename, 'w')
    red = GF(q).reducing
    F = list(GF(q))
    nums = {}
    for i in range(q):
        nums[F[i].tuplify()] = i
    file.write('%s\n' % (q))
    file.write(' '.join(map(str, F))+'\n')
    for opname, op in [('+', lambda x,y: x+y), ('-', lambda x,y: x-y), ('*', lambda x,y: x*y)]:
        file.write(opname+'\n')
        for i in range(q):
            file.write(' '.join(map(str, [nums[op(F[i],F[j]).tuplify()] for j in range(q)]))+'\n')
    file.write('inv\n')
    file.write(' '.join(map(str, [nums[(F[j]**(-1)).tuplify()] for j in range(1,q)]))+'\n')
    file.write('chosen generator\n')
    file.write('%s\n' % (nums[getGenerator(GF(q)).tuplify()]))
    file.write('reducing polynomial: %s' % (polyify(red)))
    file.close()

def writeAllFieldTables(n, filestem):
    for q in range(2,n+1):
        if len(set(primeFact(q))) > 1:
            continue
        writeFieldTable(q, '%s%s.txt' % (filestem, q))

#writeFieldTable(729, 'gftables/gftable729.txt')
for q in [29, 31, 37, 39, 41, 43, 49]:
    writeFieldTable(q*q, 'gftables/gftable%s.txt' % (q*q))
