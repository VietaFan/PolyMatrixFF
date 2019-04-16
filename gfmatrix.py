from finfield import gf0,gf1,gfx, GF_elt, GF

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
            return Matrix([[B*x for x in B.mat[i]] for i in range(self.d)])
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

