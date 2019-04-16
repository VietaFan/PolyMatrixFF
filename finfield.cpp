#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include "finfield.h"
using namespace std;

ostream& operator<<(ostream &out, vector<uint64_t> &vec) {
	out << "[";
	if (vec.size() > 0) {		
		for (int i=0; i<vec.size()-1; ++i)
			out << vec[i] << ", ";
		out << vec[vec.size()-1];
	}
	out << "]";
	return out;
}

template<class T>
ostream& operator<<(ostream &out, vector<T> &vec) {
	out << "[";
	if (vec.size() > 0) {		
		for (int i=0; i<vec.size()-1; ++i)
			out << vec[i] << ", ";
		out << vec[vec.size()-1];
	}
	out << "]";
	return out;
}

int gcd(int a, int b) {
	if (a < b) return gcd(b,a);
	if (b == 0) return a;
	if (b == 1) return 1;
	return gcd(b, a%b);
}

void readField(string filename, finfield &F) {
	string label;
	ifstream fin(filename);
	fin >> F.q;
			
	for (int i=0; i<F.q; i++) {
		fin >> F.reps[i];
	}
	int f_replen = 0;
	for (int i=0; i<F.q; i++)
		if (F.reps[i].size() > f_replen)
			f_replen = F.reps[i].size();
	for (int i=0; i<F.q; i++) {
		F.padded_reps[i] = "";
		for (int j=0; j<f_replen-F.reps[i].size(); j++) {
			F.padded_reps[i] += " ";
		}
		F.padded_reps[i] += F.reps[i];
	}
	fin >> label;
	for (int i=0; i<F.q; i++) {
		F.add[i] = new int[F.q];
		F.sub[i] = new int[F.q];
		F.mul[i] = new int[F.q];
	}
	for (int i=0; i<F.q; i++) {
		for (int j=0; j<F.q; j++) {
			fin >> F.add[i][j];
		}
	}
	fin >> label;
	for (int i=0; i<F.q; i++) {
		for (int j=0; j<F.q; j++) {
			fin >> F.sub[i][j];
		}
	}
	fin >> label;
	for (int i=0; i<F.q; i++) {
		for (int j=0; j<F.q; j++) {
			fin >> F.mul[i][j];
		}
	}
	fin >> label;
	F.inv[0] = -1;
	for (int i=1; i<F.q; i++) {
		fin >> F.inv[i];
	}
	fin >> label >> label >> F.gen;
	F.genpow[0] = 1;
	for (int i=1; i<F.q; i++) {
		F.genpow[i] = F.mul[F.genpow[i-1]][F.gen];
	}
	fin >> label >> label >> F.reducing;
	fin.close();
}
	
// assumes A and B are square matrices of the same dimension
void add(matrix &A, matrix &B, matrix &result, finfield &F) {
	result.d = A.d;
	for (int i=0; i<A.d; ++i)
		for (int j=0; j<A.d; ++j)
			result.mat[i][j] = F.add[A.mat[i][j]][B.mat[i][j]];
}

// assumes A and B are square matrices of the same dimension
void sub(matrix &A, matrix &B, matrix &result, finfield &F) {
	result.d = A.d;
	for (int i=0; i<A.d; ++i)
		for (int j=0; j<A.d; ++j)
			result.mat[i][j] = F.sub[A.mat[i][j]][B.mat[i][j]];
}

// k is an element of the finite field F
void scmul(int k, matrix &A, matrix &result, finfield &F) {
	result.d = A.d;
	for (int i=0; i<A.d; ++i)
		for (int j=0; j<A.d; ++j)
			result.mat[i][j] = F.mul[k][A.mat[i][j]];
}

void zeros(matrix &A, int d) {
	A.d = d;
	for (int i=0; i<d; ++i)
		for (int j=0; j<d; ++j)
			A.mat[i][j] = 0;
}

void eye(matrix &A, int d) {
	zeros(A, d);
	for (int i=0; i<d; ++i)
		A.mat[i][i] = 1;
}

// assumes A and B are square matrices of the same dimension
void mul(matrix &A, matrix &B, matrix &result, finfield &F) {
	result.d = A.d;
	int accum;
	for (int i=0; i<A.d; i++)
		for (int j=0; j<A.d; j++) {
			accum = 0;
			for (int k=0; k<A.d; k++)
				accum = F.add[accum][F.mul[A.mat[i][k]][B.mat[k][j]]];
			result.mat[i][j] = accum;
		}
}

void copy(matrix &M, matrix &result) {
	result.d = M.d;
	for (int i=0; i<M.d; i++) {
		for (int j=0; j<M.d; j++) {
			result.mat[i][j] = M.mat[i][j];
		}
	}
}

// matrix exponentiation
void matexp(matrix &A, int k, matrix &result, finfield &F) {
	eye(result, A.d);
	matrix B, C;
	copy(A, B);
	for (int j=1; j<=k; j<<=1) {
		if (k & j) {
			copy(result, C);
			mul(B, C, result, F);
		}
		mul(B, B, C, F);
		copy(C, B);
	}
}

int trace(matrix &M, finfield &F) {
	int tr = 0;
	for (int i=0; i<M.d; i++) {
		tr = F.add[tr][M.mat[i][i]];
	}
	return tr;
}

void printMat(matrix &M, finfield &F) {
	for (int i=0; i<M.d; i++) {
		cout << "[";
		for (int j=0; j<M.d; j++) {
			cout << F.padded_reps[M.mat[i][j]];
			if (j < M.d-1) 
				cout << " ";
		}
		cout << "]\n";
	}
}

int m_size(int d, finfield &F) {
	int nmats = 1;
	for (int i=0; i<d*d; i++)
		nmats *= F.q;
	return nmats;
}

int gl_size(int d, finfield &F) {
	int qn = 1;
	for (int i=0; i<d; i++)
		qn *= F.q;
	int gl_size = 1;
	int subt = 1;
	for (int i=0; i<d; i++) {
		gl_size *= (qn-subt);
		subt *= F.q;
	}
	return gl_size;
}

void nthmat(uint64_t n, int d, matrix &M, finfield &F) {
	M.d = d;
	for (int i=0; i<d; i++) {
		for (int j=0; j<d; j++) {
			M.mat[i][j] = n % F.q;
			n /= F.q;
		}
	}
}

uint64_t matindex(matrix &M, finfield &F) {
	uint64_t n = 0, mult = 1;
	for (int i=0; i<M.d; i++) {
		for (int j=0; j<M.d; j++) {
			n += mult*M.mat[i][j];
			mult *= F.q;
		}
	}
	return n;
}

bool iszero(matrix &M) {
	for (int i=0; i<M.d; i++) 
		for (int j=0; j<M.d; j++)
			if (M.mat[i][j]) return false;
	return true;
}

// from https://rosettacode.org/wiki/Reduced_row_echelon_form#C.2B.2B
// translated from the original Python and adapted for finite fields
void rref(matrix &M, finfield &F) {
    int lead = 0, i, temp, lv, lv_inv;
    for (int r=0; r<M.d; r++) {
        if (lead >= M.d)
            return;
        i = r;
        while (M.mat[i][lead] == 0) {
            i++;
            if (i == M.d) {
                i = r;
                lead++;
                if (M.d == lead) {
                    return;
				}
			}
		}
		for (int col=0; col<M.d; col++) {
			temp = M.mat[i][col];
			M.mat[i][col] = M.mat[r][col];
			M.mat[r][col] = temp;
		}
        lv_inv = F.inv[M.mat[r][lead]];
        for (int col=0; col<M.d; col++) {
			M.mat[r][col] = F.mul[M.mat[r][col]][lv_inv];
		}
        for (int i=0; i<M.d; i++) {
            if (i != r) {
                lv = M.mat[i][lead];
                for (int col=0; col<M.d; col++) {
					M.mat[i][col] = F.sub[M.mat[i][col]][F.mul[lv][M.mat[r][col]]];
				}
			}
		}
		lead++;
	}
}

bool inv(matrix &A, matrix &Minv, finfield &F) {
    int lead = 0, i, temp, lv, lv_inv;
    matrix M, I, diff;
    copy(A, M);
    eye(Minv, M.d);
    eye(I, M.d);
    for (int r=0; r<M.d; r++) {
        if (lead >= M.d) {
			sub(M, I, diff, F);
            return iszero(diff);
        }
        i = r;
        while (M.mat[i][lead] == 0) {
            i++;
            if (i == M.d) {
                i = r;
                lead++;
                if (M.d == lead) {
					sub(M, I, diff, F);
                    return iszero(diff);
				}
			}
		}
		for (int col=0; col<M.d; col++) {
			temp = M.mat[i][col];
			M.mat[i][col] = M.mat[r][col];
			M.mat[r][col] = temp;
			temp = Minv.mat[i][col];
			Minv.mat[i][col] = Minv.mat[r][col];
			Minv.mat[r][col] = temp;
		}
        lv_inv = F.inv[M.mat[r][lead]];
        for (int col=0; col<M.d; col++) {
			M.mat[r][col] = F.mul[M.mat[r][col]][lv_inv];
			Minv.mat[r][col] = F.mul[Minv.mat[r][col]][lv_inv];
		}
        for (int i=0; i<M.d; i++) {
            if (i != r) {
                lv = M.mat[i][lead];
                for (int col=0; col<M.d; col++) {
					M.mat[i][col] = F.sub[M.mat[i][col]][F.mul[lv][M.mat[r][col]]];
					Minv.mat[i][col] = F.sub[Minv.mat[i][col]][F.mul[lv][Minv.mat[r][col]]];
				}
			}
		}
		lead++;
	}
	return true;
}
   
bool is_ppow(int q) {
	if (q < 2) {
		return false;
	}
	int p = 2;
	while (q%p) {
		p++;
	}
	while (!(q%p)) {
		q /= p;
	}
	return q == 1;
}

void get_ppows(int *ppow, int maxq) {
	int pos = 0;
	for (int q=2; q<=maxq; q++) {
		if (is_ppow(q)) {
			ppow[pos] = q;
			pos++;
		}
	}
}

void loadFields(string fpath, finfield *GF, int maxq) {
	for (int q=2; q<=maxq; q++) {
		if (!is_ppow(q)) {
			continue;
		}
		ostringstream oss;
		oss << fpath << q << ".txt";
		readField(oss.str(), GF[q]);
	}
}

void loadFieldPair(string fpath, int q, finfield &Fq, finfield &Fq2) {
	ostringstream oss;
	oss << fpath << q << ".txt";
	readField(oss.str(), Fq);
	ostringstream oss2;
	oss2 << fpath << (q*q) << ".txt";
	readField(oss2.str(), Fq2);
}
