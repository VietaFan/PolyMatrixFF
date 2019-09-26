#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>  
#include "rcftable.h"
using namespace std;

finfield GF[50];
int ppow[50];

// returns false if such an embedding is impossible.
bool embed(matrix &M, matrix &M_ext, finfield &F, finfield &F_ext) {
	if ((F_ext.q-1)%(F.q-1) != 0) {
		zeros(M_ext, M.d);
		M_ext.d = 0;
		return false;
	}
	int ratio = (F_ext.q-1)/(F.q-1);
	int emb_map[MAX_FIELD_SIZE];
	emb_map[0] = 0;
	emb_map[1] = 1;
	for (int i=1; i<F.q-1; i++) {
		emb_map[F.genpow[i]] = F_ext.genpow[ratio*i];
	}
	M_ext.d = M.d;
	for (int i=0; i<M.d; i++) {
		for (int j=0; j<M.d; j++) {
			M_ext.mat[i][j] = emb_map[M.mat[i][j]];
		}
	}
	return true;
}

// assumes A is a matrix over F
void eig(matrix &A, finfield &F, vector<int> &eigvals, vector<vect> &eigvecs) {
	int N = 1;
	for (int i=0; i<A.d; i++) {
		N *= F.q;
	}
	vector<bool> avail(2*N, false);
	vect x, b, lx, diff;
	vector<int> candidates;
	
	int bound = 1, fact = N / F.q;
	int code;
	for (int j=0; j<A.d; j++) {
		for (int i=0; i<bound; i++) {
			code = (i*F.q+1)*fact; // [0 0 0 1 .....] - first nonzero entry = 1, then j more entries following
			candidates.push_back(code);
			avail[code] = true; 
		}
		bound *= F.q;
		fact /= F.q;
	}
	
	int lambda;
	vect u,v;
	vector<vect> ev_muls;
	for (int c: candidates) {
		if (!avail[c]) continue;
		nthvec(c, A.d, x, F);
		mul(A, x, b, F);
		
		// search for nonzero values in both, look for ratio
		for (int i=0; i<A.d; i++) {
			if (b.vals[i] > 0 && x.vals[i] > 0) {
				lambda = F.mul[b.vals[i]][F.inv[x.vals[i]]];
				scmul(lambda, x, lx, F);
				sub(lx, b, diff, F);
				if (iszero(diff)) {
					// add the new eigenvalue and eigenvector
					eigvals.push_back(lambda);
					eigvecs.push_back(x);
					ev_muls.clear();
					for (int a=1; a<F.q; a++) {
						scmul(a, x, u, F);
						ev_muls.push_back(u);
					}
					// get rid of the new things that are already linear combinations of the existing eigenvectors: these cannot be eigenvectors
					for (int c2: candidates) { 
						nthvec(c2, A.d, v, F);
						for (int j=0; j<F.q-1; j++) {
							add(v, ev_muls[j], u, F);
							if (!avail[vecindex(u, F)]) {
								avail[c2] = false;
							}
						}
					}
				}
				break;
			}
		}
	}
}

void eig_ext(matrix &A, finfield &F, finfield &F_ext, vector<int> &eigvals, vector<vect> &eigvecs) {
	matrix A_ext;
	embed(A, A_ext, F, F_ext);
	eig(A_ext, F_ext, eigvals, eigvecs);
}				

int main(int nargs, char **args) {
	int d=2,q=8;
	if (nargs > 1) {
		q = atoi(args[1]);
	}
	if (nargs > 2) {
		d = atoi(args[2]);
	}
	loadFields("gftables/gftable", GF, 49);
	get_ppows(ppow, 50);
	//int N = m_size(d, GF[q]);
	
	matrix A;
	zeros(A, d);
	//A.mat[0][0] = 2;
	A.mat[0][0] = 3;
	//A.mat[0][1] = 0;
	A.mat[1][0] = 3;
	//A.mat[1][1] = 3;
	printMat(A, GF[q]);
	
	vector<int> eigvals;
	vector<vect> eigvecs;
	eig(A, GF[q], eigvals, eigvecs);
	for (int i=0; i<eigvals.size(); i++) {
		cout << GF[q].reps[eigvals[i]] << ": [";
		for (int j=0; j<d; j++) {
			cout << GF[q].reps[eigvecs[i].vals[j]];
			if (j < d-1) cout << " ";
			else cout << "]\n";
		}
	}
	return 0;
}
