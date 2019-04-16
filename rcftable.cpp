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
   
void getRCFs2(vector<int> &rcfs, vector<int> &rcftypes, finfield &F) {
	rcfs.clear();
	rcftypes.clear();
	int q = F.q;
	matrix M;
	zeros(M, 2);
	/* [a 0]
	 * [0 a] */
	for (int a=0; a<q; a++) {
		M.mat[0][0] = a;
		M.mat[1][1] = a;
		rcfs.push_back(matindex(M, F));
		rcftypes.push_back(0);
	}
	/* [0 b]
	 * [1 a] */
    zeros(M, 2);
	for (int a=0; a<q; a++) {
		for (int b=0; b<q; b++) {
			M.mat[0][1] = a;
			M.mat[1][0] = 1;
			M.mat[1][1] = b;
			rcfs.push_back(matindex(M, F));
			rcftypes.push_back(1);
		}
	}
}
	
void getRCFs3(vector<int> &rcfs, vector<int> &rcftypes, finfield &F) {
	rcfs.clear();
	rcftypes.clear();
	int q = F.q;
	matrix M;
	zeros(M, 3);
	/* [a 0 0]
	 * [0 a 0] 
	 * [0 0 a] */
	for (int a=0; a<q; a++) {
		M.mat[0][0] = a;
		M.mat[1][1] = a;
		M.mat[2][2] = a;
		rcfs.push_back(matindex(M, F));
		rcftypes.push_back(0);
	}
	/* [a 0 0]
	 * [0 0 -ab]
	 * [0 1 a+b]*/
    zeros(M, 3);
	for (int a=0; a<q; a++) {
		for (int b=0; b<q; b++) {
			M.mat[0][0] = a;
			M.mat[1][2] = F.sub[0][F.mul[a][b]];
			M.mat[2][1] = 1;
			M.mat[2][2] = F.add[a][b];
			rcfs.push_back(matindex(M, F));
			rcftypes.push_back(1);
		}
	}
	/* [0 0 c]
	 * [1 0 b]
	 * [0 1 a] */
	zeros(M, 3);
	for (int a=0; a<q; a++) {
		for (int b=0; b<q; b++) {
			for (int c=0; c<q; c++) {
				M.mat[0][2] = c;
				M.mat[1][0] = 1;
				M.mat[1][2] = b;
				M.mat[2][1] = 1;
				M.mat[2][2] = a;
				rcfs.push_back(matindex(M, F));
				rcftypes.push_back(2);
			}
		}
	}			
}

// we have that invtable[matindex(M)] = 0 if M is not invertible
// only okay if q^(d^2) < 2^31
void make_invtable(int d, finfield &F, int *invtable) {
	int nmats = 1;
	for (int i=0; i<d*d; i++)
		nmats *= F.q;
	matrix M;
	M.d = d;
	matrix Minv;
	for (int i=0; i<nmats; i++) {
		nthmat(i, d, M, F);
		invtable[i] = inv(M, Minv, F);
		if (invtable[i]) { // leave it as 0 if it's not invertible - Minv is meaningless then
			invtable[i] = matindex(Minv, F);
		}
	}
}

void make_rcftable(int d, finfield &F, int *invtable, int *rcftable, int *simtable) {
	if (d > 3) {
		cout << "ERROR: too big\n";
		return;
	}
	vector<int> rcfs, rcftypes;
	if (d == 2) {
		getRCFs2(rcfs, rcftypes, F);
	} else if (d == 3) {
		getRCFs3(rcfs, rcftypes, F);
	}
	int glsize = gl_size(d, F);
	int nmats = m_size(d, F);
	matrix *glnq = new matrix[glsize];
	matrix *glnq_invs = new matrix[glsize];
	int *invindex = new int[glsize];
	int k = 0;
	for (int i=0; i<nmats; i++) {
		if (invtable[i]) {
			nthmat(i, d, glnq[k], F);
			nthmat(invtable[i], d, glnq_invs[k], F);
			invindex[k] = invtable[i];
			k++;
		}
	}
	matrix R, SR, X;
	int matndx;
	for (int r=0; r<rcfs.size(); r++) {
		nthmat(rcfs[r], d, R, F);
		for (int s=0; s<glsize; s++) {
			mul(glnq[s], R, SR, F);
			mul(SR, glnq_invs[s], X, F);
			matndx = matindex(X, F);
			rcftable[matndx] = r;
			simtable[matndx] = invindex[s];
		}
	}
	delete[] glnq;
	delete[] glnq_invs;
	delete[] invindex;
}

void init_tables(int d, finfield &F, int **invtable, int **rcftable, int **simtable, 
				 vector<int> &rcfs, vector<int> &rcftypes) {
	int N = m_size(d, F);
	*invtable = new int[N];
	*rcftable = new int[N];
	*simtable = new int[N];
	if (d == 2) {
		getRCFs2(rcfs, rcftypes, F);
	} else if (d == 3) {
		getRCFs3(rcfs, rcftypes, F);
	}
	make_invtable(d, F, *invtable);
	make_rcftable(d, F, *invtable, *rcftable, *simtable);
}

void free_tables(int *invtable, int *rcftable, int *simtable) {
	delete[] invtable;
	delete[] rcftable;
	delete[] simtable;
}
