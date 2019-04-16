#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#define MAX_FIELD_SIZE 64
#define MAX_MAT_SIZE 5
using namespace std;

ostream& operator<<(ostream &ostr, const vector<uint64_t> &v) {
	ostr << "[";
	for (int i=0; i<v.size(); i++) {
		ostr << v[i];
		if (i < v.size()-1) {
			ostr << ", ";
		}
	}
	ostr << "]";
	return ostr;
}

// finite field GF(q), q = p^r <= MAX_FIELD_SIZE
struct finfield {
	int q;
	string reps[MAX_FIELD_SIZE];
	string padded_reps[MAX_FIELD_SIZE];
	int add[MAX_FIELD_SIZE][MAX_FIELD_SIZE];
	int sub[MAX_FIELD_SIZE][MAX_FIELD_SIZE];
	int mul[MAX_FIELD_SIZE][MAX_FIELD_SIZE];
	int inv[MAX_FIELD_SIZE];
	int gen;
	int genpow[MAX_FIELD_SIZE];
	string reducing;
};

// square matrix of dimension d (can be over any finite field)
struct matrix {
	int d;
	int mat[MAX_MAT_SIZE][MAX_MAT_SIZE];
};

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
		for (int j=0; j<f_replen-F.reps[i].size(); j++) {
			F.padded_reps[i] += " ";
		}
		F.padded_reps[i] += F.reps[i];
	}
	fin >> label;
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

// counts the number of solutions to X^2-Y^2 = 0 in M_d(F)
uint64_t countQuadSolsDirect(int d, finfield &F) {
	matrix X, Y, temp1, temp2, ans;
	int nmats = 1;
	for (int i=0; i<d*d; i++) {
		nmats *= F.q;
	}
	uint64_t nsols = 0;
	for (int i=0; i<nmats; i++) {
		nthmat(i, d, X, F); 
		mul(X, X, temp1, F);
		for (int j=0; j<nmats; j++) {
			nthmat(j, d, Y, F);
			mul(Y, Y, temp2, F);
			sub(temp1, temp2, ans, F);
			if (iszero(ans)) {
				nsols++;
			}
		}
	}
	return nsols;
}

// counts the number of solutions to X^2-Y^2 = 0 in M_d(F)
// by finding the square of everything and multiplying class pairs
uint64_t countQuadSols2(int d, finfield &F) {
	matrix X, sq;
	unordered_map<uint64_t, uint64_t> sqrtcount;
	uint64_t nmats = 1;
	for (int i=0; i<d*d; i++) {
		nmats *= F.q;
	}
	uint64_t sq_index;
	for (uint64_t i=0; i<nmats; i++) {
		nthmat(i, d, X, F); 
		mul(X, X, sq, F);
		sq_index = matindex(sq, F);
		if (!sqrtcount.count(sq_index)) {
			sqrtcount[sq_index] = 0;
		}
		sqrtcount[sq_index]++;
	}
	uint64_t nsols = 0;
	for (int i=0; i<nmats; i++) {
		nsols += sqrtcount[i]*sqrtcount[i];
	}
	return nsols;
}

// counts the number of solutions to X^2-Y^2 = 0 in M_d(F)
// by finding the square of everything and multiplying class pairs
uint64_t countKSols2(int d, finfield &F, int k) {
	matrix X, sq;
	unordered_map<uint64_t, uint64_t> sqrtcount;
	uint64_t nmats = 1;
	for (int i=0; i<d*d; i++) {
		nmats *= F.q;
	}
	uint64_t sq_index;
	for (uint64_t i=0; i<nmats; i++) {
		nthmat(i, d, X, F); 
		matexp(X, k, sq, F);
		sq_index = matindex(sq, F);
		if (!sqrtcount.count(sq_index)) {
			sqrtcount[sq_index] = 0;
		}
		sqrtcount[sq_index]++;
	}
	uint64_t nsols = 0;
	for (int i=0; i<nmats; i++) {
		nsols += sqrtcount[i]*sqrtcount[i];
	}
	return nsols;
}


// directly counts the number of solutions to [[X,Y],Z] = 0
uint64_t countSolsDirect(int d, finfield &F) {
	matrix X, Y, Z;
	matrix temp1, temp2, temp3, temp4, temp5, ans;
	int nmats = 1;
	for (int i=0; i<d*d; i++) {
		nmats *= F.q;
	}
	uint64_t nsols = 0;
	for (int i=0; i<nmats; i++) {
		nthmat(i, d, X, F);
		for (int j=0; j<nmats; j++) {
			nthmat(j, d, Y, F);
			for (int k=0; k<nmats; k++) {
				nthmat(k, d, Z, F);
				mul(X, Y, temp1, F);
				mul(Y, X, temp2, F);
				sub(temp1, temp2, temp3, F); //temp3 = [X,Y]
				mul(temp3, Z, temp4, F); 
				mul(Z, temp3, temp5, F);
				sub(temp4, temp5, ans, F); // ans = [temp3, Z]
				if (iszero(ans)) {
					nsols++;
				}
			}
		}
	}
	return nsols;
}

int numCommuting(matrix &X, finfield &F) {
	uint64_t nmats = 1;
	for (int i=0; i<X.d*X.d; i++) {
		nmats *= F.q;
	}
	matrix Y, XY, YX, comm;
	int nsol = 0;
	for (int i=0; i<nmats; i++) {
		nthmat(i, X.d, Y, F);
		mul(X,Y,XY,F);
		mul(Y,X,YX,F);
		sub(XY,YX,comm,F);
		if (iszero(comm))
			nsol++;
	}
	return nsol;
}

int numCommutingGL(matrix &X, finfield &F, int *invtable) {
	uint64_t nmats = 1;
	for (int i=0; i<X.d*X.d; i++) {
		nmats *= F.q;
	}
	matrix Y, XY, YX, comm;
	int nsol = 0;
	for (int i=0; i<nmats; i++) {
		if (!invtable[i]) continue;
		nthmat(i, X.d, Y, F);
		mul(X,Y,XY,F);
		mul(Y,X,YX,F);
		sub(XY,YX,comm,F);
		if (iszero(comm))
			nsol++;
	}
	return nsol;
}

uint64_t numBracketPairs(matrix &R, finfield &F) {
	int n = m_size(R.d, F);
	matrix *mats = new matrix[n];
	matrix temp1, temp2, temp3, diff;
	uint64_t count = 0;
	for (int i=0; i<n; i++) {
		nthmat(i, R.d, mats[i], F);
	}
	for (int i=0; i<n; i++) {
		for (int j=0; j<n; j++) {
			mul(mats[i],mats[j],temp1,F);
			mul(mats[j],mats[i],temp2,F);
			sub(temp1,temp2,temp3,F);
			sub(temp3, R, diff, F);
			if (iszero(diff)) {
				count++;
			}
		}
	}
	delete[] mats;
	return count;
}

finfield GF[50];
int ppow[] = {2,3,4,5,7,8,9,11,13,16,17,19,23,25,27,29,31,32,37,41,43,47,49};

void loadFields() {
	for (int i=0; i<23; i++) {
		ostringstream oss;
		oss << "gftables/gftable" << ppow[i] << ".txt";
		readField(oss.str(), GF[ppow[i]]);
	}
}	

void showCommutingCounts(int d, int q) {
	cout << "numbers of matrices with each number of commuting seconds in M_" << d << "(F_" << q << ")...\n";
	matrix M;
	int nmats = 1;
	for (int i=0; i<d*d; i++) {
		nmats *= q;
	}
	
	int *commcts = new int[nmats];
	for (int i=0; i<nmats; i++) {
		nthmat(i, d, M, GF[q]);
		commcts[i] = numCommuting(M, GF[q]);
	}
	vector<int> cctvals;
	unordered_map<int, int> cct_counts;
	for (int i=0; i<nmats; i++) {
		if (!cct_counts.count(commcts[i])) {
			cct_counts[commcts[i]] = 1;
			cctvals.push_back(commcts[i]);
		} else {
			cct_counts[commcts[i]]++;
		}
	}
	sort(cctvals.begin(), cctvals.end(), greater<int>());
	cout << "#comm count\n";
	for (int i=0; i<cctvals.size(); i++) {
		cout << cctvals[i] << " " << cct_counts[cctvals[i]] << endl;
	}
	delete commcts;
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
}


// counts the number of solutions to X^2-Y^2 = 0 in GL_d(F)
// by finding the square of everything and multiplying class pairs
uint64_t countInvKSols2(int d, finfield &F, int k) {
	matrix X, Y, sq;
	unordered_map<uint64_t, uint64_t> sqrtcount;
	uint64_t nmats = 1;
	for (int i=0; i<d*d; i++) {
		nmats *= F.q;
	}
	uint64_t sq_index;
	int inv_idx;
	for (uint64_t i=0; i<nmats; i++) {
		nthmat(i, d, X, F);
		//printMat(X, F);
		//cout << inv(X,Y,F) << endl; 
		inv_idx = inv(X,Y,F);
		if (inv_idx == 0) continue;
		nthmat(i, d, X, F);
		//cout << "hello\n";
		matexp(X, k, sq, F);
		sq_index = matindex(sq, F);
		if (!sqrtcount.count(sq_index)) {
			sqrtcount[sq_index] = 0;
		}
		sqrtcount[sq_index]++;
	}
	uint64_t nsols = 0;
	for (int i=0; i<nmats; i++) {
		nsols += sqrtcount[i]*sqrtcount[i];
	}
	return nsols;
}

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

int get_period(matrix &X, int *rcftable, finfield &F) {
	unordered_set<int> rcfpos;
	matrix Y,Z;
	copy(X,Y);
	while (!rcfpos.count(rcftable[matindex(Y, F)])) {
		rcfpos.insert(rcftable[matindex(Y, F)]);
		mul(X,Y,Z,F);
		copy(Z,Y);
	}
	int form = rcftable[matindex(Y,F)];
	int period = 1;
	mul(X,Y,Z,F);
	copy(Z,Y);
	while (rcftable[matindex(Y,F)] != form) {
		mul(X,Y,Z,F);
		copy(Z,Y);
		period++;
	}
	return period;
}

void show_rcf_cycles(int d, finfield &F) {
	int *invtable, *rcftable, *simtable;
	vector<int> rcfs, rcftypes;
	init_tables(d, F, &invtable, &rcftable, &simtable, rcfs, rcftypes);
	matrix X, Y, Z;
	for (int i=0; i<rcfs.size(); i++) {
		nthmat(rcfs[i], d, X, F);
		copy(X,Y);
		//cout << "R_" << i << " period = " << get_period(X, rcftable, F) << endl;
		if (inv(X,Z,F)) continue;
		printMat(X, F);
		cout << (inv(X, Z, F)?"invertible":"not invertible") << endl;
		cout << "R_" << i << " -> ";
		for (int exp=2; exp<20; exp++) {
			mul(X,Y,Z,F);
			copy(Z,Y);
			//nthmat(rcfs[rcftable[matindex(Y, F)]], d, Z, F); // Z = rcf(Y)
			cout << "R_" << rcftable[matindex(Y, F)] << " -> ";
			//printMat(Z, F);
		}
		cout << "...\n";
	}
}

uint64_t count_sols(int d, finfield &F, bool verbose) {
	int *invtable, *rcftable, *simtable;
	vector<int> rcfs, rcftypes;
	init_tables(d, F, &invtable, &rcftable, &simtable, rcfs, rcftypes);
	vector<uint64_t> c, z, s;
	matrix *rcf_mats = new matrix[rcfs.size()];
	for (int r=0; r<rcfs.size(); r++) {
		nthmat(rcfs[r], d, rcf_mats[r], F);
	}
	for (int r=0; r<rcfs.size(); r++) {
		c.push_back(numCommuting(rcf_mats[r], F));
		z.push_back(numCommutingGL(rcf_mats[r], F, invtable));
	}
	int nmats = m_size(d, F);
	matrix *mats = new matrix[nmats];
	for (int i=0; i<nmats; i++) {
		nthmat(i, d, mats[i], F);
	}
	vector<vector<bool> > chi;
	for (int r=0; r<rcfs.size(); r++) {
		chi.push_back(vector<bool>(nmats, 0));
	}
	matrix temp1, temp2, X;
	for (int r=0; r<rcfs.size(); r++) {
		for (int i=0; i<nmats; i++) {
			mul(rcf_mats[r], mats[i], temp1, F);
			mul(mats[i], rcf_mats[r], temp2, F);
			sub(temp1, temp2, X, F);
			chi[r][matindex(X, F)] = 1;
		}
	}
	
	uint64_t next_sval;
	for (int r=0; r<rcfs.size(); r++) {
		next_sval = 0;
		for (int i=0; i<nmats; i++) {
			next_sval += chi[r][i]*c[rcftable[i]];
		}
		s.push_back(next_sval);
	}
	
	uint64_t nsol = 0;
	uint64_t glsize = gl_size(d, F);
	for (int r=0; r<rcfs.size(); r++) {
		if (verbose) {
			cout << "c(R) = " << c[r] << ", |GL(n,F_q)|/z(R) = " << glsize/z[r] << ", s(R) = " << s[r] << ", # sol from R = " << glsize*c[r]*s[r]/z[r] << "\n";
			printMat(rcf_mats[r], F);
		}
		nsol += glsize*c[r]*s[r]/z[r];
	}
	
	delete[] rcf_mats;
	delete[] mats;
	return nsol;
}

void getClassSizes(int d, finfield &F, vector<int> &rcfs, int *invtable, vector<uint64_t> &csizes) {
	uint64_t glsize = gl_size(d, F);
	matrix X;
	for (int r=0; r<rcfs.size(); r++) {
		nthmat(rcfs[r], d, X, F);
		csizes.push_back(glsize/numCommutingGL(X, F, invtable));
	}
}

uint64_t countQuadSols3(int d, finfield &F) {
	int *invtable, *rcftable, *simtable;
	vector<int> rcfs, rcftypes;
	init_tables(d, F, &invtable, &rcftable, &simtable, rcfs, rcftypes);
	vector<uint64_t> csizes;
	getClassSizes(d, F, rcfs, invtable, csizes);
	vector<int> sqclass;
	matrix R, R2;
	for (int r=0; r<rcfs.size(); r++) {
		nthmat(rcfs[r], d, R, F);
		mul(R, R, R2, F);
		sqclass.push_back(rcftable[matindex(R2, F)]);
	}
	unordered_map<int, uint64_t> start_counts;
	for (int r=0; r<rcfs.size(); r++) {
		if (!start_counts.count(sqclass[r]))
			start_counts[sqclass[r]] = 0;
		start_counts[sqclass[r]] += csizes[r];
	}
	uint64_t nsol = 0;
	for (int r=0; r<rcfs.size(); r++) {
		if (!start_counts.count(r))
			continue;
		nsol += start_counts[r] * start_counts[r]/csizes[r];
	}
	return nsol;
}

uint64_t countNilpotent(int d, finfield &F, int k) {
	uint64_t msize = m_size(d, F);
	matrix M, Mk;
	uint64_t nilct=0;
	for (uint64_t i=0; i<msize; i++) {
		nthmat(i, d, M, F);
		matexp(M, k, Mk, F);
		if (iszero(Mk)) {
			nilct++;
			printMat(M,F);
			cout << endl;
		}
	}
	return nilct;
}

uint64_t countKthSols3(int d, finfield &F, vector<int> &rcfs, int *invtable, int *rcftable, int k) {
	vector<uint64_t> csizes;
	getClassSizes(d, F, rcfs, invtable, csizes);
	vector<int> sqclass;
	matrix R, Rk, Z;
	for (int r=0; r<rcfs.size(); r++) {
		nthmat(rcfs[r], d, R, F);
		copy(R, Rk);
		for (int i=1; i<k; i++) {
			mul(R, Rk, Z, F);
			copy(Z, Rk);
		}
		sqclass.push_back(rcftable[matindex(Rk, F)]);
	}
	unordered_map<int, uint64_t> start_counts;
	for (int r=0; r<rcfs.size(); r++) {
		if (!start_counts.count(sqclass[r]))
			start_counts[sqclass[r]] = 0;
		start_counts[sqclass[r]] += csizes[r];
	}
	uint64_t nsol = 0;
	for (int r=0; r<rcfs.size(); r++) {
		if (!start_counts.count(r))
			continue;
		nsol += start_counts[r] * start_counts[r]/csizes[r];
	}
	return nsol;
}

void countAllKthSols3(int d, finfield &F, vector<int> &rcfs, int *invtable, int *rcftable, int k, vector<uint64_t> &nsol) {
	vector<uint64_t> csizes;
	getClassSizes(d, F, rcfs, invtable, csizes);
	vector<uint64_t> sqclass[10000];
	matrix R, Rk, Z;
	for (int r=0; r<rcfs.size(); r++) {
		nthmat(rcfs[r], d, R, F);
		copy(R, Rk);
		for (int i=1; i<k; i++) {
			mul(R, Rk, Z, F);
			copy(Z, Rk);
			sqclass[i+1].push_back(rcftable[matindex(Rk, F)]);
		}
	}
	unordered_map<int, uint64_t> start_counts[10000];
	for (int r=0; r<rcfs.size(); r++) {
		for (int i=2; i<=k; i++) {
			if (!start_counts[i].count(sqclass[i][r]))
				start_counts[i][sqclass[i][r]] = 0;
			start_counts[i][sqclass[i][r]] += csizes[r];
		}
	}
	nsol.clear();
	for (int i=1; i<=k+1; i++) {
		nsol.push_back(0);
	}
	for (int r=0; r<rcfs.size(); r++) {
		for (int i=2; i<=k; i++) {
			if (!start_counts[i].count(r))
				continue;
			nsol[i] += start_counts[i][r] * start_counts[i][r]/csizes[r];
		}
	}
}

int gcd(int a, int b) {
	if (a < b) return gcd(b,a);
	if (b == 0) return a;
	if (b == 1) return 1;
	return gcd(b, a%b);
}

int main(int nargs, char **args) {
	loadFields();
	for (int i=0; i<23; i++) {
		cout << ppow[i] << " " << countInvKSols2(2, GF[ppow[i]], 2) << endl;
	}
	/*matrix M;
	nthmat(34371865, 3, M, GF[8]);
	cout << "M = \n";
	printMat(M, GF[8]);
	
	matrix M2, M3;
	copy(M, M2);
	for (int i=1; i<17; i++) {
		mul(M, M2, M3, GF[8]);
		copy(M3, M2);
	}
	cout << "M^17 = \n";
	printMat(M2, GF[8]);
	
	matrix M4;
	matexp(M, 17, M4, GF[8]);
	cout << "M^17 = \n";
	printMat(M4, GF[8]);*/
	//cout << countNilpotent(2, GF[3], 5);
	//cout << countKSols2(3, GF[7], 7) << endl;
	
	//for (int i=0; i<23; i++) {
	//	cout << 3 << " " << ppow[i] << " " << countKSols2(3, GF[ppow[i]], ppow[i]) << endl; //<< " " << (ppow[i]*ppow[i]*ppow[i]*(ppow[i]*ppow[i]+ppow[i]-1)) << endl;
	//}
	
	/*
	int dvals[] = {2,2,2,2,2,2,3,3,3,4};
	int qvals[] = {2,3,4,5,7,8,2,3,4,2};
	for (int i=6; i<7; i++) {
		cout << 3 << " " << ppow[i];
		for (int k=2; k<6; k++) {
			cout << " " << countKSols2(3, GF[ppow[i]], k);
		}
		cout << endl;			
	}
	/*for (int k=2; k<100; k++) {
		cout << 2 << " " << 2 << " " << k << " " << gcd(k, gl_size(2, GF[2])) << " " << countKSols2(2, GF[2], k) << endl;
	}
	
	/*for (int k=2; k<100; k++) {
		cout << "k = " << k << ": " << countKthSols3(2, GF[2], k) << endl;
	}*/
	/*int q = 2;
	int d = 3;
	int K = 3000;
	
	int *invtable, *rcftable, *simtable;
	vector<int> rcfs, rcftypes;
	init_tables(d, GF[q], &invtable, &rcftable, &simtable, rcfs, rcftypes);
	
	int t;
	vector<uint64_t> counts;
	unordered_map<uint64_t, vector<uint64_t> > solct_map;
	vector<uint64_t> vals;
	/*for (int k=2; k<101; k++) {
		t = countKthSols3(d, GF[q], rcfs, invtable, rcftable, k);
		if (!solct_map.count(t)) {
			solct_map[t] = vector<int>();
			vals.push_back(t);
		}
		solct_map[t].push_back(k);
	}
	sort(vals.begin(), vals.end());
	
	countAllKthSols3(d, GF[q], rcfs, invtable, rcftable, K, counts);
	for (int k=2; k<K+1; k++) {
		if (!solct_map.count(counts[k])) {
			solct_map[counts[k]] = vector<uint64_t>();
			vals.push_back(counts[k]);
		}
		solct_map[counts[k]].push_back(k);
	}
	sort(vals.begin(), vals.end());
	for (int i=0; i<vals.size(); i++) {
		cout << vals[i] << " " << solct_map[vals[i]][0] << endl;//(solct_map[vals[vals.size()-1]][0] % solct_map[vals[i]][0]) << " "  << solct_map[vals[i]] << endl;
	}
	cout << vals.size() << endl;
	
	
	show_rcf_cycles(3,GF[3]);
	/*
	for (int i=0; i<10; i++) {
		cout << "M_" << dvals[i] << "(F_" << qvals[i] << "): ";
		cout << "conjectured: " << countKthSols3(dvals[i], GF[qvals[i]], 2);
		cout << " actual: " << countQuadSols2(dvals[i], GF[qvals[i]]) << endl;
	}
	//cout << countQuadSols2(5,GF[2]) << endl;
	//show_rcf_cycles(3,GF[2]);
	/*
	for (int i=0; i<14; i++) {
		uint64_t nsol = count_sols(2, GF[ppow[i]], false);
		cout << "(d,q) = (2," << ppow[i] << "); # solutions = " << nsol << endl;
	}
	/*
	for (int i=0; i<4; i++) {
		cout << "(d,q) = (3," << ppow[i] << "); # solutions = " << count_sols(3, GF[ppow[i]], true) << endl;
	}
	*/
	//cout << "conjectured # solutions: " << count_sols(2, GF[4]) << endl;
	//cout << "actual # solutions: " << countSolsDirect(2, GF[4]) << endl;
	/*matrix M;
	nthmat(34371865, 3, M, GF[8]);
	cout << "M = \n";
	printMat(M, GF[8]);
	matrix M2;
	copy(M, M2);
	rref(M2, GF[8]);
	cout << "rref(M) = \n";
	printMat(M2, GF[8]);
	matrix Minv;
	inv(M, Minv, GF[8]);
	cout << "inv(M) = \n";
	printMat(Minv, GF[8]);
	matrix prod;
	mul(M, Minv, prod, GF[8]);
	cout << "M * inv(M) = \n";
	printMat(prod, GF[8]);
	mul(Minv, M, prod, GF[8]);
	cout << "inv(M) * M = \n";
	printMat(prod, GF[8]);
	
	int *invtable, *rcftable, *simtable;
	vector<int> rcfs, rcftypes;
	//init_tables(3, GF[3], invtable, rcftable, simtable, rcfs, rcftypes);
	int q = 2; int d = 2;
	getRCFs2(rcfs, rcftypes, GF[q]);
	invtable = new int[m_size(d, GF[q])];
	make_invtable(d, GF[q], invtable);
	uint64_t nsol = 0;
	int c,z,s;
	for (int i=0; i<rcfs.size(); i++) {
		nthmat(rcfs[i], d, M, GF[q]);
		if (trace(M, GF[q]) != 0) continue;
		c = numCommuting(M, GF[q]);
		z = numCommutingGL(M, GF[q], invtable);
		s = numBracketPairs(M, GF[q]);
		cout << "type(R) = " << rcftypes[i] << ", c(R) = " << c << ", z(R) = " << z << ", s(R) = " << s << ", R =" <<  endl;
		nsol += gl_size(d, GF[q])*c*s/z;
		printMat(M, GF[q]);
	}
	delete[] invtable;
	cout << "conjectured # triples = " << nsol << endl;
	//cout << "actual # triples = " << countSolsDirect(d, GF[q]) << endl;
	return 0;
	int x;
	matrix A, Ainv, Arcf, Asim;
	while (true) {
		cin >> x;
		nthmat(x, 3, A, GF[3]);
		cout << "matrix #" << x << " in M_3(GF(3)) is A = \n";
		printMat(A, GF[3]);
		nthmat(invtable[x], 3, Ainv, GF[5]);
		cout << "A^(-1) = \n";
		printMat(Ainv, GF[3]);
		cout << "rcf(A) = \n";
		nthmat(rcfs[rcftable[x]], 3, Arcf, GF[3]);
		printMat(Arcf, GF[3]);
		cout << "sim(A) = \n";
		nthmat(simtable[x], 3, Asim, GF[3]);
		printMat(Asim, GF[3]);
	}
	free_tables(invtable, rcftable, simtable);
	
	/*nthmat(1872313, 3, M, GF[8]);
	cout << "M = \n";
	printMat(M, GF[8]);
	cout << matindex(M, GF[8]) << endl;
	matrix I;
	eye(I, 3);
	matrix M2;
	add(M, I, M2, GF[8]);
	cout << "M+I = \n";
	printMat(M2, GF[8]);
	cout << GF[8].reducing << "\n";*/
	
	//for (int i=6; i<14; i++) {
	//	cout << ppow[i] << " " << countQuadSols2(3, GF[ppow[i]]) << endl;
	//}
	//for (int d=2; d<6; ++d)
	//	cout << countQuadSols2(d,GF[2]) << endl;
	
	//for (int i=0; i<14; i++) {
	//	cout << ppow[i] << " " << countQuadSols2(2, GF[ppow[i]]) << endl;
	//}
	//cout << countQuadSolsDirect(2, GF[8]) << " " << countQuadSols2(2, GF[8]) << endl;
	/*for (int i=0; i<16; i++) {
		nthmat(i, 2, M, GF[2]);
		printMat(M, GF[2]);
		cout << "above matrix commutes with " << numCommuting(M, GF[2]) << " matrices\n";
	}*/
	
	//showCommutingCounts(2, 2);
	//showCommutingCounts(3, 2);
	//showCommutingCounts(3, 3);
	//showCommutingCounts(2, 7);
}

