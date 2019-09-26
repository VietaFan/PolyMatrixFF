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
	cout << nmats << endl;
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



int main(int nargs, char **args) {
	int d=2,q=2;
	if (nargs > 1) {
		q = atoi(args[1]);
	}
	if (nargs > 2) {
		d = atoi(args[2]);
	}
	loadFields("gftables/gftable", GF, 49);
	get_ppows(ppow, 50);
	int N = m_size(d, GF[q]);
	
	//for (int i=0; i<23; i++) {
	//	cout << ppow[i] << " " << countInvKSols2(2, GF[ppow[i]], 2) << endl;
	//}
	
	int *invtable, *rcftable, *simtable;
	vector<int> rcfs, rcftypes;
	init_tables(d, GF[q], &invtable, &rcftable, &simtable, rcfs, rcftypes);
	
	unordered_map<int,int> rcfcts;
	for (int i=0; i<N; i++) {
		rcfcts[rcfs[rcftable[i]]]++;
	}
	
	matrix R;
	for (int rcf: rcfs) {
		nthmat(rcf, d, R, GF[q]);
		printMat(R, GF[q]);
		cout << "count: " << rcfcts[rcf] << endl;
	}
	
	/*matrix M;
	for (int rcf: rcfs) {
		nthmat(rcf, d, M, GF[q]);
		printMat(M, GF[q]);
		cout << endl;
	}
	*/
	matrix X, Y, X2, ainvX2, Y2, diff, XY, YX;
	cout << "d = " << d << ", q = " << q << endl;

	for (int a=1; a<q; a++) {
		cout << "a = " << GF[q].reps[a] << endl;
		for (int rcf: rcfs) {
			nthmat(rcf, d, X, GF[q]);
			
			mul(X, X, X2, GF[q]);
			scmul(GF[q].inv[a], X2, ainvX2, GF[q]);
			// ainvX2 = a^(-1)X^2
			
			int nsols = 0, ncsols = 0;
			for (int yidx=0; yidx<N; yidx++) {
				nthmat(yidx, d, Y, GF[q]);
				mul(Y, Y, Y2, GF[q]);
				sub(ainvX2, Y2, diff, GF[q]);
				if (iszero(diff)) {
					nsols++;
					mul(X, Y, XY, GF[q]);
					mul(Y, X, YX, GF[q]);
					sub(XY, YX, diff, GF[q]);
					if (iszero(diff)) ncsols++;
				}
			}
			
			if (nsols > ncsols) {			
				cout << "X = \n";
				printMat(X, GF[q]);
				
				cout << "# solutions Y: " << nsols << endl;
				cout << "# commuting solutions Y: " << ncsols << endl;
			}
		}
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
	printMat(M4, GF[8]);
	cout << countNilpotent(2, GF[3], 5);
	cout << countKSols2(3, GF[3], 3) << endl;
	/*
	for (int i=0; i<23; i++) {
		cout << 3 << " " << ppow[i] << " " << countKSols2(3, GF[ppow[i]], ppow[i]) << endl; //<< " " << (ppow[i]*ppow[i]*ppow[i]*(ppow[i]*ppow[i]+ppow[i]-1)) << endl;
	}
	
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

