#ifndef _FINFIELD_H
#define _FINFIELD_H
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#define MAX_FIELD_SIZE 2000 
#define MAX_MAT_SIZE 5
using namespace std;


ostream& operator<<(ostream &out, vector<uint64_t> &vec);
template<class T>
ostream& operator<<(ostream &out, vector<T> &vec);

// finite field GF(q), q = p^r <= MAX_FIELD_SIZE
struct finfield {
	int q;
	string reps[MAX_FIELD_SIZE];
	string padded_reps[MAX_FIELD_SIZE];
	int* add[MAX_FIELD_SIZE];
	int* sub[MAX_FIELD_SIZE];
	int* mul[MAX_FIELD_SIZE];
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

void readField(string filename, finfield &F);
	
// assumes A and B are square matrices of the same dimension
void add(matrix &A, matrix &B, matrix &result, finfield &F);
void sub(matrix &A, matrix &B, matrix &result, finfield &F);
void mul(matrix &A, matrix &B, matrix &result, finfield &F);
void matexp(matrix &A, int k, matrix &result, finfield &F);
// k is an element of the finite field F
void scmul(int k, matrix &A, matrix &result, finfield &F);

void zeros(matrix &A, int d);
void eye(matrix &A, int d);
bool iszero(matrix &M);
void copy(matrix &M, matrix &result);

int trace(matrix &M, finfield &F);
void rref(matrix &M, finfield &F);
bool inv(matrix &A, matrix &Minv, finfield &F);

void nthmat(uint64_t n, int d, matrix &M, finfield &F);
uint64_t matindex(matrix &M, finfield &F);
void printMat(matrix &M, finfield &F);

int m_size(int d, finfield &F);
int gl_size(int d, finfield &F);

void loadFields(string fpath, finfield *GF, int maxq);
void loadFieldPair(string fpath, int q, finfield &Fq, finfield &Fq2);
void get_ppows(int *ppow, int maxq);
#endif
