#include <iostream>
#include <cstdlib>
#include <vector>
#include <unordered_map>
#include "finfield.h"
using namespace std;
template<typename Key, typename Value>
using umap = unordered_map<Key, Value>;

void get_root_cts(int q, int k, vector<int> &cts) {
	cts.clear();
	for (int i=0; i<q; i++) {
		cts.push_back(0);
	}
	finfield Fq;
	loadField("gftables/gftable", q, Fq);
	umap<uint64_t, int> scalar_of; // scmul_idx[matindex(x*eye(2), Fq)] = x
	matrix I, M;
	eye(I, 2);
	for (int x=0; x<q; x++) {
		scmul(x, I, M, Fq);
		scalar_of[matindex(M, Fq)] = x;
	}		
	uint64_t N = m_size(2, Fq), nk;
	matrix A, Ak;
	for (uint64_t n=0; n<N; n++) {
		nthmat(n, 2, A, Fq);
		if (Fq.sub[Fq.mul[A.mat[0][0]][A.mat[1][1]]][Fq.mul[A.mat[0][1]][A.mat[1][0]]] == 0) {
			continue; // not in GL(2,Fq)
		}
		matexp(A, k, Ak, Fq);
		nk = matindex(Ak, Fq);
		if (scalar_of.count(nk)) {
			cts[scalar_of[nk]]++;
		}
	}
}

int main() {
	int qvals[] = {3,5,7,9,11,13,17,19};
	for (int q: qvals) {
		cout << "q = " << q << endl;
		finfield Fq;
		loadField("gftables/gftable", q, Fq);
		vector<int> kvals;
		//get_kvals(q, kvals);
		kvals.push_back((q+1)/2);
		for (int k: kvals) {
			vector<int> sol_cts;
			vector<int> iskpow(q, 0);
			for (int t=1; t<q; t++) {
				iskpow[Fq.genpow[(t*k)%(q-1)]] = 1;
			}
			if (k%2 == 0) {
				for (int t=1; t<q; t++) {
					if (iskpow[Fq.genpow[(t*k/2)%(q-1)]] == 0) {
						iskpow[Fq.genpow[(t*k/2)%(q-1)]] = 2;
					}
				}
			}
			get_root_cts(q, k, sol_cts);
			cout << k << ": " << iskpow << endl;
			cout << k << ": " << sol_cts << endl;
		}
	}
	return 0;
}
