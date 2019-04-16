#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <cmath>
#include "finfield.h"
#include "glfqchar.h"
#include <Eigen/Dense>

using namespace std; 

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


void get_kth_root_cts(int d, finfield &F, int k, vector<matrix> &elts, vector<uint64_t> &counts) {
	matrix X, Y, Xk;
	unordered_map<uint64_t, uint64_t> krt_ct;
	uint64_t nmats = 1;
	for (int i=0; i<d*d; i++) {
		nmats *= F.q;
	}
	uint64_t kpow_idx;
	int inv_idx;
	for (uint64_t i=0; i<nmats; i++) {
		nthmat(i, d, X, F);
		inv_idx = inv(X,Y,F);
		if (inv_idx == 0) continue;
		nthmat(i, d, X, F);
		matexp(X, k, Xk, F);
		kpow_idx = matindex(Xk, F);
		if (!krt_ct.count(kpow_idx)) {
			krt_ct[kpow_idx] = 0;
		}
		krt_ct[kpow_idx]++;
	}
	uint64_t elt_idx;
	for (int i=0; i<elts.size(); i++) {
		elt_idx = matindex(elts[i], F);
		if (!krt_ct.count(elt_idx)) {
			counts.push_back(0);
		} else {
			counts.push_back(krt_ct[elt_idx]);
		}
	}
}

bool decompose_kr(int k, finfield &Fq, Eigen::MatrixXd &char_mat, vector<matrix> &reprs) {
	vector<uint64_t> krcts;
	get_kth_root_cts(2, Fq, k, reprs, krcts);
	Eigen::VectorXd v(2*(Fq.q*Fq.q-1));
	for (int i=0; i<Fq.q*Fq.q-1; i++) {
		v(2*i) = krcts[i];
		v(2*i+1) = 0;
	}
	//VectorXd x = char_mat.colPivHouseholderQr().solve(v);
	Eigen::VectorXd x = char_mat.fullPivLu().solve(v);
	Eigen::RowVectorXd y = x.transpose();
	for (int i=0; i<Fq.q*Fq.q-1; i++) {
		y(i) = fix_roundoff(y(i));
	}
	//cout << "eigen-computed coeffs for k = " << k << endl;
	//cout << y << endl;
	for (int i=0; i<Fq.q*Fq.q-1; i++) {
		if (y(i) < -1e-12 || fabs(y(i)-round(y(i))) > 1e-12) {
			return false;
		}
	}
	double pred_nsol = 0;
	for (int i=0; i<Fq.q*Fq.q-1; i++) {
		pred_nsol += y(i)*y(i);
	}
	pred_nsol *= gl_size(2, Fq);
	if (pred_nsol != countInvKSols2(2, Fq, k)) {
		cout << "the total number of solutions is not |GL|*(m1^2+...+mk^2)\n";
	}
	return true;
	/*cout << "k = " << k << endl;
	cout << y << endl;
	double pred_nsol = 0;
	for (int i=0; i<Fq.q*Fq.q-1; i++) {
		pred_nsol += y(i)*y(i);
	}
	pred_nsol *= gl_size(2, Fq);
	cout << "predicted # solutions: " << pred_nsol << endl;*/
}	

bool decompose_kr2(int k, finfield &Fq, vector<matrix> &reprs, 
  vector<vector<cnum> > &irrep_chars, vector<cnum> &components) {
	vector<uint64_t> krcts;
	get_kth_root_cts(2, Fq, k, reprs, krcts);
	vector<cnum> chi(Fq.q*Fq.q-1);
	for (int i=0; i<krcts.size(); i++) {
		chi[i].re = krcts[i];
		chi[i].im = 0;
	}
	vector<uint64_t> clsizes;
	get_class_sizes(Fq.q, clsizes);
	components.clear();
	uint64_t gls = gl_size(2,Fq);
	for (int i=0; i<Fq.q*Fq.q-1; i++) {
		components.push_back(inner_product(chi, irrep_chars[i], clsizes, gls));
	}
	for (cnum z: components) {
		if (fabs(z.im) > 1e-12 || fabs(z.re-round(z.re)) > 1e-12 || z.re < -1e-12) {
			if (fabs(z.re-round(z.re)) < 0.001 && fabs(z.re-round(z.re)) > 1e-12) {
				cout << "q=" << Fq.q << ", k=" << k << ", margin: " << fabs(z.re-round(z.re)) << endl;
			}
			if (fabs(z.im) < 0.001 && fabs(z.im) > 1e-12) {
				cout << "q=" << Fq.q << ", k=" << k << ", margin: " << fabs(z.im) << endl;
			}
			if (z.re > -0.001 && z.re < 1e-12) {
				cout << "q=" << Fq.q << ", k=" << k << ", margin: " << -z.re << endl;
			}
			return false;
		}
	}
	uint64_t pred_nsol = 0;
	for (int i=0; i<Fq.q*Fq.q-1; i++) {
		pred_nsol += (uint64_t) (components[i].re*components[i].re+1e-12);
	}
	pred_nsol *= gl_size(2, Fq);
	if (pred_nsol != countInvKSols2(2, Fq, k)) {
		cout << "the total number of solutions is not |GL|*(m1^2+...+mk^2)\n";
		cout << "predicted # solutions = " << pred_nsol << endl;
		cout << "actual # solutions = " << countInvKSols2(2, Fq, k) << endl;
		cout << "q = " << Fq.q << ", k = " << k << endl;
	}
	return true;
}

void get_kseries(finfield &Fq, vector<uint64_t> &ks) {
	int gls = gl_size(2,Fq);
	for (int k=1; k<=gls; k++) {
		if (gls%k == 0) {
			ks.push_back(k);
		}
	}
}

ostream& operator<<(ostream &out, vector<cnum> &vec) {
	out << "[";
	if (vec.size() > 0) {		
		for (int i=0; i<vec.size()-1; ++i)
			out << vec[i] << ", ";
		out << vec[vec.size()-1];
	}
	out << "]";
	return out;
}

int main(int argn, char **argv) {
	//loadFields("gftables/gftable", GF, 289);
	//get_ppows(ppow, 289);
	int qvals[] = {2,3,5,7,9};
	finfield GFq,GFq2;
	for (int q: qvals) {
		cout << "q = " << q << "\n------\n";
		loadFieldPair("gftables/gftable", q, GFq, GFq2);
		vector<matrix> reprs;
		vector<gl2q_class> classes;
		get_classes(q, classes);
		for (int i=0; i<classes.size(); i++) {
			matrix M;
			get_representative(classes[i], M, GFq, GFq2);
			reprs.push_back(M);
			//printMat(M, GFq);
			//cout << endl;
		}
		vector<uint64_t> clsizes;
		get_class_sizes(q, clsizes);
		vector<uint64_t> kvals;
		get_kseries(GFq, kvals);
		vector<vector<cnum> > irrep_chars;
		get_irrep_chars(q, irrep_chars);
		
		vector<cnum> components;
		vector<uint64_t> chars, non_chars;
		for (int k: kvals) {
			if (decompose_kr2(k, GFq, reprs, irrep_chars, components)) {
				chars.push_back(k);
				//cout << k << ": " << components << endl;
			} else {
				non_chars.push_back(k);
				//cout << k << ": not a character\n";
			}
			cout << k << ":" << (k == chars[chars.size()-1] ? "y" : "n") << "\t [ ";
			for (cnum z: components) {
				if (fabs(z.im) > 1e-12 || fabs(z.re-round(z.re)) > 1e-12 || z.re < -1e-12) {
					cout << "X ";
				} else {
					cout << z << " ";
				}
			}
			cout << "]\n";			
		}
		cout << "computing with inner product...\n";
		cout << "the # of kth roots is a character for gcd(k,|GL_2(F_" << q << ")|) in " << chars << endl;
		vector<uint64_t> v;
		for (int k: chars) {
			v.push_back(k%((q+1)/2));
		}
		cout << v << endl;
		cout << "it is not a character for gcd(k, |GL_2(F_" << q << ")|) in " << non_chars << endl;
		v.clear();
		for (int k: non_chars) {
			v.push_back(k%((q+1)/2));
		}
		cout << v << endl;
	}
		/*
		Eigen::MatrixXd M;
		get_char_mat(q, M);
		//cout << M << endl;
		
		vector<uint64_t> krcts;
		
		for (int i=0; i<(q*q-1); i++) {
			for (int j=0; j<(q*q-1); j++) {
				double dot_prod = 0;
				for (int k=0; k<2*(q*q-1); k+=2) {
					dot_prod += M(k, i)*M(k, j)*clsizes[k/2];
					dot_prod += M(k+1, i)*M(k+1, j)*clsizes[k/2]; // it's times the conjugate, so we add even though we're multiplying the imaginary parts
				}
				dot_prod /= gl_size(2,GF[q]);
				if (fabs(dot_prod-(i==j)) > 1e-14) {
					cout << "the irreducible characters are not an orthonormal basis!!\n";
				}
				//cout << fix_roundoff(dot_prod) << " ";
			}
			//cout << endl;
		}
		
		vector<int> kvals;
		get_kseries(q, kvals);
		
		vector<uint64_t> chars, non_chars;
		for (int k: kvals) {
			krcts.clear();
			get_kth_root_cts(2, GF[q], k, reprs, krcts);
			uint64_t sanity_check = 0, nsol = 0;
			for (int i=0; i<clsizes.size(); i++) {
				sanity_check += krcts[i]*clsizes[i];
				nsol += krcts[i]*krcts[i]*clsizes[i];
			}
			//cout << krcts << endl;
			//cout << clsizes << endl;
			if (sanity_check != gl_size(2,GF[q])) {
				cout << "the total number of kth roots is not right\n";
				cout << sanity_check << " " << gl_size(2,GF[q]) << endl;
			}
			
			if (decompose_kr(k,GF[q],M,reprs)) {
				chars.push_back(k);
			} else {
				non_chars.push_back(k);
			}
			//cout << "actual # solutions = " << nsol << endl;
		}
		cout << "q = " << q << endl;
		cout << "the # of kth roots is a character for gcd(k,|GL_2(F_q)|) in " << chars << endl;
		vector<uint64_t> v;
		for (int k: chars) {
			v.push_back(k%((q+1)/2));
		}
		cout << v << endl;
		cout << "it is not a character for gcd(k, |GL_2(F_q)|) in " << non_chars << endl;
		v.clear();
		for (int k: non_chars) {
			v.push_back(k%((q+1)/2));
		}
		cout << v << endl;
		
		
		//	cout << "inner product-computed coeffs for k=" << k << endl;
		//	cout << components << endl;
		}
		/*for (int i=0; i<q*q-1; i++) {
			for (int j=0; j<q*q-1; j++) {
				cout << fix_roundoff(inner_product(irrep_chars[i], irrep_chars[j], clsizes, gl_size(2,GF[q])).re) << " ";
			}
			cout << endl;
		}
		
		cout << endl << endl;
		for (int i=0; i<q*q-1; i++) {
			for (int j=0; j<q*q-1; j++) {
				cout << fix_roundoff(inner_product(irrep_chars[i], irrep_chars[j], clsizes, gl_size(2,GF[q])).im) << " ";
			}
			cout << endl;
		}
		
		cout << "computing with inner product...\n";
		cout << "the # of kth roots is a character for gcd(k,|GL_2(F_q)|) in " << chars2 << endl;
		cout << "it is not a character for gcd(k, |GL_2(F_q)|) in " << non_chars2 << endl;
		
		//cout << gl_size(2,GF[q]) << endl;
		//for (int k=1; k<15; k++) {
			//decompose_kr(k, GF[q], M, reprs);
		//}
	}*/
}
