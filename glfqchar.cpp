/* code for constructing the character table of GL(2,F_q)
 * 
 * see http://www.math.umd.edu/~jda/characters/characters.pdf
 * and http://www.math.umd.edu/~jda/characters/gl2/ 
 */

#include <cstdlib>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "finfield.h"
#include "glfqchar.h"
#define PI 3.141592653589793
using namespace std;

double l2magn(cnum x) {
	return sqrt(x.re*x.re+x.im*x.im);
}
ostream& operator<<(ostream &out, cnum &vec) {
	if (fabs(vec.re) < 1e-10) {
		if (fabs(vec.im) < 1e-10) {
			out << "0";
			return out;
		} else {
			out << vec.im << "i";
			return out;
		}
	} else if (fabs(vec.im) < 1e-10) {
		out << vec.re;
		return out;
	} else {
		out << vec.re << "+" << vec.im << "i";
	}
	return out;
}
/*
 * g: generator of GF[q^2]*
 * 
 * type 1: 0 <= x_gpow < q-1, x = g**(x_gpow*(q+1))
 * type 2: 0 <= x_gpow < q-1, x = g**(x_gpow*(q+1))
 * type 3: 0 <= x_gpow < y_gpow < q-1
 *          x = g**(x_gpow*(q+1)), y = g**(y_gpow*(q+1))
 * type 4: 0 <= z_gpow < (q^2-1)/2, z%(q+1) != 0, z = g**(z_gpow)
 * 		    (it's isomorphic to its conjugate, so we'll assume 
 * 		     the exponent is in the lower half. g**(q^2-1) = 1,
 * 	         so we have g**z = g**(q^2-1-z))
 * 		    also z_gpow = (q+1)n is handled earlier
 * 
 * type 1: (q-1) conjugacy classes
 * type 2: (q-1) conjugacy classes
 * type 3: (q-1)(q-2)/2 conjugacy classes
 * type 4: q(q-1)/2 conjugacy classes
 * -----------------------------------
 * q^2-1 conjugacy classes in total
 *
 * see UMD site for explanation
 */
 
 // sets M to be a matrix in GL(2,F) that is in the conjugacy class cl
void get_representative(gl2q_class &cl, matrix &M, finfield &Fq, finfield &Fq2) {
	switch (cl.type) {
		case 1:
			zeros(M,2);
			M.mat[0][0] = Fq.genpow[cl.x_gpow];
			M.mat[1][1] = Fq.genpow[cl.x_gpow];
			break;
		case 2:
			zeros(M,2);
			M.mat[0][0] = Fq.genpow[cl.x_gpow];
			M.mat[1][1] = Fq.genpow[cl.x_gpow];
			M.mat[0][1] = 1;
			break;
		case 3:
			zeros(M, 2);
			M.mat[0][0] = Fq.genpow[cl.x_gpow];
			M.mat[1][1] = Fq.genpow[cl.y_gpow];
			break;
		case 4:
			zeros(M, 2);
			int q = Fq.q;
			// special case: to be able to take square roots as in the UMD construction
			// we need the field to be of odd characteristic
			//
			// There's just one representative of case 4 for GL(2,2); we give it explicitly here.
			if (q == 2) {
				M.mat[0][1] = 1;
				M.mat[1][0] = 1;
				M.mat[1][1] = 1;
				break;
			}
			int Delta = Fq.gen;
			int delta = Fq2.genpow[(q+1)/2];
			int z;
			vector<int> qpows, q2pows; // contains the elements of Fq in Fq2
			qpows.push_back(0);
			q2pows.push_back(0);
			// we identify Fq.gen with Fq2.gen^(q+1)
			for (int x=0; x<q-1; x++) {
				qpows.push_back(Fq.genpow[x]);
				q2pows.push_back(Fq2.genpow[x*(q+1)]);
			}
			for (int x=0; x<q; x++) {
				for (int y=0; y<q; y++) {
					z = Fq2.add[q2pows[x]][Fq2.mul[q2pows[y]][delta]];
					if (z == Fq2.genpow[cl.z_gpow]) {
						M.mat[0][0] = qpows[x];
						M.mat[0][1] = Fq.mul[qpows[y]][Delta];
						M.mat[1][0] = qpows[y];
						M.mat[1][1] = qpows[x];
						goto complete;
					}
				}
			}
			complete:						
			break;
	}
}

void get_classes(int q, vector<gl2q_class> &classes) {
	for (int i=0; i<q-1; i++) {
		gl2q_class cl;
		cl.type = 1;
		cl.q = q;
		cl.x_gpow = i;
		classes.push_back(cl);
	}
	for (int i=0; i<q-1; i++) {
		gl2q_class cl;
		cl.type = 2;
		cl.q = q;
		cl.x_gpow = i;
		classes.push_back(cl);
	}
	for (int i=0; i<q-1; i++) {
		for (int j=i+1; j<q-1; j++) {
			gl2q_class cl;
			cl.type = 3;
			cl.q = q;
			cl.x_gpow = i;
			cl.y_gpow = j;
			classes.push_back(cl);
		}
	}
	vector<bool> used;
	for (int i=0; i<q*q-1; i++) {
		used.push_back(0);
	}
	for (int i=0; i<q*q-1; i++) {
		if (i%(q+1) == 0 || used[i]) continue;
		gl2q_class cl;
		cl.type = 4;
		cl.q = q;
		cl.z_gpow = i;
		classes.push_back(cl);
		used[i] = 1;
		used[(i*q)%(q*q-1)] = 1;
	}
	//cout << classes.size() << endl;
}

void get_reps(int q, vector<gl2q_repn> &reps) {
	for (int i=0; i<q-1; i++) {
		for (int j=i+1; j<q-1; j++) {
			gl2q_repn rep;
			rep.type = 1;
			rep.q = q;
			rep.a = i;
			rep.b = j;
			reps.push_back(rep);
		}
	}
	for (int i=0; i<q-1; i++) {
		gl2q_repn rep;
		rep.type = 2;
		rep.q = q;
		rep.a = i;
		reps.push_back(rep);
	}
	for (int i=0; i<q-1; i++) {
		gl2q_repn rep;
		rep.type = 3;
		rep.q = q;
		rep.a = i;
		reps.push_back(rep);
	}
	vector<bool> used;
	for (int i=0; i<q*q-1; i++) {
		used.push_back(0);
	}
	for (int i=0; i<q*q-1; i++) {
		if (i%(q+1) == 0 || used[i]) continue;
		gl2q_repn rep;
		rep.type = 4;
		rep.q = q;
		rep.c = i;
		reps.push_back(rep);
		used[i] = 1;
		used[(i*q)%(q*q-1)] = 1;
	}
}

void init_roots_of_unity(int n, double *real, double *imag) {
	double omega = 2*PI/n;
	for (int t=0; t<n; t++) {
		real[t] = cos(omega*t);
		imag[t] = sin(omega*t);
	}
}

// unity_re and unity_im should be initialized to hold the (q^2-1)th roots of unity
// exp(2*pi*i*n/(q^2-1)) = unity_re[n]+i*unity_im[n] for 0 <= n < q^2-1
// see http://www.math.umd.edu/~jda/characters/gl2/ for character table
void character(gl2q_repn &rep, gl2q_class &cl, double *unity_re, 
		double *unity_im, double &real, double &imag) {
	int q = rep.q;
	int mod = q*q-1;
	int ax_res,bx_res,ay_res,by_res;
	switch (rep.type*10+cl.type) {
		case 11:
			ax_res = (cl.x_gpow*rep.a*(q+1))%mod;
			bx_res = (cl.x_gpow*rep.b*(q+1))%mod;
			real = (q+1)*( unity_re[ax_res]*unity_re[bx_res]-unity_im[ax_res]*unity_im[bx_res] );
			imag = (q+1)*( unity_re[ax_res]*unity_im[bx_res]+unity_im[ax_res]*unity_re[bx_res] );
			break;
		case 12:
			ax_res = (cl.x_gpow*rep.a*(q+1))%mod;
			bx_res = (cl.x_gpow*rep.b*(q+1))%mod;
			real = ( unity_re[ax_res]*unity_re[bx_res]-unity_im[ax_res]*unity_im[bx_res] );
			imag = ( unity_re[ax_res]*unity_im[bx_res]+unity_im[ax_res]*unity_re[bx_res] );
			break;
		case 13:
			ax_res = (cl.x_gpow*rep.a*(q+1))%mod;
			bx_res = (cl.x_gpow*rep.b*(q+1))%mod;
			ay_res = (cl.y_gpow*rep.a*(q+1))%mod;
			by_res = (cl.y_gpow*rep.b*(q+1))%mod;
			real = unity_re[ax_res]*unity_re[by_res]-unity_im[ax_res]*unity_im[by_res]
				  +unity_re[ay_res]*unity_re[bx_res]-unity_im[ay_res]*unity_im[bx_res];
			imag = unity_re[ax_res]*unity_im[by_res]+unity_im[ax_res]*unity_re[by_res]
				  +unity_re[ay_res]*unity_im[bx_res]+unity_im[ay_res]*unity_re[bx_res];
			break;			
		case 14:
			real = 0.0;
			imag = 0.0;
			break;
		case 21:
			real = q*unity_re[(2*cl.x_gpow*rep.a*(q+1))%mod];
			imag = q*unity_im[(2*cl.x_gpow*rep.a*(q+1))%mod];
			break;
		case 22:
			real = 0.0;
			imag = 0.0;
			break;
		case 23:
			real = unity_re[((cl.x_gpow+cl.y_gpow)*rep.a*(q+1))%mod];
			imag = unity_im[((cl.x_gpow+cl.y_gpow)*rep.a*(q+1))%mod];
			break;
		case 24:
			real = -unity_re[(cl.z_gpow*rep.a*(q+1))%mod];
			imag = -unity_im[(cl.z_gpow*rep.a*(q+1))%mod];
			break;
		case 31:
			real = unity_re[(2*cl.x_gpow*rep.a*(q+1))%mod];
			imag = unity_im[(2*cl.x_gpow*rep.a*(q+1))%mod];
			break;
		case 32:
			real = unity_re[(2*cl.x_gpow*rep.a*(q+1))%mod];
			imag = unity_im[(2*cl.x_gpow*rep.a*(q+1))%mod];
			break;
		case 33:
			real = unity_re[((cl.x_gpow+cl.y_gpow)*rep.a*(q+1))%mod];
			imag = unity_im[((cl.x_gpow+cl.y_gpow)*rep.a*(q+1))%mod];
			break;
		case 34:
			real = unity_re[(cl.z_gpow*rep.a*(q+1))%mod];
			imag = unity_im[(cl.z_gpow*rep.a*(q+1))%mod];
			break;
		case 41:
			real = (q-1)*unity_re[((q+1)*cl.x_gpow*rep.c)%mod];
			imag = (q-1)*unity_im[((q+1)*cl.x_gpow*rep.c)%mod];
			break;
		case 42:
			real = -unity_re[((q+1)*cl.x_gpow*rep.c)%mod];
			imag = -unity_im[((q+1)*cl.x_gpow*rep.c)%mod];
			break;
		case 43:
			real = 0.0;
			imag = 0.0;
			break;
		case 44:
			real = -unity_re[(cl.z_gpow*rep.c)%mod]-unity_re[(q*cl.z_gpow*rep.c)%mod];
			imag = -unity_im[(cl.z_gpow*rep.c)%mod]-unity_im[(q*cl.z_gpow*rep.c)%mod];
			break;
	}
}

void get_character(gl2q_repn &rep, vector<gl2q_class> &cl, 
		double *unity_re, double *unity_im, vector<cnum> &chvals) {
	for (int i=0; i<cl.size(); i++) {
		character(rep, cl[i], unity_re, unity_im, chvals[i].re, chvals[i].im);
	}
}

void get_irrep_chars(int q, vector<vector<cnum> > &irrep_chars) {
	vector<gl2q_class> classes;
	vector<gl2q_repn> reps;
	get_classes(q, classes);
	get_reps(q, reps);
	double *unity_re = new double[q*q-1];
	double *unity_im = new double[q*q-1];
	init_roots_of_unity(q*q-1, unity_re, unity_im);
	for (gl2q_repn rep: reps) {
		vector<cnum> chvals(q*q-1);
		get_character(rep, classes, unity_re, unity_im, chvals);
		irrep_chars.push_back(chvals);
	}
	delete[] unity_re;
	delete[] unity_im;
}

cnum inner_product(vector<cnum> &chi1, vector<cnum> &chi2, vector<uint64_t> &clsizes, uint64_t glsize) {
	cnum ans;
	ans.re = 0;
	ans.im = 0;
	for (int i=0; i<clsizes.size(); i++) {
		ans.re += clsizes[i]*(chi1[i].re*chi2[i].re+chi1[i].im*chi2[i].im); // we're adding b/c it's the conjugate
		ans.im += clsizes[i]*(chi1[i].im*chi2[i].re-chi1[i].re*chi2[i].im);
	}
	ans.re /= glsize;
	ans.im /= glsize;
	fix_roundoff(ans.re);
	fix_roundoff(ans.im);
	return ans;
}

double fix_roundoff(double d) {
	if (fabs(d) < 1e-12) {
		return 0;
	}
	double r = round(d);
	if (fabs(d-r) < 1e-12) {
		return r;
	} else {
		return d;
	}
}
// returns a (2*(q^2-1)) x (q^2-1) MatrixXd of doubles
// position (2*i,j) has the real part of the character of the jth irrep on the ith conjugacy class
// position (2*i+1,j) has the imaginary part of the character of the jth irrep on the ith conjugacy class
void get_char_mat(int q, Eigen::MatrixXd &M) {
	int Q = q*q-1;
	M = Eigen::MatrixXd(2*Q, Q);
	double *unity_re = new double[Q],
		   *unity_im = new double[Q];
	init_roots_of_unity(Q, unity_re, unity_im);
	vector<gl2q_class> conjcl;
	vector<gl2q_repn> irreps;
	get_classes(q, conjcl);
	get_reps(q, irreps);
	double re,im;
	for (int i=0; i<Q; i++) {
		for (int j=0; j<Q; j++) {
			character(irreps[j], conjcl[i], unity_re, unity_im, re, im);
			M(2*i,j) = fix_roundoff(re);
			M(2*i+1,j) = fix_roundoff(im);
		}
	}
	delete[] unity_re;
	delete[] unity_im;
}

void get_class_sizes(int q, vector<uint64_t> &class_nums) {
	for (int i=0; i<q-1; i++) {
		class_nums.push_back(1);
	}
	for (int i=0; i<q-1; i++) {
		class_nums.push_back(q*q-1);
	}
	for (int i=0; i<(q-1)*(q-2)/2; i++) {
		class_nums.push_back(q*q+q);
	}
	for (int i=0; i<q*(q-1)/2; i++) {
		class_nums.push_back(q*q-q);
	}
}
