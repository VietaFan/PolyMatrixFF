/* code for constructing the character table of GL(2,F_q)
 * 
 * see http://www.math.umd.edu/~jda/characters/characters.pdf
 * and http://www.math.umd.edu/~jda/characters/gl2/ 
 */
#ifndef _GLFQCHAR_H
#define _GLFQCHAR_H
#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <algorithm>
#include <cmath>
#include <Eigen/Dense>
#define PI 3.141592653589793
using namespace std;

struct cnum {
	double re, im;
};

double l2magn(cnum x);
ostream& operator<<(ostream &out, cnum &vec);
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
struct gl2q_class {
	int q, type;
	int x_gpow, y_gpow, z_gpow;
};

/*
 * g: generator of GF[q^2]*
 * h = g^{q+1}: (a generator of GF[q]*)
 * omega = exp((2*pi*i)/(q^2-1))
 * 
 * alpha(h^n) = omega^(a*n*(q+1))
 * beta(h^n) = omega^(b*n*(q+1))
 * chi(g^n) = omega^(cn)
 * mu((x,y)) = alpha(x)*beta(y)
 * mu^w((x,y)) = alpha(y)*beta(x)
 *  
 * type 1: (q-1)(q-2)/2 representations - uses alpha,beta,mu,mu^w
 *        0 <= a < b < q-1
 * type 2: q-1 representations - uses alpha
 * 		  0 <= a < q-1
 * type 3: q-1 representations - uses alpha
 * 		  0 <= a < q-1
 * type 4: q(q-1)/2 representations
 *        0 <= c < (q^2-1)/2, c%(q+1) != 0
 * 
 * see UMD site for explanation
 */

struct gl2q_repn {
	int q, type;
	int a, b, c;
};

void get_representative(gl2q_class &cl, matrix &M, finfield &Fq, finfield &Fq2);
void get_classes(int q, vector<gl2q_class> &classes);
void get_reps(int q, vector<gl2q_repn> &reps);
void init_roots_of_unity(int n, double *real, double *imag);
double fix_roundoff(double d);

// unity_re and unity_im should be initialized to hold the (q^2-1)th roots of unity
// exp(2*pi*i*n/(q^2-1)) = unity_re[n]+i*unity_im[n] for 0 <= n < q^2-1
// see http://www.math.umd.edu/~jda/characters/gl2/ for character table
void character(gl2q_repn &rep, gl2q_class &cl, double *unity_re, 
				double *unity_im, double &real, double &imag);
void get_character(gl2q_repn &rep, vector<gl2q_class> &cl, 
		double *unity_re, double *unity_im, vector<cnum> &chvals);
void get_char_mat(int q, Eigen::MatrixXd &M);
void get_class_sizes(int q, vector<uint64_t> &class_nums);
void get_irrep_chars(int q, vector<vector<cnum> > &irrep_chars);
cnum inner_product(vector<cnum> &chi1, vector<cnum> &chi2, vector<uint64_t> &clsizes, uint64_t glsize);
#endif
