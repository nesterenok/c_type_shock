
#include <stdlib.h>
#include <iostream>
#include <math.h>
#include <limits>
#include "special_functions.h"

using namespace std;

// Error function; 
// Press W.H. et al. Numerical recipes in C, Cambridge University Press, 3 ed., 2007;
double err_func::f(double x)
{
	if (x >= 0.) return 1. - fccheb(x);
	else return fccheb(-x) - 1.;
}

double err_func::fc(double x)
{
	if (x >= 0.) return fccheb(x);
	else return 2. - fccheb(-x);
}

// Evaluate function using stored Chebyshev coefficients; 
double err_func::fccheb(double z)
{
	int j;
	double t, ty, tmp, d = 0., dd = 0.;
	if (z < 0.) {
		cout << endl << "Math error: erfccheb requires nonnegative arguments;";
		exit(1);
	}
	t = 2./(2. + z);
	ty = 4.*t - 2.;
	for (j = ncof-1; j > 0; j--)
	{
		tmp = d;
		d = ty*d - dd + cof[j];
		dd = tmp;
	}
	return t*exp(-z*z + 0.5*(cof[0] + ty*d) - dd);
}

// Inverse of complementary error function. Returns x such that erfc(x) = p for argument p in [0,2];
double err_func::inverfc(double p)
{
	double x, err, t, pp;
	if (p >= 2.) return -100.;
	if (p <= 0.) return 100.;

	pp = (p < 1.) ? p : 2. - p;
	t = sqrt(-2.*log(pp/2.));
	x = -0.70711*((2.30753 + t*0.27061)/(1. + t*(0.99229 + t*0.04481)) - t);

	for (int j = 0; j < 2; j++) {
		err = fc(x) - pp;
		x += err/(1.12837916709551257 *exp(-x*x) - x*err);
	}
	return (p < 1. ? x : -x);
}

const double err_func::cof[28] = { -1.3026537197817094, 6.4196979235649026e-1, 1.9476473204185836e-2, -9.561514786808631e-3,
	-9.46595344482036e-4, 3.66839497852761e-4, 4.2523324806907e-5, -2.0278578112534e-5, -1.624290004647e-6, 1.303655835580e-6, 
	1.5626441722e-8, -8.5238095915e-8, 6.529054439e-9, 5.059343495e-9, -9.91364156e-10, -2.27365122e-10, 9.6467911e-11, 
	2.394038e-12, -6.886027e-12, 8.94487e-13, 3.13092e-13, -1.12708e-13, 3.81e-16, 7.106e-15, -1.523e-15, -9.4e-17, 1.21e-16,
	-2.8e-17};

// Exponential integral E_n(x);
//  Press W.H. et al. Numerical recipes in C, Cambridge University Press, 3 ed., 2007;
double expint(const int n, const double x)
{
	static const int MAXIT = 100;
	static const double EULER = 0.577215664901533, 
		EPS = numeric_limits<double>::epsilon(), // desired relative error;
		BIG = numeric_limits<double>::max()*EPS; // number near the largest representable floating-point number;

	int i, ii, nm1 = n-1;
	double a, b, c, d, del, fact, h, psi, ans;
	
	if (n < 0 || x < 0. || (x == 0. && (n == 0 || n == 1))) {
		cout << endl << "Math error: Bad arguments in function expint();";
		exit(1);
	}
	if (n == 0) ans = exp(-x)/x;
	else {
		if (x == 0.) ans = 1./nm1;
		else {
			if (x > 1.) 
			{
				b = x + n;
				c = BIG;
				d = 1./b;
				h = d;
				for (i = 1; i <= MAXIT; i++) 
				{
					a = -i*(nm1 + i);
					b += 2.;
					d = 1./(a*d + b);  // denominators cannot be zero;
					c = b + a/c;
					del = c*d;
					h *= del;

					if (fabs(del - 1.) <= EPS) {	
						ans = h*exp(-x);
						return ans;
					}
					if (i == MAXIT) 
						ans = h*exp(-x);
				}
				cout << endl << "Math warning: continued fraction failed in function expint();" << endl;
			}
			else
			{
				ans = (nm1 != 0 ? 1./nm1 : -log(x) - EULER);
				fact = 1.;
				for (i = 1; i <= MAXIT; i++) 
				{
					fact *= -x/i;
					if (i != nm1) del = -fact/(i - nm1);
					else 
					{
						psi = -EULER;
						for (ii = 1; ii <= nm1; ii++) psi += 1./ii;
						del = fact *(-log(x) + psi);
					}
					ans += del;
					if (fabs(del) < fabs(ans) *EPS) return ans;
				}
				cout << endl << "Math warning: series failed in function expint();" << endl;
			}
		}
	}
	return ans;
}

// Logarithm of the gamma function.
// The code was downloaded from https://www.astro.umd.edu/~ricotti/NEWWEB/teaching/ASTR415/NR3/code/
// Press W.H. et al. Numerical recipes in C, Cambridge University Press, 3 ed., 2007;
double gammln(const double xx) 
{
	int j;
	double x, tmp, y, ser;
	static const double cof[14] = {57.1562356658629235, -59.5979603554754912, 14.1360979747417471, 
		-0.491913816097620199, .339946499848118887e-4, .465236289270485756e-4,-.983744753048795646e-4,
		.158088703224912494e-3, -.210264441724104883e-3, .217439618115212643e-3, -.164318106536763890e-3,
	.844182239838527433e-4, -.261908384015814087e-4, .368991826595316234e-5};
	
	if (xx <= 0) {
		throw("bad arg in gammln");
	}
	y = x = xx;
	tmp = x + 5.24218750000000000;
	tmp = (x + 0.5)*log(tmp) - tmp;
	ser = 0.999999999999997092;
	
	for (j = 0; j < 14; j++) ser += cof[j]/++y;
	return tmp + log(2.5066282746310005*ser/x);
}
