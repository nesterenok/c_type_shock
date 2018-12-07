//
// 06.03.2017. Check for errors.
// 11.09.2017. Check for errors. Collisions of CII with He were added, by assuming k_He = 0.38*k_H. 
//		The CII-e collision strength was assumed constant at T < 1000 K down to 3 K;

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <stdio.h>
#include <stdlib.h>

#include <cmath>
#include <fstream>
#include <iostream>
#include <memory.h>

#include "coll_rates_ions.h"
#include "utils.h"
#include "constants.h"

#define MAX_TEXT_LINE_WIDTH 240
#define SOURCE_NAME "coll_rates_ions.cpp"
using namespace std;

// Data files must have general structure;
ion_electron_coll_data::ion_electron_coll_data(const string fname, const energy_diagram *ion_di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, up, low;
	ifstream input;

	input.open(fname.c_str(), std::ios_base::in);
	if (!input) { 
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << fname << endl;
		exit(1);
	}
	// three comment lines:
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	input >> nb_lev >> jmax;
	
	if (nb_lev > ion_di->nb_lev)
		nb_lev = ion_di->nb_lev;

	imax = nb_lev*(nb_lev-1)/2;
	jmax++; // the number of tempratures is increased to include the point T = 0 K;

	tgrid = new double [jmax];
	
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, imax*jmax*sizeof(double));

	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, imax*jmax*sizeof(double));

	tgrid[0] = 0.;
	for (j = 1; j < jmax; j++) {
		input >> tgrid[j];
	}
	for (i = 0; i < imax; i++) 
	{
		input >> up >> low;
		for (j = 1; j < jmax; j++) 
		{
			// collision rate coefficient is assumed to be 0 at T = 0 K;
			input >> coeff[i][j];
			coeff[i][j] *= COLL_STRENGTH_TO_RATE_COEFF/( pow(tgrid[j], 0.5) *ion_di->lev_array[up].g);
		}
	}
	input.close();
	calc_coeff_deriv();

	if (verbosity) {
		cout << "  data have been read from file " << fname << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

ion_neutral_coll_data::ion_neutral_coll_data(const string fname, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, up, low;
	ifstream input;

	input.open(fname.c_str(), std::ios_base::in);
	if (!input) { 
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << fname << endl;
		exit(1);
	}
	// Three comment lines:
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	input >> nb_lev >> jmax;
	
	imax = nb_lev*(nb_lev-1)/2;
	jmax++; // the number of tempratures is increased to include the point T = 0 K;

	tgrid = new double [jmax];
	
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, imax*jmax*sizeof(double));
	
	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, imax*jmax*sizeof(double));

	tgrid[0] = 0.;
	for (j = 1; j < jmax; j++) {
		input >> tgrid[j];
	}
	for (i = 0; i < imax; i++) 
	{
		input >> up >> low;
		for (j = 1; j < jmax; j++) {
			input >> coeff[i][j];
		}
	}
	input.close();
	calc_coeff_deriv();

	if (verbosity) {
		cout << "  data have been read from file " << fname << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

//
// Data on OI collisional transitions
//

OI_electrons_berrington_data::OI_electrons_berrington_data(int verbosity)
{
	int j, nb_d;
	double b, sc, t;

	nb_d = 15; // nb of temperature points in an interval of order magnitude;
	sc = pow(10., 1./nb_d);

	nb_lev = 5; // nb of OI levels;
	imax = nb_lev*(nb_lev-1)/2;
	// 10 < T < 10000 K;
	jmax = (int) (log10(10000./10.)*nb_d) + 2;
	
	tgrid = new double [jmax];

	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, imax*jmax*sizeof(double));

	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, imax*jmax*sizeof(double));
	
	tgrid[0] = 0.;
	tgrid[1] = 10.; // temperature is in K;
	for (j = 2; j < jmax; j++) {
		tgrid[j] = sc*tgrid[j-1];
	}

	for (j = 1; j < jmax; j++)
	{
		b = COLL_STRENGTH_TO_RATE_COEFF/pow(tgrid[j], 0.5);
		t = 1.e-4*tgrid[j];

		// Original data - Berrington, J. Phys. B 21, p. 1083 (1988);
		if (tgrid[j] < 1.e+4) 
		{
			// 1 -> 0
			coeff[0][j] = b/3. *( 0.041*pow(t, 0.69) + 0.064*pow(t, 1.72) );

			// 2 -> 0
			coeff[1][j] = b	*( 0.0136*pow(t, 0.61) + 0.0186*pow(t, 1.49) );

			// 2 -> 1
			coeff[2][j] = b	*( 0.00166*pow(t, 0.71) + 0.0288*pow(t, 1.97) );
		}
		else 
		{ // The fit formulae change at temperatures > 10000 K (Pequignot, A&A 231, p. 499, 1990): 
			coeff[0][j] = b/3. *0.106*pow(t, 1.1);
			coeff[1][j] = b *0.0321 *t;
			coeff[2][j] = b *0.0283*pow(t, 1.5);
		}

		// Original data - Berrington & Burke, Planet Space Science 29, p. 377 (1981);
		// The fit is assumed to be acceptable at any temperature Pequignot, A&A 231, p. 499 (1990);
		// Pequignot (1990) provided total rate to the 3P term levels including sublevels, it was taken into account by Draine (2011); 
		// 1D_2 level has spin 2, 1S_0 level - spin 0
		// 3 -> 2
		coeff[5][j] = b/5.*0.0476*pow(t, 1.43)/(1. + 0.605*pow(t, 1.105));

		// 3 -> 0
		coeff[3][j] = coeff[5][j] *5.;
		// 3 -> 1
		coeff[4][j] = coeff[5][j] *3.;

		// 4 -> 2
		coeff[8][j] = b*0.00653 *pow(t, 1.5)/(1. + 0.8*pow(t, 1.125));

		// 4 -> 0
		coeff[6][j] = coeff[8][j] *5.;
		// 4 -> 1
		coeff[7][j] = coeff[8][j] *3.;

		// 4 -> 3
		coeff[9][j] = b *0.116*pow(t, 0.53)/(1. + 0.111*pow(t, 0.16));
	}
	calc_coeff_deriv();

	if (verbosity) {
		cout << "  data for OI-e collisions (Berrington & Burke 1982; Berrington, 1988)" << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

OI_h_abrahamsson_data::OI_h_abrahamsson_data(int verbosity)
{
	int i, j, nb_d;
	// data by Abrahamsson & Krems (2007):
	double sc, a1[8] = {4.581, -156.118, 2679.979, -78996.962, 1308323.468, -13011761.861, 71010784.971, -162826621.855},
		a2[8] = {3.297, -168.382, 1844.099, -68362.889, 1376864.737, -17964610.169, 134374927.808, -430107587.886},
		a3[8] = {3.437, 17.443, -618.761, 3757.156, -12736.468, 22785.266, -22759.228, 12668.261}; 
	// data by Krems, Jamieson, Dalgarno (2006), temperature (K), collision rate (cm3 s-1):
	double b1[14] = {100., 200., 300., 500., 700., 1000., 2000., 3000., 4000., 5000., 6000., 7000., 8000., 10000},
		b2[14] = {2.83e-13, 4.85e-13, 5.84e-13, 6.84e-13, 7.41e-13, 8.e-13, 9.38e-13, 1.03e-12, 1.09e-12, 1.12e-12, 
		1.13e-12, 1.12e-12, 1.1e-12, 1.05e-12};

	nb_d = 15;
	sc = pow(10., 1./nb_d);

	nb_lev = 4;
	imax = nb_lev*(nb_lev-1)/2;
	// the minimal temperature is 30 K; the maximal temperature - 10000 K;
	jmax = (int) (log10(10000./30.)*nb_d) + 3;
	
	tgrid = new double [jmax];

	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, imax*jmax*sizeof(double));

	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, imax*jmax*sizeof(double));

	tgrid[0] = 0.;
	tgrid[1] = 30.;
	for (j = 2; j < jmax; j++) {
		tgrid[j] = sc*tgrid[j-1];
	}

	for (j = 1; j < jmax; j++)
	{
		// The approximations given by Abrahamsson & Krems (2007) is valid at < 1000 K; only for lowest three levels;
		if (tgrid[j] < 1000.) 
		{
			// 1 -> 0
			coeff[0][j] = a1[0];
			for (i = 1; i < 8; i++) {
				coeff[0][j] += a1[i] *pow(tgrid[j], -0.75*i);
			}
			coeff[0][j] = 1.e-11*exp(coeff[0][j]) *5./3.*exp(158.265*CM_INVERSE_TO_KELVINS/tgrid[j]);

			// 2 -> 0
			coeff[1][j] = a2[0];
			for (i = 1; i < 8; i++) {
				coeff[1][j] += a2[i] *pow(tgrid[j], -0.75*i);
			}
			coeff[1][j] = 1.e-11*exp(coeff[1][j]) *5.*exp(226.977*CM_INVERSE_TO_KELVINS/tgrid[j]);

			// 2 -> 1
			coeff[2][j] = a3[0];
			for (i = 1; i < 8; i++) {
				coeff[2][j] += a3[i] *pow(tgrid[j], -0.5*i);
			}
			coeff[2][j] = 1.e-11*exp(coeff[2][j]) *3.*exp((226.977-158.265)*CM_INVERSE_TO_KELVINS/tgrid[j]);
		}
		else {
			coeff[0][j] = coeff[0][j-1];
			coeff[1][j] = coeff[1][j-1];
			coeff[2][j] = coeff[2][j-1];
		}
	}
	// The data by Krems, Jamieson, Dalgarno (2006); for transitions connecting fourth level with the lowest three;
	for (j = 1; j < jmax; j++)
	{
		i = 1;
		while (b1[i] < tgrid[j]) {
			i++;
		}
		if (tgrid[j] < 1000.) {
			coeff[5][j] = b2[i-1] + (tgrid[j] - b1[i-1])*(b2[i] - b2[i-1])/(b1[i] - b1[i-1]); 
		}
		else 
		{ // the accuracy of approximation is 1% (Krems et al., ApJ 647, p. 1531, 2006)
			sc = tgrid[j]/6000.;
			coeff[5][j] = (1.74*sc + 0.06)*1.e-12*exp(-0.47*sc)/sqrt(sc);
		}
		// 3 -> 0
		coeff[3][j] = 5./9.*coeff[5][j];
		// 3 -> 1
		coeff[4][j] = coeff[5][j] /3.;
		// 3 -> 2
		coeff[5][j] /= 9.;
	}
	calc_coeff_deriv();

	if (verbosity) {
		cout << "  data for OI-H collisions (Abrahamsson & Krems, 2007; Krems et al., 2006)" << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

OI_oh2_data::OI_oh2_data(int verbosity)
{
	int j, nb_d;
	double sc;

	nb_d = 15;
	sc = pow(10., 1./nb_d);

	nb_lev = 3;
	imax = nb_lev*(nb_lev-1)/2;
	// the minimal temperature (for which the rates are calculated here) is 10 K; the maximal temperature - 1500 K;
	jmax = (int) (log10(1500./10.)*nb_d) + 2;
	
	tgrid = new double [jmax];
	
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, imax*jmax*sizeof(double));

	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, imax*jmax*sizeof(double));

	tgrid[0] = 0.;
	tgrid[1] = 10.;
	for (j = 2; j < jmax; j++) {
		tgrid[j] = sc*tgrid[j-1];
	}
	for (j = 1; j < jmax; j++)
	{
		// Glover & Jappsen (2007) approximations (they referred Flower as private communication);
		// 1 -> 0
		coeff[0][j] = 2.7e-11*pow(tgrid[j], 0.362);
		// 2 -> 0
		coeff[1][j] = 5.49e-11*pow(tgrid[j], 0.317);
		// 2 -> 1
		coeff[2][j] = 2.74e-14*pow(tgrid[j], 1.06);
	}
	calc_coeff_deriv();

	if (verbosity) {
		cout << "  data for OI-oH2 collisions (approx. Glover & Jappsen, 2007; Jaquet et al., 1992)" << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << " (original data 20-1500 K)" << endl;
	}
}

OI_ph2_data::OI_ph2_data(int verbosity)
{
	int j, nb_d;
	double sc;

	nb_d = 15;
	sc = pow(10., 1./nb_d);

	nb_lev = 3;
	imax = nb_lev*(nb_lev-1)/2;
	// the minimal temperature (for which the rates are calculated here) is 10 K; the maximal temperature - 1500 K;
	jmax = (int) (log10(1500./10.)*nb_d) + 2;
	
	tgrid = new double [jmax];
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, imax*jmax*sizeof(double));

	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, imax*jmax*sizeof(double));

	tgrid[0] = 0.;
	tgrid[1] = 10.;
	for (j = 2; j < jmax; j++) {
		tgrid[j] = sc*tgrid[j-1];
	}
	for (j = 1; j < jmax; j++)
	{
		// Glover & Jappsen (2007) approximations;
		// 1 -> 0
		coeff[0][j] = 3.46e-11*pow(tgrid[j], 0.316);
		// 2 -> 0
		coeff[1][j] = 7.07e-11*pow(tgrid[j], 0.268);
		// 2 -> 1
		coeff[2][j] = 3.33e-15*pow(tgrid[j], 1.36);
	}
	calc_coeff_deriv();

	if (verbosity) {
		cout << "  data for OI-pH2 collisions (approx. Glover & Jappsen, 2007; Jaquet et al., 1992)" << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << " (original data 20-1500 K)" << endl;
	}
}

OI_collisions::OI_collisions(const string &path, const energy_diagram *ion_di, int verbosity)
	: collisional_transitions()
{
	if (verbosity) {
		cout << "OI collisional rate coefficients are being initializing..." << endl;
	}
	coll_data.push_back( new OI_ph2_data(verbosity) );
	coll_data.push_back( new OI_oh2_data(verbosity) );
	coll_data.push_back( new ion_neutral_coll_data(path + "coll_ions/coll_OI_he_monteiro1987.txt", verbosity) );
	coll_data.push_back( new OI_h_abrahamsson_data(verbosity) );
	nb1 = (int) coll_data.size();

	// data by Bell et al. (1998) includes three lowerest levels of OI;
	coll_data.push_back( new ion_electron_coll_data(path + "coll_ions/coll_OI_e_bell1998.txt", ion_di, verbosity) );
	
	// data by Berrington & Burke (1981), Berrington (1988), approximated by Pequignot (1990), includes 5 levels of OI;
	coll_data.push_back( new OI_electrons_berrington_data(verbosity) );
	nb2 = (int) coll_data.size();

	// data on collisions with H+ must be here;
	nb3 = (int) coll_data.size();

	max_temp = new double [nb3];
	for (int i = 0; i < nb3; i++) {
		max_temp[i] = coll_data[i]->get_max_temp();
	}
}

void OI_collisions::set_gas_param(double temp_neutrals, double temp_el, double he_conc, double ph2_conc, double oh2_conc, 
	double h_conc, double el_conc, double *&concentration, int *&indices) const
{
	collisional_transitions::set_gas_param(temp_neutrals, temp_el, he_conc, ph2_conc, oh2_conc, h_conc, el_conc, concentration, 
		indices);

	concentration[0] = ph2_conc;
	concentration[1] = oh2_conc;
	concentration[2] = he_conc;
	concentration[3] = h_conc;
}

//
// Data on CI collisional transitions
//

CI_e_johnson_data::CI_e_johnson_data(int verbosity)
{
	int i, j, nb_d;
	double sc;
	double a0[6] = {-9.25141, -7.69735, -7.4387, 444.6, 350.609, 386.186}, 
		a1[6] = {-0.773782, -1.30743, -0.57443, -227.913, -187.474, -202.192},
		a2[6] = {0.361184, 0.697638, 0.358264, 42.5952, 36.1803, 38.5049}, 
		a3[6] = {-0.0150892, -0.111338, -0.0418166, -3.4762, -3.03283, -3.19268},
		a4[6] = {-0.000656325, 0.00705277, 0.00235272, 0.105085, 0.0938138, 0.0978573};

	nb_d = 15;
	sc = pow(10., 1./nb_d);

	nb_lev = 3;
	imax = nb_lev*(nb_lev-1)/2;
	// the minimal temperature is 7.5 K; the maximal temperature - 10000 K;
	jmax = (int) (log10(10000./7.5)*nb_d) + 3;
	
	tgrid = new double [jmax];
	
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, imax*jmax*sizeof(double));

	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, imax*jmax*sizeof(double));

	tgrid[0] = 0.;
	tgrid[1] = 7.5;
	for (j = 2; j < jmax; j++) {
		tgrid[j] = sc*tgrid[j-1];
	}

	for (j = 1; j < jmax; j++)
	{
		sc = log(tgrid[j]);
		for (i = 0; i < 3; i++) 
		{
			if (tgrid[j] < 1000.)
				coeff[i][j] = exp(a0[i] + a1[i]*sc + a2[i]*sc*sc + a3[i]*sc*sc*sc + a4[i]*sc*sc*sc*sc);
			else coeff[i][j] = exp(a0[i+3] + a1[i+3]*sc + a2[i+3]*sc*sc + a3[i+3]*sc*sc*sc + a4[i+3]*sc*sc*sc*sc);
		}
		// 1 -> 0
		coeff[0][j] *= COLL_STRENGTH_TO_RATE_COEFF/( pow(tgrid[j], 0.5) *3.);
		// 2 -> 0
		coeff[1][j] *= COLL_STRENGTH_TO_RATE_COEFF/( pow(tgrid[j], 0.5) *5.);
		// 2 -> 1
		coeff[2][j] *= COLL_STRENGTH_TO_RATE_COEFF/( pow(tgrid[j], 0.5) *5.);
	}
	calc_coeff_deriv();

	if (verbosity) {
		cout << "  data for CI-e collisions (Johnson et al., 1987)" << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

// The approximations given by Abrahamsson & Krems (2007) is valid at < 1000 K; only for lowest three levels;
// see also approximations by Hollenbach, McKee, ApJ 342, p.306 (1989) based on the data by Launay, Roueff, A&A 56, p.289 (1977); 
CI_h_abrahamsson_data::CI_h_abrahamsson_data(int verbosity)
{
	int i, j, nb_d;
	double sc, a1[9] = {3.6593, 56.6023, -802.9765, 5025.1882, -17874.4255, 38343.6655, -49249.4895, 34789.3941, -10390.9809},
		a2[9] = {10.8377, -173.4153, 2024.0272, -13391.6549, 52198.5522, -124518.3586, 178182.5823, -140970.6106, 47504.5861},
		a3[9] = {15.8996, -201.303, 1533.6164, -6491.0083, 15921.9239, -22691.1632, 17334.7529, -5517.936, 0.}; 
	
	nb_d = 15;
	sc = pow(10., 1./nb_d);

	nb_lev = 3;
	imax = nb_lev*(nb_lev-1)/2;
	// the minimal temperature is 5 K; the maximal temperature - 1000 K;
	jmax = (int) (log10(1000./5.)*nb_d) + 3;
	
	tgrid = new double [jmax];
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, imax*jmax*sizeof(double));

	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, imax*jmax*sizeof(double));

	tgrid[0] = 0.;
	tgrid[1] = 5.;
	for (j = 2; j < jmax; j++) {
		tgrid[j] = sc*tgrid[j-1];
	}

	for (j = 1; j < jmax; j++)
	{
		// 1 -> 0
		coeff[0][j] = a1[0];
		for (i = 1; i < 9; i++) {
			coeff[0][j] += a1[i] *pow(tgrid[j], -0.25*i);
		}
		coeff[0][j] = 1.e-11*exp(coeff[0][j]) *exp(16.4*CM_INVERSE_TO_KELVINS/tgrid[j])/3.;

		// 2 -> 0
		coeff[1][j] = a2[0];
		for (i = 1; i < 9; i++) {
			coeff[1][j] += a2[i] *pow(tgrid[j], -0.33333333*i);
		}
		coeff[1][j] = 1.e-11*exp(coeff[1][j]) *exp(43.4*CM_INVERSE_TO_KELVINS/tgrid[j])/5.;

		// 2 -> 1
		coeff[2][j] = a3[0];
		for (i = 1; i < 9; i++) {
			coeff[2][j] += a3[i] *pow(tgrid[j], -0.25*i);
		}
		coeff[2][j] = 1.e-11*exp(coeff[2][j]) *exp((43.4-16.4)*CM_INVERSE_TO_KELVINS/tgrid[j])*3./5.;
	}
	calc_coeff_deriv();

	if (verbosity) {
		cout << "  data for CI-H collisions (Abrahamsson & Krems, 2007)" << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

CI_collisions::CI_collisions(const string &path, const energy_diagram *ion_di, int verbosity)
	: collisional_transitions()
{
	if (verbosity) {
		cout << "CI collisional rate coefficients are being initializing..." << endl;
	}
	coll_data.push_back( new ion_neutral_coll_data(path + "coll_ions/coll_CI_he_staemmler1991.txt", verbosity) );
	coll_data.push_back( new ion_neutral_coll_data(path + "coll_ions/coll_CI_ph2_schroder1991.txt", verbosity) );
	coll_data.push_back( new ion_neutral_coll_data(path + "coll_ions/coll_CI_oh2_schroder1991.txt", verbosity) );
	coll_data.push_back( new CI_h_abrahamsson_data(verbosity) );
	nb1 = (int) coll_data.size();

	coll_data.push_back( new CI_e_johnson_data(verbosity) );
	nb2 = (int) coll_data.size();
	
	// The data on collisions with H+ must be here;
	nb3 = (int) coll_data.size();

	max_temp = new double [nb3];
	for (int i = 0; i < nb3; i++) {
		max_temp[i] = coll_data[i]->get_max_temp();
	}
}

void CI_collisions::set_gas_param(double temp_neutrals, double temp_el, double he_conc, double ph2_conc, double oh2_conc, 
	double h_conc, double el_conc, double *&concentration, int *&indices) const
{
	collisional_transitions::set_gas_param(temp_neutrals, temp_el, he_conc, ph2_conc, oh2_conc, h_conc, el_conc, concentration, 
		indices);

	concentration[0] = he_conc;
	concentration[1] = ph2_conc;
	concentration[2] = oh2_conc;
	concentration[3] = h_conc;
}

//
// Data on CII collisional transitions
//

CII_e_tayal2008_data::CII_e_tayal2008_data(int verbosity)
{
	int j, nb_d;
	double sc;
	
	nb_d = 15;
	sc = pow(10., 1./nb_d);

	nb_lev = 2;
	imax = nb_lev*(nb_lev-1)/2;
	// the maximal temperature - 300000 K (Tayal, 2008);
	jmax = (int) (log10(100000./3.)*nb_d) + 2;
	
	tgrid = new double [jmax];

	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, imax*jmax*sizeof(double));

	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, imax*jmax*sizeof(double));

	tgrid[0] = 0.;
	tgrid[1] = 3.;
	for (j = 2; j < jmax; j++) {
		tgrid[j] = sc*tgrid[j-1];
	}
	
	for (j = 1; j < jmax; j++) {
		// Note: collision strength is finite at T -> 0 and collision coefficient becomes infinite;
		coeff[0][j] = (1.55 + 1.25e-4*tgrid[j])/(1. + 0.35*pow(1.e-4*tgrid[j], 1.25)) 
			*0.25*COLL_STRENGTH_TO_RATE_COEFF/pow(tgrid[j], 0.5);
	}
	calc_coeff_deriv();

	if (verbosity) {
		cout << "  data for CII-e collisions (Tayal, 2008)" << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

CII_h2_flower1977::CII_h2_flower1977(bool is_para_h2, int verbosity)
{
	int j, nb_d;
	double sc;
	
	nb_d = 15;
	sc = pow(10., 1./nb_d);

	nb_lev = 2;
	imax = nb_lev*(nb_lev-1)/2;
	// The maximal temperature - 250 K (Flower & Launay, 1977);
	jmax = (int) (log10(250./10.)*nb_d) + 3;
	
	tgrid = new double [jmax];
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, imax*jmax*sizeof(double));

	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, imax*jmax*sizeof(double));

	tgrid[0] = 0.;
	tgrid[1] = 10.;
	for (j = 2; j < jmax; j++) {
		tgrid[j] = sc*tgrid[j-1];
	}
	
	for (j = 1; j < jmax; j++)
	{
		if (is_para_h2)	
			coeff[0][j] = 4.25e-10*pow(0.01*tgrid[j], 0.124 - 0.018*log(0.01*tgrid[j]) );
		else coeff[0][j] = 5.14e-10*pow(0.01*tgrid[j], 0.095 + 0.023*log(0.01*tgrid[j]) );
	}
	calc_coeff_deriv();

	if (verbosity) {
		cout << "  data for CII-H2 collisions (Flower & Launay, 1977)" << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}


CII_h2_lique2013::CII_h2_lique2013(bool is_para_h2, int verbosity)
{
	const int nb = 7;
	int j;
	double a, at[nb] = {10., 20., 50., 100., 200., 300., 500.},
		aj0[nb] = {4.36e-10, 4.53e-10, 4.63e-10, 4.67e-10, 4.59e-10, 4.48e-10, 4.36e-10}, // H2 j = 0 -> j = 0
		aj1[nb] = {5.29e-10, 5.33e-10, 5.37e-10, 5.45e-10, 5.62e-10, 5.71e-10, 5.79e-10}, // H2 j = 1 -> j = 1
		aj2[nb] = {4.43e-10, 4.48e-10, 4.58e-10, 4.72e-10, 4.99e-10, 5.20e-10, 5.43e-10}, // H2 j = 2 -> j = 2
		aj3[nb] = {0.73e-10, 0.73e-10, 0.74e-10, 0.74e-10, 0.74e-10, 0.74e-10, 0.74e-10}; // H2 j = 2 -> j = 0

	nb_lev = 2;
	imax = nb_lev*(nb_lev-1)/2;
	jmax = nb+1;
	
	tgrid = new double [jmax];
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, imax*jmax*sizeof(double));

	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, imax*jmax*sizeof(double));

	tgrid[0] = 0.;
	for (j = 1; j < jmax; j++) {
		tgrid[j] = at[j-1];
	}
	
	for (j = 1; j < jmax; j++)
	{
		if (is_para_h2)	{
			a = 5.*exp(-354.35*CM_INVERSE_TO_KELVINS/at[j-1]); // local thermodynamic population is considered;
			coeff[0][j] = (aj0[j-1] + (aj2[j-1] + aj3[j-1])*a)/(1. + a);
		}
		else {
			coeff[0][j] = aj1[j-1];
		}
	}
	calc_coeff_deriv();

	if (verbosity) {
		cout << "  data for CII-H2 collisions (Lique et al. 2013)" << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}


CII_collisions::CII_collisions(const string &path, const energy_diagram *ion_di, int verbosity)
	: collisional_transitions()
{
	bool is_para_h2;

	if (verbosity) {
		cout << "CII collisional rate coefficients are being initializing..." << endl;
	}
	coll_data.push_back( new CII_h2_lique2013(is_para_h2 = true, verbosity) );
	coll_data.push_back( new CII_h2_lique2013(is_para_h2 = false, verbosity) );
	coll_data.push_back( new ion_neutral_coll_data(path + "coll_ions/coll_CII_h_barinovs2005.txt", verbosity) );
	nb1 = (int) coll_data.size();

	coll_data.push_back( new CII_e_tayal2008_data(verbosity) );
	nb2 = (int) coll_data.size();
	
	// The data on collisions with H+ must be here;
	nb3 = (int) coll_data.size();

	max_temp = new double [nb3];
	for (int i = 0; i < nb3; i++) {
		max_temp[i] = coll_data[i]->get_max_temp();
	}
}

void CII_collisions::set_gas_param(double temp_neutrals, double temp_el, double he_conc, double ph2_conc, double oh2_conc, double h_conc,
	double el_conc, double *&concentration, int *&indices) const
{
	collisional_transitions::set_gas_param(temp_neutrals, temp_el, he_conc, ph2_conc, oh2_conc, h_conc, el_conc, concentration, 
		indices);

	concentration[0] = ph2_conc;
	concentration[1] = oh2_conc;
	// Draine B.T., "Interstellar medium" (2011), p.501, collisions with He can be estimated by 0.38* H rate;
	concentration[2] = h_conc + 0.38*he_conc;
}
