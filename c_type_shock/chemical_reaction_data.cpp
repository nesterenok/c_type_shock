
//
// 14.03.2017. Check for errors.
// 18.09.2017. Check for errors.
// 31.01.2018. Check for errors.

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include "utils.h"
#include "integration.h"
#include "special_functions.h"
#include "constants.h"
#include "parameters.h"
#include "photoelectric_emission.h"
#include "chemical_reaction_data.h"

#include <stdio.h>
#include <stdlib.h> 

#include <cfloat>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sstream>

#define MAX_TEXT_LINE_WIDTH 320 // long lines
#define SOURCE_NAME "chemical_reaction_data.cpp"
using namespace std;

//
// Classes for reaction cross sections
//
process_cross_section::process_cross_section(double m1, double m2) : mass1(m1), mass2(m2), en_min(0.)
{
	reduced_mass = mass1*mass2/(mass1 + mass2);
	mass_sum = mass1 + mass2;
}

process_cross_section::process_cross_section(double m) : mass1(m), mass2(1.e+99), en_min(0.) 
{
	reduced_mass = mass1;
	mass_sum = mass2; 
}

chem_reaction_cross_section_1::chem_reaction_cross_section_1(double m1, double m2, double a, double b, double c)
	: process_cross_section(m1, m2)
{
	p = b + 0.5;
	en_min = c*BOLTZMANN_CONSTANT; // must be in erg;
	sigma = a* 0.25* sqrt(2.*M_PI*reduced_mass) *exp(-gammln(p+1)) *pow(BOLTZMANN_CONSTANT, -b);
}

// given energy is in erg:
double chem_reaction_cross_section_1::get(double en) const
{
	if (en > en_min)
		return sigma*pow(en - en_min, p); // in cm2*erg
	else return 0.;
}

chem_reaction_cross_section_2::chem_reaction_cross_section_2(double m1, double m2, double ioniz_potential, double oscill_strength)
	: process_cross_section(m1, m2)
{
	double q;
	
	q = 13.6*EV_TO_ERGS/ioniz_potential;
	en_min = 2.*ioniz_potential; // in erg;
	// for the speed up of the calculations the cs is divided by en_min; 
	sigma = 7.67e-19*q*q*oscill_strength/en_min; // in cm2/erg 
}

double chem_reaction_cross_section_2::get(double en) const
{
	if (en > en_min)
		return sigma*en*(en - en_min); // in cm2*erg
	else return 0.;
}

chem_reaction_cross_section_table::chem_reaction_cross_section_table(const std::string &data_path, const std::string &name, double m1, double m2)
	: process_cross_section(m1, m2)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i;
	double energy, cs;
	
	string file_name;
	ifstream input;

	file_name = data_path + name;
	input.open(file_name.c_str(), ios_base::in);

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with cross section data " << file_name << endl;
		exit(1);
	}
	// comment lines are read:
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	input >> nb_cs;

	en_arr = new double [nb_cs];
	cs_arr = new double [nb_cs];

	for (i = 0; i < nb_cs; i++) 
	{
		input >> energy >> cs;
		// the energy data in the file have dimension eV, conversion to erg:
		en_arr[i] = energy*EV_TO_ERGS;
		// the cross section data in the file have dimension 1.e-16 cm2:
		cs_arr[i] = cs*1.e-16;
	}
	input.close();
}

chem_reaction_cross_section_table::~chem_reaction_cross_section_table()
{
	delete [] en_arr;
	delete [] cs_arr;
}

double chem_reaction_cross_section_table::get(double energy) const
{
	int l;
	double answ;
	// if values in question are out of range of the array with a dimension dim, returned index < 0 or index >= dim-1
	locate_index(en_arr, nb_cs, energy, l);

	if (l < 0)
		answ = 0.;
	else if (l >= nb_cs-1) 
		answ = cs_arr[nb_cs-1];
	else answ = cs_arr[l] + (energy - en_arr[l]) *(cs_arr[l+1] - cs_arr[l])/(en_arr[l+1] - en_arr[l]);

	return answ*energy; // cm2*erg
}

// Sputtering yield
sputtering_yield::sputtering_yield(double m_p, double m_t, double b_en) 
	: process_cross_section(m_p), bind_en(b_en*BOLTZMANN_CONSTANT)
{
	double a, rm, sm;
	
	sm = m_p + m_t;
	rm = m_p*m_t/sm;

	a = 4.*rm/sm;       // efficiency factor is set to be equal 1.;
	en_min = bind_en/a; // here binding energy is in erg;
	
	if (en_min < 4.*bind_en)
		en_min = 4.*bind_en;

	norm = 2.*8.3e-4*a*a/(bind_en*bind_en); // the factor 2. is due to angle averaging (Draine, Salpeter, 1979);
	p = a/(30.*bind_en);
}

// energy is in erg,
double sputtering_yield::get(double energy) const
{
	return norm*energy*(energy - en_min)*(energy - en_min)/(1. + pow(p*energy, 1.333333));
}

//
// Class with rate coefficient data 
//
reaction_rate_data::reaction_rate_data() : nb_v1(0), nb_v2(0), var1(0), var2(0), f_arr(0), int_err(1.e-5)
{;}

reaction_rate_data::reaction_rate_data(const reaction_rate_data & obj)
{
	int_err = obj.int_err;
	nb_v1 = obj.nb_v1;
	nb_v2 = obj.nb_v2;
	
	if (nb_v1 > 0 && nb_v2 > 0) 
	{
		var1 = new double [nb_v1];
		var2 = new double [nb_v2];
		
		memcpy(var1, obj.var1, nb_v1*sizeof(double));
		memcpy(var2, obj.var2, nb_v2*sizeof(double));
		
		f_arr = alloc_2d_array<double>(nb_v1, nb_v2);
		memcpy(*f_arr, *obj.f_arr, nb_v1*nb_v2*sizeof(double));
	}
	else {
		var1 = var2 = 0;
		f_arr = 0;
	}
}

void reaction_rate_data::delete_data()
{
	delete [] var1;
	delete [] var2;
	if (f_arr != 0) 
		free_2d_array(f_arr);
}

void reaction_rate_data::calc_data(process_cross_section* cross_section)
{
	int i, j, nb;
	double a, at, s, b, en_min;
	// temperature in K,
	double tmin = 3., tmax = 3.e+5, smin = 0.01, smax = 50;

	ch_func1 f1;
	ch_func2 f2;

	delete_data();

	nb = 30;
	a = pow(10., 1./nb);
	
	nb_v1 = (int) (nb*log10(tmax/tmin)) + 1;
	nb_v2 = (int) (nb*log10(smax/smin)) + 2; 

	var1 = new double [nb_v1];
	var2 = new double [nb_v2];

	f_arr = alloc_2d_array<double>(nb_v1, nb_v2);
	memset(*f_arr, 0, nb_v1*nb_v2*sizeof(double));

	var1[0] = tmin;
	for (i = 1; i < nb_v1; i++) {
		var1[i] = var1[i-1]*a;
	}
	
	var2[0] = 0.;
	var2[1] = smin;
	for (i = 2; i < nb_v2; i++) {
		var2[i] = var2[i-1]*a;
	}

	en_min = cross_section->en_min;
	f1.cs = f2.cs = cross_section;
	
	for (i = 0; i < nb_v1; i++) 
	{
		// the lower limit of the integration:
		at = en_min/(BOLTZMANN_CONSTANT*var1[i]);
		f1.p1 = f2.p1 = var1[i]*BOLTZMANN_CONSTANT;
		
		// The case s = 0 (no drift) needs special treating,
		// improper integration, the integrand function may go to infinite value at the lower boundary:
		f_arr[i][0] = sqrt(8./(M_PI *f2.p1 *cross_section->reduced_mass)) 
			*qromo<ch_func2>(f2, at, at + 15., int_err);
	
		at = sqrt(at);
		for (j = 1; j < nb_v2; j++) 
		{
			s = var2[j];
			if (s > at + 6.) {
				a = s - 6.;
				b = s + 6.;
			}
			else if (s > at) {
				a = at;
				b = s + 6.;
			}
			else {
				a = at;
				b = at + 6.;
			}

			f1.p2 = s;
			f_arr[i][j] = qromo<ch_func1>(f1, a, b , int_err);

			if (s < 6.) // may be higher upper value of s?
			{
				f1.p2 = -s;
				f_arr[i][j] -= qromo<ch_func1>(f1, at, at + 6., int_err);
			}	
			f_arr[i][j] *= M_SQRT2/(sqrt(M_PI *f1.p1 *cross_section->reduced_mass) *s);
		}
	}
}

double reaction_rate_data::get(double v1, double v2) const
{
	int k, l;
	double t, u;
	
	// if values in question are out of range of the array with a dimension dim, returned index < 0 or index >= dim-1
	locate_index(var1, nb_v1, v1, k);
	locate_index(var2, nb_v2, v2, l);

	if (k < 0) {
		k = 0;
		t = 0.;
	}
	else if (k >= nb_v1-1) {
		k = nb_v1-2;
		t = 1.;
	}
	else t = (v1 - var1[k])/(var1[k+1] - var1[k]);

	if (l < 0) {
		l = 0;
		u = 0.;
	}
	else if (l >= nb_v2-1) {
		l = nb_v2-2;
		u = 1.;
	}
	else u = (v2 - var2[l])/(var2[l+1] - var2[l]);
	
	return f_arr[k][l]*(1.-t)*(1.-u) + f_arr[k+1][l]*t*(1.-u) + f_arr[k][l+1]*(1.-t)*u + f_arr[k+1][l+1]*u*t;
}

void reaction_rate_data::print(std::string rname)
{
	int i, j;
	string fname;
	ofstream output;

	fname = rname + ".txt";
	output.open(fname.c_str());
	output << scientific;
	output.precision(2);
	
	output << "<Sputt yield *velocity>;" << endl 
		<< "Row - reduced temperature [K], column - drift velocity normalized on reduced thermal speed;" << endl;
	output << setw(12) << " ";
	for (j = 0; j < nb_v2; j++) {
		output << left << setw(12) << var2[j];
	}
	output << endl;

	for (i = 0; i < nb_v1; i++) {
		output << left << setw(12) << var1[i];
		for (j = 0; j < nb_v2; j++) {
			output << left << setw(12) << f_arr[i][j];
		}
		output << endl;
	}
	output.close();
}

//
// The class with the data on accretion rate coefficients
//

accretion_rate_functions::accretion_rate_functions()
{
	int i, j, nb;
	double t, sc, min_rtemp, max_rtemp;
	
	constant = BOLTZMANN_CONSTANT/(ELECTRON_CHARGE*ELECTRON_CHARGE);

	nb = 10;
	sc = pow(10, 1./nb);

	// number of reduced temperature points;
	min_rtemp = 3.5e-8*3.; // in cm*K;
	max_rtemp = 1.e-3*100000.;

	nb_rtemp = (int) (nb*log(max_rtemp/min_rtemp)) + 1;
	nb_chr = 1301; // nb of charge ratio values;

	chr_arr = new double [nb_chr];
	rtemp_arr = new double [nb_rtemp];

	accr_rate_func = alloc_2d_array<double>(nb_chr, nb_rtemp);
	energy_removal_func = alloc_2d_array<double>(nb_chr, nb_rtemp);

	chr_arr[nb_chr/2] = 0.;
	for (j = 1; j < nb_chr/2 + 1; j++)
	{
		if (j <= 30) {
			chr_arr[nb_chr/2 + j] = j;
		}
		else if (j <= 100) {
			chr_arr[nb_chr/2 + j] = chr_arr[nb_chr/2 + j - 1] + 2.;
		}
		else if (j <= 300) {
			chr_arr[nb_chr/2 + j] = chr_arr[nb_chr/2 + j - 1] + 5.;
		}
		else if (j <= 500) {
			chr_arr[nb_chr/2 + j] = chr_arr[nb_chr/2 + j - 1] + 10.;
		}
		else {
			chr_arr[nb_chr/2 + j] = chr_arr[nb_chr/2 + j - 1] + 50.;
		}
	}
	for (j = 1; j < nb_chr/2 + 1; j++) {
		chr_arr[nb_chr/2 - j] = -chr_arr[nb_chr/2 + j];
	}

	for (t = min_rtemp, i = 0; i < nb_rtemp; t *= sc, i++) 
	{
		rtemp_arr[i] = t;
		accr_rate_func[0][i] = energy_removal_func[0][i] = 0.;
		for (j = 1; j < nb_chr-1; j++)
		{
			accr_rate_func[j][i] = calc_accr_rate(t, chr_arr[j]); 
			energy_removal_func[j][i] = calc_energy_removal(t, chr_arr[j]);
		}
		accr_rate_func[j][i] = energy_removal_func[j][i] = 0.;
	}
}

accretion_rate_functions::~accretion_rate_functions()
{
	delete [] chr_arr;
	delete [] rtemp_arr;

	free_2d_array(accr_rate_func);
	free_2d_array(energy_removal_func);
}

void accretion_rate_functions::save(const string &path) const
{
	int i, j;	
	string fname;
	ofstream output;

	fname = path + "accret_rate_coeff.txt";
	output.open(fname.c_str());
	output << scientific;
	output.precision(2);
	
	output << "# Reduced rate coefficient (Draine, Sutin, ApJ 320, p. 803, 1987);" << endl 
		<< "# Reduced temperature K*cm, J(), L();" << endl
		<< "# nb of temperature values, nb of charge ratio values;" << endl;
	output << left << setw(5) << nb_rtemp << setw(5) << nb_chr << endl;

	output << left << setw(10) << "";
	for (j = 0; j < nb_chr; j++) {
		output << left << setw(10) << chr_arr[j] << setw(10) << chr_arr[j];
	}
	output << endl;

	for (i = 0; i < nb_rtemp; i++) 
	{
		output << left << setw(10) << rtemp_arr[i];
		for (j = 0; j < nb_chr; j++) {
			output << left << setw(10) << accr_rate_func[j][i] << setw(10) << energy_removal_func[j][i];
		}
		if (i < nb_rtemp-1)
			output << endl;
	}
	output.close();
}

double accretion_rate_functions::calc_energy_removal(double rtemp, double charge_ratio)
{
	double a, tau;
	tau = rtemp *constant;

	if (fabs(charge_ratio) < DBL_EPSILON) {
		a = 2. + 1.5*sqrt(0.5*M_PI/tau);
	}
	else if (charge_ratio > 0.) 
	{ // tau > 0.001, error 5%; 
		double theta = charge_ratio/(1. + pow(charge_ratio, -0.5));
		a = (2. + charge_ratio/tau) *(1. + pow(1.5/tau + 3.*charge_ratio, -0.5)) *exp(-theta/tau);
	}
	else { // tau > 0.1, error 10%; 
		a = (2. - charge_ratio/tau)*(1. + pow(tau - charge_ratio, -0.5));
	}
	return a;
}

double accretion_rate_functions::calc_accr_rate(double rtemp, double charge_ratio)
{
	double a, tau;
	tau = rtemp *constant;

	if (fabs(charge_ratio) < DBL_EPSILON) {
		a = 1. + sqrt(0.5*M_PI/tau);
	}
	else if (charge_ratio > 0.)
	{ 
		double theta = charge_ratio/(1. + pow(charge_ratio, -0.5)); // accuracy 0.7%
		a = 1. + pow(4.*tau + 3.*charge_ratio, -0.5);
		a = a*a *exp(-theta/tau); // is accurate within 4%
	}
	else { // is valid for tau > 1.e-3 within 5%; a = 4 A, T = 10 K, tau = 2.e-4
		a = (1. - charge_ratio/tau)*(1. + sqrt(2./(tau - 2.*charge_ratio)));
	}
	return a;
}

double accretion_rate_functions::get_accr_rate(double rt, double charge_ratio) const
{
	int k, l;
	double t, u;
	
	locate_index(chr_arr, nb_chr, charge_ratio, k);

	if (k < 0) {
		k = 0;
		t = 0.;
	}
	else if (k > nb_chr-2) {
		k = nb_chr-2;
		t = 1.;
	}
	else t = (charge_ratio - chr_arr[k])/(chr_arr[k+1] - chr_arr[k]);
	
	locate_index(rtemp_arr, nb_rtemp, rt, l);

	if (l < 0) {
		l = 0;
		u = 0.;
	}
	else if (l > nb_rtemp-2) { 
		l = nb_rtemp-2;
		u = 1.;
	}
	else u = (rt - rtemp_arr[l])/(rtemp_arr[l+1] - rtemp_arr[l]);
	
	return accr_rate_func[k][l]*(1.-t)*(1.-u) + accr_rate_func[k+1][l]*t*(1.-u) 
		+ accr_rate_func[k][l+1]*(1.-t)*u + accr_rate_func[k+1][l+1]*u*t;
}

double accretion_rate_functions::get_energy_removal(double rt, double charge_ratio) const
{
	int k, l;
	double t, u;
	
	locate_index(chr_arr, nb_chr, charge_ratio, k);

	if (k < 0) {
		k = 0;
		t = 0.;
	}
	else if (k > nb_chr-2) {
		k = nb_chr-2;
		t = 1.;
	}
	else t = (charge_ratio - chr_arr[k])/(chr_arr[k+1] - chr_arr[k]);
	
	locate_index(rtemp_arr, nb_rtemp, rt, l);

	if (l < 0) {
		l = 0;
		u = 0.;
	}
	else if (l > nb_rtemp-2) { 
		l = nb_rtemp-2;
		u = 1.;
	}
	else u = (rt - rtemp_arr[l])/(rtemp_arr[l+1] - rtemp_arr[l]);
	
	return energy_removal_func[k][l] *(1.-t)*(1.-u) + energy_removal_func[k+1][l]*t*(1.-u) 
		+ energy_removal_func[k][l+1]*(1.-t)*u + energy_removal_func[k+1][l+1]*u*t;
}

//
// Functions
//

void construct_gas_grain_reactions(string input_file, string output_path)
{
	int i, j;
	char text_line[MAX_TEXT_LINE_WIDTH];
	double photodes_yield;

	string sp_name, fname;
	stringstream ss;
	ifstream input;
	ofstream out1, out2;
	
	fname = output_path + "reactions_adsorp_desorption.txt";
	out1.open(fname.c_str(), ios_base::out);
	out1 << scientific;
	out1.precision(3);

	// Sputtering of grain mantles:
	fname = output_path + "reactions_grain_mantle_sputt.txt";
	out2.open(fname.c_str(), ios_base::out);
	out2 << scientific;
	out2.precision(3);
	
	input.open(input_file.c_str(), ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with chemical species " << input_file << endl;
		exit(1);
	}

	i = j = 1;
	out1 << "# Adsorption-desorption reactions, only photodesorption yields are given in parameter list;";
	out2 << "# Sputtering of grain mantle species, no reaction parameters are given;";

	while (!input.eof())
	{
		// comment lines are read:
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '#');

		if (text_line[0] == '\0') // check for empty line at the file end;
			break;

		ss.clear(); // the problem of the symbol of line end;
		memmove(text_line, text_line+2, strlen(text_line)-2); // for data by Penteado et al., ApJ 844, p.71 (2017);
		
		ss.str(text_line);
		// only specimen name is needed:
		ss >> sp_name;
		
		// Adsorption:
		// no parameters are saved;
		out1 << endl << i++ << ":" << "ADS:" << sp_name << ":GRAIN:*" << sp_name 
			<< ":GRAIN:::1:0.:0.:0.:3:100000:::::";

		// Thermal desorption:
		// no parameters are saved;
		out1 << endl << i++ << ":" << "THD:*" << sp_name << ":GRAIN:" << sp_name 
			<< ":GRAIN:::1:0.:0.:0.:3:100000:::::";

		// CR desorption: 
		// reaction code - CPD; no parameters are saved;
		out1 << endl << i++ << ":" << "CPD:*" << sp_name << ":GRAIN:" << sp_name 
			<< ":GRAIN:::1:0.:0.:0.:3:100000:::::";
		
		// photodesorption
		// the unique value of 1.e-4 was adopted for all species by Ruaud et al., MNRAS 459, p. 3756 (2016);
		photodes_yield = 1.e-4; // molecules per photon;

		// the discussion of photodesorption rates:
		// Hollenbach et al. ApJ 690, p. 1497, 2009; Walsh et al. ApJ 747, p. 114, 2012;
		// Oberg et al. A&A 496, p. 281, 2009a; Oberg et al. ApJ 693, p. 1209, 2009b;
		if (sp_name != "O" && sp_name != "H2O" && sp_name != "CO" && sp_name != "CO2" && sp_name != "N2")
		{
			out1 << endl << i++ << ":PHD:*" << sp_name << ":PHOTON:" << sp_name << "::::1:" << photodes_yield << ":0.:0.:3:100000:::::";
			out1 << endl << i++ << ":CRD:*" << sp_name << ":CRPHOT:" << sp_name << "::::1:" << photodes_yield << ":0.:0.:3:100000:::::";
		}
		
		// sputtering
		// maximum temperature does not matter here, no parameters of the reaction are saved;
		out2 << endl << j++ << ":" << "SPU:*" << sp_name << ":H2:" << sp_name 
			<< ":H2:::1:0.:0.:0.:3:100000:::::";

		out2 << endl << j++ << ":" << "SPU:*" << sp_name << ":He:" << sp_name 
			<< ":He:::1:0.:0.:0.:3:100000:::::";

		// Jimenez-Serra et al., A&A 482, 549-559 (2008);
		out2 << endl << j++ << ":" << "SPU:*" << sp_name << ":CO:" << sp_name 
			<< ":CO:::1:0.:0.:0.:3:100000:::::";
	}
	// Photodesorption, individually:
	// Atomic O (Hollenbach et al., 2009);
	out1 << endl << i++ << ":PHD:*O:PHOTON:O::::1:1.e-4:0.:0.:3:100000:::\"Hollenbach_ApJ690_1497_2009\"::";
	out1 << endl << i++ << ":CRD:*O:CRPHOT:O::::1:1.e-4:0.:0.:3:100000:::\"Hollenbach_ApJ690_1497_2009\"::";

	// H2O photodesorption value according to Oberg et al., 2009b;
	// photodesorption of water is temperature dependent;
	// reaction route probabilities are taken to be equal (see also discussion by Hollenbach et al., 2009);
	out1 << endl << i++ << ":PHD:*H2O:PHOTON:H2O::::1:1.e-3:0.:0.:3:100000:::\"Oberg_ApJ693_1209_2009\"::";
	out1 << endl << i++ << ":CRD:*H2O:CRPHOT:H2O::::1:1.e-3:0.:0.:3:100000:::\"Oberg_ApJ693_1209_2009\"::";
	
	out1 << endl << i++ << ":PHD:*H2O:PHOTON:H:OH:::1:1.e-3:0.:0.:3:100000:::\"Oberg_ApJ693_1209_2009\"::";
	out1 << endl << i++ << ":CRD:*H2O:CRPHOT:H:OH:::1:1.e-3:0.:0.:3:100000:::\"Oberg_ApJ693_1209_2009\"::";

	// CO photodesorption (Oberg et al., 2009a)
	out1 << endl << i++ << ":PHD:*CO:PHOTON:CO::::1:2.7e-3:0.:0.:3:100000:::\"Oberg_A&A496_281_2009\"::";
	out1 << endl << i++ << ":CRD:*CO:CRPHOT:CO::::1:2.7e-3:0.:0.:3:100000:::\"Oberg_A&A496_281_2009\"::";
	
	// CO2 photodesorption (Oberg et al., 2009a)
	out1 << endl << i++ << ":PHD:*CO2:PHOTON:CO2::::1:1.2e-3:0:0:3:100000:::\"Oberg_A&A496_281_2009\"::";
	out1 << endl << i++ << ":CRD:*CO2:CRPHOT:CO2::::1:1.2e-3:0:0:3:100000:::\"Oberg_A&A496_281_2009\"::";
	
	out1 << endl << i++ << ":PHD:*CO2:PHOTON:CO:O:::1:1.1e-3:0:0:3:100000:::\"Oberg_A&A496_281_2009\"::";
	out1 << endl << i++ << ":CRD:*CO2:CRPHOT:CO:O:::1:1.1e-3:0:0:3:100000:::\"Oberg_A&A496_281_2009\"::";

	// N2 photodesorption (Oberg et al., 2009a), the upper limit is given in the paper;
	out1 << endl << i++ << ":PHD:*N2:PHOTON:N2::::1:2.e-4:0:0:3:100000:::\"Oberg_A&A496_281_2009\"::";
	out1 << endl << i++ << ":CRD:*N2:CRPHOT:N2::::1:2.e-4:0:0:3:100000:::\"Oberg_A&A496_281_2009\"::";

	out1.close();
	out2.close();
	input.close();
}

static const int nb_of_ions = 267;
static int ion_masses[nb_of_ions] = 
   {1,         2,         3,         4,         5,         12,        13,        14,        14,        15,
	15,        16,        16,        16,        17,        17,        17,        18,        18,        19,
	19,        20,        21,        23,        24,        24,        25,        26,        26,        27,
	27,        28,        28,        28,        28,        28,        28,        29,        29,        29, 
	29,        29,        30,        30,        30,        30,        30,        30,        31,        31,
	31,        31,        31,        31,        32,        32,        32,        32,        32,        32, 
	33,        33,        33,        33,        33,        34,        34,        35,        35,        36, 
	36,        37,        37,        38,        38,        38,        39,        39,        39,        40, 
	40,        40,        40,        41,        41,        41,        41,        42,        42,        42,
	42,        42,        42,        43,        43,        43,        43,        43,        43,        43,
	43,        43,        43,        44,        44,        44,        44,        44,        44,        44, 
	44,        44,        44,        44,        44,        44,        45,        45,        45,        45,
	45,        45,        45,        46,        46,        46,        46,        46,        46,        46,
	46,        46,        47,        47,        47,        47,        47,        47,        47,        47, 
    47,        47,        47,        48,        48,        48,        48,        49,        49,        49, 
    49,        50,        50,        51,        51,        51,        52,        52,        52,        52,
	53,        53,        53,        53,        53,        54,        54,        54,        54,        55, 
	55,        55,        55,        55,        56,        56,        56,        56,        57,        57, 
    58,        58,        59,        59,        59,        60,        60,        60,        60,        61,
    61,        61,        61,        61,        62,        62,        63,        63,        64,        64,        64, 
    64,        64,        65,        65,        65,        65,        65,        66,        66,        66,
    67,        67,        68,        68,        69,        72,        73,        74,        74,        75,
    75,        76,        76,        76,        77,        77,        77,        78,        79,        79,
    80,        80,        81,        84,        85,        86,        87,        88,        89,        90,
    96,        97,        98,        98,        99,        99,        100,       100,       101,       101, 
    108,       109,       110,       111,       112,       113,       114,       120,       121,       122, 
    122,       123,       123,       124,       125,       132};
static string ion_names[nb_of_ions] = 
   {"H+",      "H2+",     "H3+",     "He+",     "HeH+",    "C+",      "CH+",     "CH2+",    "N+",      "NH+",
	"CH3+",    "NH2+",    "O+",      "CH4+",    "OH+",     "NH3+",    "CH5+",    "NH4+",    "H2O+",    "H3O+", 
	"F+",      "HF+",     "H2F+",    "Na+",     "C2+",     "Mg+",     "C2H+",    "C2H2+",   "CN+",     "C2H3+", 
	"HCN+",    "CO+",     "Si+",     "N2+",     "H2NC+",   "HCNH+",   "C2H4+",   "SiH+",    "HCO+",    "HOC+",
	"C2H5+",   "N2H+",    "H2CO+",   "SiH2+",   "CH4N+",   "NO+",     "CH2NH2+", "CH3CH3+", "SiH3+",   "H3CO+",
	"C2H7+",   "P+",      "HNO+",    "CF+",     "S+",      "O2+",     "SiH4+",   "CH3OH+",  "H2NO+",   "PH+",
	"PH2+",    "O2H+",    "SiH5+",   "CH3OH2+", "HS+",     "H2S+",    "PH3+",    "H3S+",    "Cl+",     "C3+",
	"HCl+",    "H2Cl+",   "C3H+",    "C3H2+",   "C2N+",    "CNC+",    "CH2CCH+", "C2NH+",   "C3H3+",   "C2O+",
	"CH2CN+",  "C3H4+",   "SiC+",    "C3H5+",   "CH3CN+",  "HC2O+",   "HCSi+",   "CH2CO+",  "OCN+",    "SiCH2+",
	"C3H6+",   "SiN+",    "CH3CNH+", "C3H7+",   "CP+",     "HNSi+",   "HCNO+",   "SiCH3+",  "CH3CO+",  "HONC+", 
	"HNCO+",   "HOCN+",   "NH2CNH+", "N2O+",    "SiO+",    "HCNOH+",  "CS+",     "H2CNO+",  "H2NCO+",  "CO2+", 
	"SiCH4+",  "HCP+",    "H2OCN+",  "CH3CHO+", "HNCOH+",  "SiNH2+",  "HCS+",    "SiOH+",   "HN2O+",   "PN+",
	"CH3CHOH+","HCO2+",   "PCH2+",   "HPN+",    "HCOOH+",  "C2H5OH+", "PCH3+",   "NS+",     "H2CS+",   "NO2+",
	"H2SiO+",  "CH3OCH3+","HNS+",    "H3SiO+",  "SiF+",    "PNH2+",   "PCH4+",   "CCl+",    "PO+",     "H3CS+", 
	"CH3OCH4+","C2H5OH2+","HCOOH2+", "HPO+",    "C4+",     "PNH3+",   "SO+",     "C4H+",    "H2PO+",   "H2CCl+", 
    "HSO+",    "C3N+",    "C4H2+",   "HC3N+",   "ClO+",    "C4H3+",   "C2N2+",   "C4H4+",   "SiC2+",   "C3O+",
	"C4H5+",   "HC3O+",   "SiC2H+",  "CH2CHCN+","NCCNH+",  "SiNC+",   "C3H2O+",  "CH2CHCNH+","SiC2H2+","CCP+",
	"SiC2H3+", "H3C3O+",  "C4H7+",   "SiNCH+",  "Fe+",     "C2H5CNH+","C2S+",    "HC2P+",   "HC2S+",   "PC2H2+", 
    "CH3COCH3+","PC2H3+", "CH3CS+",  "PC2H4+",  "CH3COCH4+","COOCH4+","OCS+",    "C5+",     "SiS+",    "HOCS+", 
    "HSiS+",   "C5H+",    "HSiO2+",  "H5C2O2+", "C5H2+",   "C4N+",    "C5H3+",   "HC4N+",   "H2C4N+",  "SiC3+",   "CH3C4H+", 
    "SO2+",    "S2+",     "C5H5+",   "HSO2+",   "HS2+",    "CH3C3N+", "SiC3H+",  "H2S2+",   "CH3C3NH+","SiC3H2+",
    "NCCNCH3+","H3S2+",   "PC3H+",   "C3S+",    "HC3S+",   "C6+",     "C6H+",    "C5N+",    "C6H2+",   "HC5N+",
    "C6H3+",   "HC5NH+",  "C6H4+",   "SiC4+",   "SiC4H+",  "C6H5+",   "H3C5N+",  "C6H6+",   "C6H7+",   "C4P+",
    "PC4H+",   "C4S+",    "HC4S+",   "C7+",     "C7H+",    "C7H2+",   "C7H3+",   "C7H4+",   "C7H5+",   "CH3C5NH+",
    "C8+",     "C8H+",    "C7N+",    "C8H2+",   "HC7N+",   "C8H3+",   "C8H4+",   "H2C7N+",  "H3C7N+",  "C8H5+", 
    "C9+",     "C9H+",    "C9H2+",   "C9H3+",   "C9H4+",   "C9H5+",   "CH3C7NH+","C10+",    "C10H+",   "C10H2+", 
    "C9N+",    "HC9N+",   "C10H3+",  "H2C9N+",  "H3C9N+",  "C11+"
};

// H3CO+ is the ion of CH3O
static string diss_recomb_first_channel[nb_of_ions] =
   {"H",       "H2",      "",        "He",      "",        "C",       "CH",      "CH2",     "N",       "NH",
	"CH3",     "NH2",     "O",       "CH4",     "OH",      "NH3",     "",        "",        "H2O",     "", 
	"F",       "HF",      "",        "Na",      "C2",      "Mg",      "C2H",     "C2H2",    "CN",      "C2H3", 
	"HCN",     "CO",      "Si",      "N2",      "",        "",        "C2H4",    "SiH",     "HCO",     "",
	"C2H5",    "",        "H2CO",    "SiH2",    "",        "NO",      "",        "CH3CH3",  "SiH3",    "CH3O",
	"",        "P",       "HNO",     "",        "S",       "O2",      "SiH4",    "CH3OH",   "",        "PH",
	"PH2",     "O2H",     "",        "",        "HS",      "H2S",     "",        "",        "Cl",      "C3",
	"HCl",     "",        "C3H",     "C3H2",    "C2N",     "",        "CH2CCH",  "",        "",        "C2O",
	"CH2CN",   "CH3CCH",  "SiC",     "",        "CH3CN",   "",        "HCSi",    "CH2CO",   "OCN",     "SiCH2",
	"",        "SiN",     "",        "",        "CP",      "HNSi",    "HCNO",    "SiCH3",   "",        "HONC", 
	"HNCO",    "HOCN",    "",        "N2O",     "SiO",     "",        "CS",      "",        "",        "CO2", 
	"",        "HCP",     "",        "CH3CHO",  "",        "",        "HCS",     "",        "",        "PN",
	"",        "",        "",        "",        "HCOOH",   "C2H5OH",  "",        "NS",      "H2CS",    "NO2",
	"H2SiO",   "CH3OCH3", "",        "",        "",        "",        "",        "CCl",     "PO",      "", 
	"",        "",        "",        "HPO",     "C4",      "",        "SO",      "C4H",     "",        "", 
    "",        "C3N",     "HC4H",    "HC3N",    "ClO",     "C4H3",    "NCCN",    "",        "SiC2",    "C3O", 
    "",        "",        "SiC2H",   "CH2CHCN", "",        "SiNC",    "",        "",        "SiC2H2",  "CCP",
    "",        "",        "",        "",        "Fe",      "",        "C2S",     "HC2P",    "",        "",
    "CH3COCH3","",        "",        "",        "",        "HCOOCH3", "OCS",     "C5",      "SiS",     "",
    "",        "C5H",     "",        "",        "C5H2",    "C4N",     "",        "",        "",        "SiC3",   "CH3C4H",
    "SO2",     "S2",      "",        "",        "HS2",     "CH3C3N",  "SiC3H",   "H2S2",    "",        "",
    "",        "",        "",        "C3S",     "",        "C6",      "C6H",     "C5N",     "C6H2",    "HC5N", 
    "",        "",        "",        "SiC4",    "",        "",        "",        "C6H6",    "",        "C4P", 
    "",        "C4S",     "",        "C7",      "C7H",     "C7H2",    "",        "CH3C6H",  "",        "", 
    "C8",      "C8H",     "C7N",     "C8H2",    "HC7N",    "",        "",        "",        "",        "", 
    "C9",      "C9H",     "C9H2",    "",        "",        "",        "",        "C10",     "C10H",    "C10H2", 
    "C9N",     "HC9N",    "",        "",        "",        "C11"
};

static string diss_recomb_second_channel_1[nb_of_ions] = 
   {"",        "H",       "H2",      "",       "He",       "",        "C",       "CH",      "",        "N",
	"CH2",     "NH",      "",        "CH3",     "O",       "NH2",     "CH4",     "NH3",     "OH",      "H2O", 
	"",        "F",       "HF",      "",        "",        "",        "C2",      "C2H",     "",        "C2H2", 
	"CN",      "",        "",        "",        "HNC",     "HCN",     "C2H3",    "Si",      "CO",      "CO",
	"C2H4",    "N2",      "HCO",     "SiH",     "CH2NH",   "",        "CH2NH",   "C2H5",    "SiH2",    "H2CO",
	"CH3CH3",  "",        "NO",      "",        "",        "",        "SiH3",    "CH3O",    "HNO",     "P",
	"PH",      "O2",      "SiH4",    "CH3OH",   "S",       "HS",      "PH2",     "H2S",     "",        "",
	"Cl",      "HCl",     "C3",      "C3H",     "",        "",        "H2CCC",   "C2N",     "C3H2",    "",
	"",        "",        "",        "CH3CCH",  "CH2CN",   "C2O",     "SiC",     "",        "",        "HCSi",
	"",        "",        "CH3CN",   "CH3CHCH2","",        "SiN",     "CNO",     "SiCH2",   "CH2CO",   "CNO", 
	"OCN",     "OCN",     "NH2CN",   "",        "",        "HONC",    "",        "HCNO",    "HNCO",    "", 
	"SiCH3",   "CP",      "HOCN",    "",        "HOCN",    "HNSi",    "CS",      "SiO",     "N2O",     "",
	"CH3CHO",  "CO2",     "HCP",     "PN",      "",        "",        "",        "",        "HCS",     "",
	"",        "",        "NS",      "H2SiO",   "",        "",        "CH2PH",   "",        "",        "H2CS", 
	"CH3OCH3", "C2H5OH",  "HCOOH",   "PO",      "",        "",        "",        "C4",      "HPO",     "", 
    "SO",      "",        "C4H",     "C3N",     "",        "HC4H",    "",        "C4H3",    "",        "",
    "",        "C3O",     "SiC2",    "",        "NCCN",    "",        "",        "CH2CHCN", "SiC2H",   "",
    "SiC2H2",  "",        "",        "SiNC",    "",        "C2H5CN",  "",        "CCP",     "C2S",     "HC2P",
    "",        "",        "",        "",        "CH3COCH3","",        "",        "",        "",        "OCS",
    "SiS",     "C5",      "SiO2",    "HCOOCH3", "C5H",     "",        "C5H2",    "C4N",     "",        "",      "",
    "",        "",        "CH3C4H",  "SO2",     "S2",      "",        "SiC3",    "HS2",     "CH3C3N",  "SiC3H",
    "",        "H2S2",    "C3P",     "",        "C3S",     "",        "C6",      "",        "C6H",     "C5N", 
    "C6H2",    "HC5N",    "",        "",        "SiC4",    "",        "",        "",        "C6H6",    "", 
    "C4P",     "",        "C4S",     "",        "C7",      "C7H",     "C7H2",    "",        "CH3C6H",  "CH3C5N", 
    "",        "C8",      "",        "C8H",     "C7N",     "C8H2",    "",        "HC7N",    "",        "", 
    "",        "C9",      "C9H",     "C9H2",    "",        "",        "CH3C7N",  "",        "C10",     "C10H", 
    "",        "C9N",     "C10H2",   "HC9N",    "",    ""
};

void construct_ion_recomb_grains(string output_path)
{
	const double rate = sqrt( 8.*BOLTZMANN_CONSTANT/(M_PI*ATOMIC_MASS_UNIT) );
	int i, j, l;
	
	string fname;
	stringstream ss;
	ofstream out;

	fname = output_path + "reactions_ion_recomb_grains.txt";
	out.open(fname.c_str(), ios_base::out);
	out << scientific;
	out.precision(3);
	
	out << "# Ion recombination on grains; Wakelam & Herbst, ApJ 680, p.371, 2008" << endl
		<< "# H3CO+ is the ion of CH3O; parameter*sqrt(T[K]) = sticking coeff * branching ratio *thermal velocity[cm/s];";
	
	i = 1;
	for (l = 0; l < nb_of_ions; l++) 
	{
		j = 0;
		if (diss_recomb_first_channel[l] != "")
			j++;
		if (diss_recomb_second_channel_1[l] != "")
			j++;

		// recombination of small molecular ions is probably dissociative (Aikawa et al., ApJ 527, p. 262, 1999); 
		if ( diss_recomb_first_channel[l] != "")
		{
			// one reaction parameter: sticking coefficient * branching ratio * rate, 
			// rate*sqrt(T in K) is the average thermal speed in cm/s:
			out << endl << i++ << ":IAG:" << ion_names[l] << ":GRAIN:" << diss_recomb_first_channel[l] << ":GRAIN:::1:" 
				<< rate/(sqrt((double) ion_masses[l])*j) << ":0:0:3:1000000:::::";
		}
		if ( diss_recomb_second_channel_1[l] != "")
		{
			out << endl << i++ << ":IAG:" << ion_names[l] << ":GRAIN:" << diss_recomb_second_channel_1[l] << ":H:GRAIN::1:"
				<< rate/(sqrt((double) ion_masses[l])*j) << ":0:0:3:1000000:::::";
		}			
	}

	ss.str("");
	for (l = 0; l < nb_of_ions; l++) 
	{
		if (diss_recomb_first_channel[l] == "" && diss_recomb_second_channel_1[l] == "")
			ss << " " << ion_names[l];
	}
	out << endl << "# Ions that have special recombination reactions:" << endl << "#" << ss.str();
	out << endl << i++ << ":IAG:CF+:GRAIN:C:F:GRAIN::1:" << rate/sqrt(31.) << ":0:0:3:1000000:::::";
	out << endl << i++ << ":IAG:CNC+:GRAIN:CN:C:GRAIN::1:" << rate/sqrt(38.) << ":0:0:3:1000000:::::";
	out << endl << i++ << ":IAG:C3H6+:GRAIN:CH3CCH:H2:GRAIN::1:" << rate/sqrt(42.) << ":0:0:3:1000000:::::";
	out << endl << i++ << ":IAG:PCH3+:GRAIN:P:CH3:GRAIN::1:" << rate/sqrt(46.) << ":0:0:3:1000000:::::";
	out << endl << i++ << ":IAG:SiF+:GRAIN:Si:F:GRAIN::1:" << rate/sqrt(47.) << ":0:0:3:1000000:::::";
	out << endl << i++ << ":IAG:PNH2+:GRAIN:P:NH2:GRAIN::1:" << rate/sqrt(47.) << ":0:0:3:1000000:::::";
	out << endl << i++ << ":IAG:PNH3+:GRAIN:P:NH3:GRAIN::1:" << rate/sqrt(48.) << ":0:0:3:1000000:::::";
	out << endl << i++ << ":IAG:H2CCl+:GRAIN:CCl:H:H:GRAIN:1:" << rate/sqrt(49.) << ":0:0:3:1000000:::::";
	out << endl << i++ << ":IAG:C4H5+:GRAIN:CH3CCH:CH:GRAIN::1:" << rate/sqrt(53.) << ":0:0:3:1000000:::::";
	out << endl << i++ << ":IAG:C3H2O+:GRAIN:CO:C2H2:GRAIN::1:" << rate/sqrt(54.) << ":0:0:3:1000000:::::";
	out << endl << i++ << ":IAG:H3C3O+:GRAIN:C3O:H2:H:GRAIN:1:" << rate/sqrt(55.) << ":0:0:3:1000000:::::";
	out << endl << i++ << ":IAG:C4H7+:GRAIN:C4H3:H2:H2:GRAIN:1:" << rate/sqrt(55.) << ":0:0:3:1000000:::::";
	out << endl << i++ << ":IAG:PC2H3+:GRAIN:HC2P:H2:GRAIN::1:" << rate/sqrt(58.) << ":0:0:3:1000000:::::";
	out << endl << i++ << ":IAG:CH3CS+:GRAIN:CH3:CS:GRAIN::1:" << rate/sqrt(59.) << ":0:0:3:1000000:::::";
	out << endl << i++ << ":IAG:PC2H4+:GRAIN:P:C2H4:GRAIN::1:" << rate/sqrt(59.) << ":0:0:3:1000000:::::";
	out << endl << i++ << ":IAG:H2C4N+:GRAIN:HC3N:CH:GRAIN::1:" << rate/sqrt(64.) << ":0:0:3:1000000:::::";
	out << endl << i++ << ":IAG:NCCNCH3+:GRAIN:NCCN:CH3:GRAIN::1:" << rate/sqrt(67.) << ":0:0:3:1000000:::::";
	out << endl << i++ << ":IAG:C6H4+:GRAIN:C6H2:H2:GRAIN::1:" << rate/sqrt(76.) << ":0:0:3:1000000:::::";
	out << endl << i++ << ":IAG:C6H5+:GRAIN:C6H:H2:H2:GRAIN:1:" << rate/sqrt(77.) << ":0:0:3:1000000:::::";
	out << endl << i++ << ":IAG:H3C5N+:GRAIN:HC5N:H2:GRAIN::1:" << rate/sqrt(77.) << ":0:0:3:1000000:::::";
	out << endl << i++ << ":IAG:C8H4+:GRAIN:C8H2:H2:GRAIN::1:" << 0.1*rate << ":0:0:3:1000000:::::";
	out << endl << i++ << ":IAG:H3C7N+:GRAIN:HC7N:H2:GRAIN::1:" << rate/sqrt(101.) << ":0:0:3:1000000:::::";
	out << endl << i++ << ":IAG:C8H5+:GRAIN:C8H:H2:H2:GRAIN:1:" << rate/sqrt(101.) << ":0:0:3:1000000:::::";
	out << endl << i++ << ":IAG:C9H4+:GRAIN:C9H2:H2:GRAIN::1:" << rate/sqrt(112.) << ":0:0:3:1000000:::::";
	out << endl << i++ << ":IAG:C9H5+:GRAIN:C9H:H2:H2:GRAIN:1:" << rate/sqrt(113.) << ":0:0:3:1000000:::::";
	out << endl << i++ << ":IAG:H3C9N+:GRAIN:HC9N:H2:GRAIN::1:" << rate/sqrt(125.) << ":0:0:3:1000000:::::";

	out.close();
}

// Some reactions of H and H2 destruction by ions are present in the UMIST database for astrochemistry;
// Is not used;
void construct_hydrogen_ioniz_reactions(string output_path)
{
	int i, l;
	double rate, mass, cs;

	string fname;
	ofstream out;

	fname = output_path + "reactions_hydrogen_destruct_ions.txt";
	out.open(fname.c_str(), ios_base::out);
	out << scientific;
	out.precision(3);

	i = 1;
	out << "# H2 and H ionization by ions (Draine et al. 1983); H2 dissociation by ions (Wilgenbus et al. 2000);" << endl
		<< "# representative ion, the destruction of the ion is not taken into account.";
	// only one ion is taken into account (HCO+);
	
	l = 38;
	// the H2 molecule and H atom ionization according to Draine et al., ApJ 264, 485 (1983);
	out << endl << i++ << ":NII:H2:" << ion_names[l] << ":H2+:" << ion_names[l] << ":e-::1:" << 15.4259*EV_TO_ERGS << ":1.:0" 
		<< ":3.:1.e+6:::\"Draine_Roberge_Dalgarno_ApJ264_485_1983\"::";

	out << endl << i++ << ":NII:H:" << ion_names[l] << ":H+:" << ion_names[l] << ":e-::1:" << 13.5984*EV_TO_ERGS << ":0.665:0" 
		<< ":3.:1.e+6:::\"Draine_Roberge_Dalgarno_ApJ264_485_1983\"::";
		
	// the H2 molecule dissociation rate is calculated according to Wilgenbus et al., A&A 356, p. 1010 (2000);
	cs = M_PI*7.4e-9*7.4e-9; // in cm2, the H2 internuclear distance is 7.4e-9 cm;
	mass = ATOMIC_MASS_UNIT*2.*ion_masses[l]/(2. + ion_masses[l]);
		
	// k = a*(T/300)^b*exp(-g/T)
	rate = cs *sqrt(8.*BOLTZMANN_CONSTANT*300./(M_PI*mass));

	// dissociation energy of H2 molecule is taken to equal 4.48 eV = 52000 K;
	out << endl << i++ << ":NDI:H2:" << ion_names[l] << ":H:H:" << ion_names[l] << "::1:" << rate 
		<< ":0.5:52000.:3.:1.e+6:::\"Wilgenbus_AA356_1010_2000\"::"; 
	
	out.close();
}
