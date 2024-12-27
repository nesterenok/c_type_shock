//
// 10.03.2017. Check for errors.
// 04.09.2017. Check for errors.
// 15.04.2022. Number of comment lines in the files with momentum transfer cross sections was changed,
//    new data points were added in the H2-e file.

#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif

#include <stdio.h>
#include <stdlib.h>

#include<cmath>
#include<cstring>
#include<string>
#include<fstream>

#include "constants.h"
#include "utils.h"
#include "interpolation.h"
#include "integration.h"
#include "elastic_scattering.h"

#define MAX_TEXT_LINE_WIDTH 240
#define SOURCE_NAME "elastic_scattering.cpp"
using namespace std;

// Momentum transfer cross section for elastic scattering between neutral and charged particles (Osterbrock, ApJ 134, 270, 1961);
// <cross section*velocity>_mt = osterbrock_const *Z *sqrt(polarizability [bohr radii ^3] /reduced mass [g]);
const double osterbrock_const = 2.41*M_PI*ELECTRON_CHARGE*sqrt(BOHR_RADIUS*BOHR_RADIUS*BOHR_RADIUS);

elastic_cross_section_const::elastic_cross_section_const(double polarizability, double neutral_mass, double ion_mass)
{
	double reduced_mass = neutral_mass*ion_mass/(neutral_mass + ion_mass);
	cs_osterbrock = osterbrock_const*sqrt(polarizability/reduced_mass);
}

elastic_cross_section_powerlaw::elastic_cross_section_powerlaw(double polarizability, double neutral_mass, double ion_mass, double coll_energy_lim, double gg) 
	: elastic_cross_section_const(polarizability, neutral_mass, ion_mass), g(gg)
{
	// limiting collision energy in K, 1.871 a.m.u is the reduced mass of the system H2-HCO+;
	vel0 = sqrt(2.*coll_energy_lim*BOLTZMANN_CONSTANT/(1.871*ATOMIC_MASS_UNIT));
}

double elastic_cross_section_powerlaw::get(double vel) const
{
	if (vel < vel0)
		return cs_osterbrock;
	else return cs_osterbrock *pow(vel/vel0, 1.+g); // cs*velocity [cm3/s]
}

elastic_cross_section_table::elastic_cross_section_table(const std::string &data_path, const std::string &name)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i;
	double energy, cs, reduced_mass, g;
	
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
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	// reduced mass is necessary for conversion from energy units to velocity:
	input >> g >> reduced_mass >> nb_cs;

	vel_arr = new double [nb_cs];
	cs_arr = new double [nb_cs];

	for (i = 0; i < nb_cs; i++) 
	{
		input >> energy >> cs;
		// the energy data in the file have dimension eV, conversion to the velocity in the frame of mass centre (cm/s):
		vel_arr[i] = sqrt(2.*energy*EV_TO_ERGS/reduced_mass);
		// the cross section data in the file have dimension 1.e-16 cm2:
		cs_arr[i] = cs*1.e-16;
	}
	input.close();
}

elastic_cross_section_table::~elastic_cross_section_table()
{
	delete [] vel_arr;
	delete [] cs_arr;
}

double elastic_cross_section_table::get(double vel) const
{
	int l;
	double answ;
	// if values in question are out of range of the array with a dimension dim, returned index < 0 or index >= dim-1
	locate_index(vel_arr, nb_cs, vel, l);

	if (l < 0)
		answ = cs_arr[0];
	else if (l >= nb_cs-1) 
		answ = cs_arr[nb_cs-1];
	else answ = cs_arr[l] + (vel - vel_arr[l]) *(cs_arr[l+1] - cs_arr[l])/(vel_arr[l+1] - vel_arr[l]);

	return answ*vel;
}

elastic_scattering::elastic_scattering(double m1, double m2, int v) : mass1(m1), mass2(m2), verbosity(v)
{
	mass_sum = m1 + m2;
	reduced_mass = m1*m2 /(m1 + m2);
}

void elastic_scattering::calc_data(elastic_cross_section *cross_section)
{
	if (verbosity)
		cout << "Calculating the elastic scattering data..." << endl;
	csv = cross_section->get(0.);
}

// vel_12 = vel_1 - vel_2
void elastic_scattering::calc_source_terms(double & mom_gain1, double & mom_gain2, double & energy_gain1, 
	double & energy_gain2, double conc1, double conc2, double vel_12, double temp1, double temp2)
{
	double a = conc1 *conc2 *reduced_mass *csv;
	mom_gain1 -= a *vel_12;
	mom_gain2 += a *vel_12;

	energy_gain1 += a/mass_sum *(3.*BOLTZMANN_CONSTANT*(temp2 - temp1) + mass2 *vel_12*vel_12);
	energy_gain2 += a/mass_sum *(3.*BOLTZMANN_CONSTANT*(temp1 - temp2) + mass1 *vel_12*vel_12);
}

// Static media, there is no drift between two fluids:
void elastic_scattering::calc_source_terms(double & energy_gain1, double & energy_gain2, double conc1, double conc2, 
	double temp1, double temp2)
{
	double a = conc1 *conc2 *reduced_mass *csv *3.*BOLTZMANN_CONSTANT/mass_sum;
	energy_gain1 += a*(temp2 - temp1);
	energy_gain2 += a*(temp1 - temp2);
}

// Neutral mass is the first mass (mass1), charged mass is the second (mass2);
elastic_scatt_neutral_charged::elastic_scatt_neutral_charged(double neutral_mass, double charged_mass, double velth_min, 
	double velth_max, int verbosity) : elastic_scattering(neutral_mass, charged_mass, verbosity)
{
	int i, nb = 30;
	double a = pow(10., 1./nb), s_min = 0.0001, s_max = 100.;
	
	nb_velth = (int) (nb*log10(velth_max/velth_min)) + 1;
	nb_s = (int) (nb*log10(s_max/s_min)) + 1; 

	velth_arr = new double [nb_velth];
	s_arr = new double [nb_s];
	
	c4_arr = new double [nb_velth];
	memset(c4_arr, 0, nb_velth*sizeof(double));

	c1_arr = alloc_2d_array<double>(nb_velth, nb_s);
	memset(*c1_arr, 0, nb_velth*nb_s*sizeof(double));

	c2_arr = alloc_2d_array<double>(nb_velth, nb_s);
	memset(*c2_arr, 0, nb_velth*nb_s*sizeof(double));

	c3_arr = alloc_2d_array<double>(nb_velth, nb_s);
	memset(*c3_arr, 0, nb_velth*nb_s*sizeof(double));

	// velocity in cm/s: 
	velth_arr[0] = velth_min;
	for (i = 1; i < nb_velth; i++) {
		velth_arr[i] = velth_arr[i-1]*a;
	}
	s_arr[0] = s_min;
	for (i = 1; i < nb_s; i++) {
		s_arr[i] = s_arr[i-1]*a;
	}
}

elastic_scatt_neutral_charged::~elastic_scatt_neutral_charged()
{
	delete [] velth_arr;
	delete [] s_arr;
	delete [] c4_arr;

	free_2d_array(c1_arr);
	free_2d_array(c2_arr);
	free_2d_array(c3_arr);
}

void elastic_scatt_neutral_charged::calc_data(elastic_cross_section *cross_section)
{
	int i, j;
	double a, b, s, in1, in2, in3, in4;
	
	el_func1 f1;
	el_func2 f2;
	el_func3 f3;

	el_func4 f4;
	el_func5 f5;
	el_func6 f6;
	el_func7 f7;
	
	csv = cross_section->get(0.);
	f1.cs = f2.cs = f3.cs = f4.cs = f5.cs = f6.cs = f7.cs = cross_section;

	if (verbosity)
		cout << "Calculating the elastic scattering data..." << endl;

	for (i = 0; i < nb_velth; i++) 
	{
		f1.p1 = f2.p1 = f3.p1 = velth_arr[i];
		f4.p = f5.p = f6.p = f7.p = velth_arr[i];
		
		in1 = qromb<el_func4>(f4, 0, 5., 1.e-9);
		in2 = qromb<el_func5>(f5, 0, 5., 1.e-9);
		in3 = qromb<el_func6>(f6, 0, 5., 1.e-9);
		in4 = qromb<el_func7>(f7, 0, 5., 1.e-9);

		// the data for the calculation of zero-order approximation:
		c4_arr[i] = in2*ONEDIVBY_SQRT_PI;

		// the functors return the integral value * thermal speed
		// but the thermal velocity factor is necessary for the calculation of source terms;
		for (j = 0; j < nb_s; j++) 
		{	
			s = s_arr[j];
			if (s < 0.1)
			{		
				c1_arr[i][j] = 2.*s*(in1 + 2./3.*s*s *in2 + 2./15.*s*s*s*s *in3 + 4./315*s*s*s*s*s*s *in4);
				c2_arr[i][j] = in1 + 2.*s*s *in2 + 2./3.*s*s*s*s *in3 + 4./45.*s*s*s*s*s*s *in4;
				c3_arr[i][j] = 2.*s*(in2 + 2./3.*s*s *in3 + 2./15.*s*s*s*s *in4);
				
				c1_arr[i][j] *= exp(-s*s);
				c2_arr[i][j] *= exp(-s*s);
				c3_arr[i][j] *= exp(-s*s);
			}
			else 
			{
				b = s + 5.;
				if (s > 5.) 
					a = s - 5.;
				else a = 0.;

				f1.p2 = f2.p2 = f3.p2 = s;

				// The integral is split on two parts:
				c1_arr[i][j] = 0.5*qromb<el_func1>(f1, a, b, 1.e-7);
				c2_arr[i][j] = 0.5*qromb<el_func2>(f2, a, b, 1.e-7);
				c3_arr[i][j] = 0.5*qromb<el_func3>(f3, a, b, 1.e-7);
				
				if (s < 5.) // 
				{
					f1.p2 = f2.p2 = f3.p2 = -s;
		
					c1_arr[i][j] -= 0.5*qromb<el_func1>(f1, 0., 5., 1.e-7);
					c2_arr[i][j] += 0.5*qromb<el_func2>(f2, 0., 5., 1.e-7); // Note the sign '+' - the cosh is evaluated;
					c3_arr[i][j] -= 0.5*qromb<el_func3>(f3, 0., 5., 1.e-7);
				}
			}
			
			c1_arr[i][j] *= ONEDIVBY_SQRT_PI;
			c2_arr[i][j] *= ONEDIVBY_SQRT_PI;
			c3_arr[i][j] *= ONEDIVBY_SQRT_PI;
		}
	}
}

void elastic_scatt_neutral_charged::save_data(const string & fname)
{
	int i, j;
	ofstream output;
	
	output.open(fname.c_str(), ios_base::out);
	output << scientific;
	output.precision(4);

	output << "x axis is for s parameter, y axis - reduced thermal velocity;" << endl;
	output << "Integral: v_th*exp(-s*s)/sqrt(pi)* int_0^inf x^2 *exp(-x^2) *sinh(2xs) *sigma(x*v_th)" << endl;
	output << left << setw(14) << "";
	for (j = 0; j < nb_s; j++) {
		output << left << setw(14) << s_arr[j];
	}
	output << endl;

	for (i = 0; i < nb_velth; i++) 
	{
		output << left << setw(14) << velth_arr[i];
		for (j = 0; j < nb_s; j++) {
			output << left << setw(14) << c1_arr[i][j];
		}
		output << endl;
	}
	output << endl;
	
	output << "Integral: v_th*exp(-s*s)/sqrt(pi)* int_0^inf x^3 *exp(-x^2) *cosh(2xs) *sigma(x*v_th)" << endl;
	for (i = 0; i < nb_velth; i++) 
	{
		output << left << setw(14) << velth_arr[i];
		for (j = 0; j < nb_s; j++) {
			output << left << setw(14) << c2_arr[i][j];
		}
		output << endl;
	}
	output << endl;
	
		output << "Integral: v_th*exp(-s*s)/sqrt(pi)* int_0^inf x^4 *exp(-x^2) *sinh(2xs) *sigma(x*v_th)" << endl;
	for (i = 0; i < nb_velth; i++) 
	{
		output << left << setw(14) << velth_arr[i];
		for (j = 0; j < nb_s; j++) {
			output << left << setw(14) << c3_arr[i][j];
		}
		output << endl;
	}
	output.close();
}

// Linear interpolation is used;
void elastic_scatt_neutral_charged::get_integrals(double &c1, double &c2, double &c3, double vel, double s) const
{
	int k, l;
	double t, u;
	
	// if values in question are out of range of the array with a dimension dim, returned index < 0 or index >= dim-1
	locate_index(velth_arr, nb_velth, vel, k);
	locate_index(s_arr, nb_s, s, l);

	if (k < 0) {
		k = 0;
		t = 0.;
	}
	else if (k >= nb_velth-1) {
		k = nb_velth-2;
		t = 1.;
	}
	else t = (vel - velth_arr[k])/(velth_arr[k+1] - velth_arr[k]);
	
	if (l < 0) {
		l = 0;
		u = 0.;
	}
	else if (l >= nb_s-1) {
		l = nb_s-2;
		u = 1.;
	}
	else u = (s - s_arr[l])/(s_arr[l+1] - s_arr[l]);
	
	c1 = c1_arr[k][l]*(1.-t)*(1.-u) + c1_arr[k+1][l]*t*(1.-u) + c1_arr[k][l+1]*(1.-t)*u + c1_arr[k+1][l+1]*u*t;
	c2 = c2_arr[k][l]*(1.-t)*(1.-u) + c2_arr[k+1][l]*t*(1.-u) + c2_arr[k][l+1]*(1.-t)*u + c2_arr[k+1][l+1]*u*t;
	c3 = c3_arr[k][l]*(1.-t)*(1.-u) + c3_arr[k+1][l]*t*(1.-u) + c3_arr[k][l+1]*(1.-t)*u + c3_arr[k+1][l+1]*u*t;
}

void elastic_scatt_neutral_charged::get_integrals(double &c4, double vel) const
{
	int k;
	locate_index(velth_arr, nb_velth, vel, k);

	if (k < 0)
		c4 = c4_arr[0];
	else if (k >= nb_velth-1)
		c4 = c4_arr[nb_velth-1];
	else c4 = c4_arr[k] + (vel - velth_arr[k])*(c4_arr[k+1] - c4_arr[k])/(velth_arr[k+1] - velth_arr[k]);
}

// vel_12 = vel_n - vel_ch
void elastic_scatt_neutral_charged::calc_source_terms(double & mom_gain_n, double & mom_gain_ch, double & energy_gain_n, 
	double & energy_gain_ch, double conc_n, double conc_ch, double vel_12, double temp_n, double temp_ch)
{
	double a, b, d, vel_th, s, c1, c2, c3;
	// reduced thermal speed:
	vel_th = sqrt(2.*BOLTZMANN_CONSTANT*(mass1 *temp_ch + mass2 *temp_n)/(mass1*mass2));
	s = fabs(vel_12/vel_th);
	
	// The factors exp(-s*s)*vel_th/sqrt(pi) are taken into account in the calculation of c-arrays;
	if (s > 1.e-4)
	{
		get_integrals(c1, c2, c3, vel_th, s);
	
		a = conc_n*conc_ch *reduced_mass/s;
		b = s*c2 - 0.5*c1;
		d = 2.*a*b*vel_12/(s*s); // momentum loss of neutral component;
	
		mom_gain_n -= d;
		mom_gain_ch += d;

		a *= 4.*BOLTZMANN_CONSTANT/mass_sum;

		energy_gain_n += a*( (temp_ch - temp_n)*c3 + mass2*temp_n/reduced_mass *b );
		energy_gain_ch += a*( (temp_n - temp_ch)*c3 + mass1*temp_ch/reduced_mass *b );
	}
	else 
	{ // the zero-order approximation, relative motion between charged particles and neutrals is small;
		get_integrals(c1, vel_th);

		a = 8.*conc_n*conc_ch *reduced_mass*c1;
		d = a*vel_12/3.;

		mom_gain_n -= d;
		mom_gain_ch += d;

		a *= BOLTZMANN_CONSTANT*(temp_ch - temp_n)/mass_sum;
		energy_gain_n += a;
		energy_gain_ch -= a;
	}
}

void elastic_scatt_neutral_charged::calc_source_terms(double & energy_gain_n, double & energy_gain_ch, double conc_n, double conc_ch, double temp_n, 
		double temp_ch)
{
	double a, vel_th, c1;
	// reduced thermal velocity:
	vel_th = sqrt(2.*BOLTZMANN_CONSTANT*(mass1 *temp_ch + mass2 *temp_n)/(mass1*mass2));
	
	get_integrals(c1, vel_th);
	a = 8.*BOLTZMANN_CONSTANT*(temp_ch - temp_n)*conc_n*conc_ch *reduced_mass*c1 /mass_sum;
	
	energy_gain_n += a;
	energy_gain_ch -= a;
}
