
//
// 06.03.2017. Check for errors.
// 12.09.2017. Minor changes. Check for errors.

#include <stdio.h>
#include <stdlib.h>

#include <fstream>
#include <iostream>
#include <iomanip>
#include <memory.h>
#include <cmath>
#include <cfloat>

#include "coll_rates_co.h"
#include "utils.h"
#include "constants.h"

#define MAX_TEXT_LINE_WIDTH 240
#define SOURCE_NAME "coll_rates_co.cpp"
using namespace std;

co_h2_coll_data::co_h2_coll_data(const string path, const energy_diagram *di, bool is_ortho_h2, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, n1, n2, l, v1, v2, j1, j2, nb_lines;
	double a;
	
	string fname;
	ifstream input;

	if (is_ortho_h2) 
		fname = path + "coll_co/coll_co_oh2.txt";
	else fname = path + "coll_co/coll_co_ph2.txt";	

	input.open(fname.c_str(), ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << fname << endl;
		exit(1);
	}
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	
	input >> nb_lines >> jmax;
	jmax++; // one point is reserved for 0 K;

	// the nb of rotational levels is 41 (J <= 40) for these data, the nb in the list including vibrationally excited state - 64;
	nb_lev = di->nb_lev;
	// nb of rows in the array:
	imax = nb_lev*(nb_lev-1)/2;

	tgrid = new double [jmax];
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, jmax *imax *sizeof(double));
	
	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, imax*jmax*sizeof(double));
	
	tgrid[0] = 0.;
	for (j = 1; j < jmax; j++) {
		input >> tgrid[j];
	}

	for (i = 0; i < nb_lines; i++)
	{
		// upper - lower level;
		input >> v1 >> j1 >> v2 >> j2;
		n1 = di->get_nb(v1, j1);
		n2 = di->get_nb(v2, j2);

		if (n1 != -1 && n2 != -1)
		{
			l = n1*(n1-1)/2 + n2;
			for (j = 1; j < jmax; j++) {
				input >> coeff[l][j];
			}
		}
		else {
			for (j = 1; j < jmax; j++) {
				input >> a;
			}
		}
	}
	input.close();
	calc_coeff_deriv();
	
	if (verbosity) {
		cout << "  data have been read from file " << fname << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

co_he_coll_data::co_he_coll_data(const string path, const energy_diagram *di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, v1, v2, j1, j2, nb_lines;
	string fname;
	ifstream input;

	fname = path + "coll_co/coll_co_he.txt";
	input.open(fname.c_str(), ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << fname << endl;
		exit(1);
	}
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	
	input >> nb_lines >> jmax;
	jmax++; // one point is reserved for 0 K;

	nb_lev = 15;
	imax = nb_lev*(nb_lev-1)/2;
	
	tgrid = new double [jmax];
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, jmax *imax *sizeof(double));
	
	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, imax*jmax*sizeof(double));
	
	tgrid[0] = 0.;
	for (j = 1; j < jmax; j++) {
		input >> tgrid[j];
	}

	for (i = 0; i < nb_lines && i < imax; i++)
	{
		// upper - lower level;
		input >> v1 >> j1 >> v2 >> j2;
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

co_h_coll_data::co_h_coll_data(const string path, const energy_diagram *di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, v1, v2, j1, j2, l, l1, l2, nb_lines;
	double a;
	
	string fname;
	ifstream input;
	
	// maximal level in the data is (v,j) = (0,45), it's number is 77 in the total level list;
	nb_lev = di->nb_lev;
	imax = nb_lev*(nb_lev-1)/2;

	fname = path + "coll_co/coll_co_h_vibr0-0.txt";
	input.open(fname.c_str(), ios_base::in);
	
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << fname << endl;
		exit(1);
	}
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	
	input >> nb_lines >> jmax;
	jmax++; // one point is reserved for 0 K;

	tgrid = new double [jmax];

	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, jmax *imax *sizeof(double));
	
	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, imax*jmax*sizeof(double));
			
	tgrid[0] = 0.;
	for (j = 1; j < jmax; j++) {
		input >> tgrid[j];
	}
		
	for (i = 0; i < nb_lines; i++)
	{
		// upper - lower level;
		input >> v1 >> j1 >> v2 >> j2;
			
		l1 = di->get_nb(v1, j1);
		l2 = di->get_nb(v2, j2);

		if (l1 != -1 && l2 != -1)
		{
			l = l1*(l1-1)/2 + l2;
			for (j = 1; j < jmax; j++) {
				input >> coeff[l][j];
			}
		}
		else {
			for (j = 1; j < jmax; j++) {
				input >> a;
			}
		}		
	}
	input.close();	
	
	calc_coeff_deriv();
	if (verbosity) {
		cout << "  data have been read from file " << fname << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

co_h_vibr_coll_data::co_h_vibr_coll_data(const string path, const energy_diagram *di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, v1, v2, j1, j2, l, l1, l2, nb_lines;
	double a;
	
	string fname;
	ifstream input;
	
	nb_lev = di->nb_lev;
	imax = nb_lev*(nb_lev-1)/2;

	fname = path + "coll_co/coll_co_h_vibr.txt";
	input.open(fname.c_str(), ios_base::in);
	
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << fname << endl;
		exit(1);
	}
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	
	input >> nb_lines >> jmax;
	jmax++; // one point is reserved for 0 K;

	tgrid = new double [jmax];

	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, jmax *imax *sizeof(double));
	
	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, imax*jmax*sizeof(double));
			
	tgrid[0] = 0.;
	for (j = 1; j < jmax; j++) {
		input >> tgrid[j];
	}
		
	for (i = 0; i < nb_lines; i++)
	{
		// upper - lower level;
		input >> v1 >> j1 >> v2 >> j2;
			
		l1 = di->get_nb(v1, j1);
		l2 = di->get_nb(v2, j2);

		if (l1 != -1 && l2 != -1)
		{
			l = l1*(l1-1)/2 + l2;
			for (j = 1; j < jmax; j++) {
				input >> coeff[l][j];
			}
		}
		else {
			for (j = 1; j < jmax; j++) {
				input >> a;
			}
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
// The class calculates collisional rates
//

co_collisions::co_collisions(const string &data_path, const energy_diagram* co_di, int verbosity)
{
	bool is_ortho_h2;
	if (verbosity) 
		cout << "CO collisional rate coefficients are being initializing..." << endl;
	
	nb_lev = co_di->nb_lev;

	coll_data.push_back( new co_he_coll_data(data_path, co_di, verbosity) );
	coll_data.push_back( new co_h2_coll_data(data_path, co_di, is_ortho_h2 = false, verbosity) );
	coll_data.push_back( new co_h2_coll_data(data_path, co_di, is_ortho_h2 = true, verbosity) );
	coll_data.push_back( new co_h_coll_data(data_path, co_di, verbosity) );
 	coll_data.push_back( new co_h_vibr_coll_data(data_path, co_di, verbosity) );
	nb1 = (int) coll_data.size();

	// the data on electron collisions must be here;
	nb2 = (int) coll_data.size();

	// the data on H+ collisions must be here;
	nb3 = (int) coll_data.size();

	max_temp = new double [nb3];
	for (int i = 0; i < nb3; i++) {
		max_temp[i] = coll_data[i]->get_max_temp();
	}
}

void co_collisions::set_gas_param(double temp_neutrals, double temp_el, double he_conc, double ph2_conc, double oh2_conc, 
	double h_conc, double el_conc, double *&concentration, int *&indices) const
{
	collisional_transitions::set_gas_param(temp_neutrals, temp_el, he_conc, ph2_conc, oh2_conc, h_conc, el_conc, concentration, 
		indices);

	concentration[0] = he_conc;
	concentration[1] = ph2_conc;
	concentration[2] = oh2_conc;
	concentration[3] = h_conc;
	concentration[4] = h_conc;
}

// The energy of the first level is higher, up_lev.nb > low_lev.nb;
void co_collisions::get_rate_neutrals(const energy_level &up_lev, const energy_level &low_lev, double &down_rate, 
	double &up_rate, double temp_neutrals, const double *concentration, const int *indices) const
{
	up_rate = down_rate = 0.;
	if (up_lev.v == 0 && low_lev.v == 0) 
	{
		if (up_lev.nb < coll_data[0]->nb_lev) {
			down_rate = coll_data[0]->get_rate(up_lev.nb, low_lev.nb, indices[0], (temp_neutrals < max_temp[0]) ? temp_neutrals : max_temp[0]) *concentration[0];
		}
	
		down_rate += coll_data[1]->get_rate(up_lev.nb, low_lev.nb, indices[1], (temp_neutrals < max_temp[1]) ? temp_neutrals : max_temp[1]) *concentration[1] 
			+ coll_data[2]->get_rate(up_lev.nb, low_lev.nb, indices[2], (temp_neutrals < max_temp[2]) ? temp_neutrals : max_temp[2]) *concentration[2];
	
		down_rate += coll_data[3]->get_rate(up_lev.nb, low_lev.nb, indices[3], (temp_neutrals < max_temp[3]) ? temp_neutrals : max_temp[3]) *concentration[3];
	}
	else {
		down_rate += coll_data[4]->get_rate(up_lev.nb, low_lev.nb, indices[4], (temp_neutrals < max_temp[4]) ? temp_neutrals : max_temp[4]) *concentration[4];
	}

	if (down_rate > MIN_COLLISION_RATE)
		up_rate = down_rate *exp((low_lev.energy - up_lev.energy)*CM_INVERSE_TO_KELVINS/temp_neutrals) *up_lev.g /((double) low_lev.g);
	else down_rate = 0.;
}

//
// Functions
//

void merge_co_h_coll_data(const std::string & path)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, n, m, v1, v2, j1, j2, k, l, imax, nb_temp2, nb_lines, nb_temp, nb_lev;
	double a, mass;
	double *temp_arr, *temp_arr2, **coeff_arr;

	string fname;
	ifstream input;
	ofstream output;

	nb_lev = 132; // maximal level considered: v = 2, j = 30
	imax = nb_lev*(nb_lev-1)/2;

	mass = 28.*ATOMIC_MASS_UNIT;
	molecule co_mol("co", 1, mass, 0.);

	co_diagram *co_di =
		new co_diagram(path, co_mol, nb_lev);
	
	fname = path + "coll_co/coll_co_h_vibr1-0.txt";
	input.open(fname.c_str(), ios_base::in);

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> nb_lines >> nb_temp;

	temp_arr = new double [nb_temp];
	coeff_arr = alloc_2d_array<double>(imax, nb_temp);
	memset(*coeff_arr, 0, imax*nb_temp*sizeof(double));

	for(j = 0; j < nb_temp; j++) {
		input >> temp_arr[j];
	}
	for (i = 0; i < nb_lines; i++)
	{
		input >> v1 >> j1 >> v2 >> j2;
		
		k = co_di->get_nb(v1, j1);
		l = co_di->get_nb(v2, j2);

		if (k != -1 && l != -1)
		{
			n = k*(k-1)/2 + l;
			for (j = 0; j < nb_temp; j++) {
				input >> coeff_arr[n][j];
			}
		}
		else
		{
			for (j = 0; j < nb_temp; j++) {
				input >> a;
			}
		}
	}
	input.close();

	fname = path + "coll_co/coll_co_h_vibr2-1.txt";
	input.open(fname.c_str(), ios_base::in);

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> nb_lines >> nb_temp;

	for (j = 0; j < nb_temp; j++) {
		input >> a;
	}
	for (i = 0; i < nb_lines; i++)
	{
		input >> v1 >> j1 >> v2 >> j2;
		
		k = co_di->get_nb(v1, j1);
		l = co_di->get_nb(v2, j2);

		if (k != -1 && l != -1)
		{
			n = k*(k-1)/2 + l;
			for (j = 0; j < nb_temp; j++) {
				input >> coeff_arr[n][j];
			}
		}
		else
		{
			for (j = 0; j < nb_temp; j++) {
				input >> a;
			}
		}
	}
	input.close();

	fname = path + "coll_co/coll_co_h_vibr0-0.txt";
	input.open(fname.c_str(), ios_base::in);

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	
	input >> nb_lines >> nb_temp2;
	temp_arr2 = new double [nb_temp2];

	for (j = 0; j < nb_temp2; j++) {
		input >> temp_arr2[j];
	}
	for (i = 0; i < nb_lines; i++)
	{
		input >> v1 >> j1 >> v2 >> j2;
		
		k = co_di->get_nb(0, j1);
		l = co_di->get_nb(0, j2);

		if (k != -1 && l != -1)
		{
			n = k*(k-1)/2 + l;
			j = 0;
			for (m = 0; m < nb_temp2; m++) 
			{
				if (rounding(temp_arr[j] - temp_arr2[m]) == 0) {
					input >> coeff_arr[n][j];
					j++;
				}
				else input >> a;
			}
		}
		else
		{
			for (m = 0; m < nb_temp2; m++) {
				input >> a;
			}
		}

		k = co_di->get_nb(1, j1);
		l = co_di->get_nb(1, j2);

		if (k != -1 && l != -1)
		{
			m = k*(k-1)/2 + l;
			for (j = 0; j < nb_temp; j++) {
				coeff_arr[m][j] = coeff_arr[n][j];
			}
		}

		k = co_di->get_nb(2, j1);
		l = co_di->get_nb(2, j2);

		if (k != -1 && l != -1)
		{
			m = k*(k-1)/2 + l;
			for (j = 0; j < nb_temp; j++) {
				coeff_arr[m][j] = coeff_arr[n][j];
			}
		}
	}
	input.close();

	fname = path + "coll_co/coll_co_h_vibr.txt";
	output.open(fname.c_str(), ios_base::out);

	output << "# CO-H collision rate coefficients; compilated data;" << endl
		<< "# Song et al., J. Chem. Phys. 142, 204303 (2015); Song et al., ApJ 813, p. 96, 2015; Walker et al., ApJ 811, p. 27 (2015);" << endl
		<< "# up v,j -> low v',j'; rates in cm3/s; nb of levels " << nb_lev << endl;
	output << left << setw(8) << imax << nb_temp << endl;

	output << left << setw(20) << "";
	for (j = 0; j < nb_temp; j++) {
		output << left << setw(13) << temp_arr[j];
	}
	output << endl;

	k = 1;
	l = 0;
	for (i = 0; i < imax; i++)
	{
		output << left << setw(5) << co_di->lev_array[k].v << setw(5) << rounding(co_di->lev_array[k].j)
			<< setw(5) << co_di->lev_array[l].v << setw(5) << rounding(co_di->lev_array[l].j);

		for (j = 0; j < nb_temp; j++) {
			output << left << setw(13) << coeff_arr[i][j];
		}
		if (i < imax-1)
			output << endl;

		l++;
		if (k == l) {
			k++;
			l = 0;
		}
	}
	output.close();

	delete [] temp_arr; 
	delete [] temp_arr2;
	free_2d_array(coeff_arr);
}
