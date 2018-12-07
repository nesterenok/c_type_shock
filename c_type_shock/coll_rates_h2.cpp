//
// 04.03.2017. Check for errors. Errors were found in the data file for e-H2 collisions.
// 10.09.2017. Check for errors. Error was found in the H2-e collision data.
//		Note: there is a problem when number of H2 levels is small;

#include <stdlib.h>
#include <memory.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <math.h>
#include <cfloat>
#include <sstream>

#include "constants.h"
#include "utils.h"
#include "linear_algebra.h"
#include "coll_rates_h2.h"

#define MAX_TEXT_LINE_WIDTH 3000
#define SOURCE_NAME "coll_rates_h2.cpp"
using namespace std;

//
// Classes with collisional coefficient data
//

h2_oh2_flower_data::h2_oh2_flower_data(const std::string &data_path, const energy_diagram *h2_di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, li, lf, f, nb, nb_lines;
	double **temp;

	string file_name;
	ifstream input;

	if (rounding(h2_di->mol.spin) == 0) {
		file_name = data_path + "coll_h2/coll_ph2_oh2.txt";
	}
	else if (rounding(h2_di->mol.spin) == 1) {
		file_name = data_path + "coll_h2/coll_oh2_oh2.txt";
	}
	else {
		// the file, which contain rate coefficients for both ortho- and para-H2:
		file_name = data_path + "coll_h2/coll_h2_oh2.txt";
	}
	
	input.open(file_name.c_str(), std::ios_base::in);
	if (!input) { 
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
		exit(1);
	}

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> nb_lev;
	
	jmax = 10; // the T = 0 K point is included;
	imax = nb_lev*(nb_lev-1)/2;
	
	tgrid = new double [jmax];

	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, jmax *imax *sizeof(double));

	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, jmax *imax *sizeof(double));

	temp = alloc_2d_array<double>(nb_lev*nb_lev, jmax);
	memset(*temp, 0, nb_lev*nb_lev *jmax *sizeof(double));

	nb_lines = nb_lev*nb_lev;
	tgrid[0] = 0.; // 0 K point, what one can say about rate coefficients at small T?

	for (j = 1; j < jmax; j++) {
		input >> tgrid[j];
	}
	for (i = 0; i < nb_lines; i++)
	{
		input >> li >> lf >> f >> f;
		for (j = 1; j < jmax; j++) {
			input >> temp[i][j];
		}
	}
	input.close();
	
	for (li = 1; li < nb_lev; li++) {
		for (lf = 0; lf < li; lf++)
		{
			i = li *nb_lev + lf;
			f = lf *nb_lev + li;
			nb = li*(li-1)/2 + lf;

			if (li < h2_di->nb_lev) {
				for (j = 1; j < jmax; j++)
				{
					coeff[nb][j] = 0.5*(temp[i][j] + temp[f][j] *exp( (h2_di->lev_array[li].energy - h2_di->lev_array[lf].energy) *CM_INVERSE_TO_KELVINS/tgrid[j] ) 
						*h2_di->lev_array[lf].g /h2_di->lev_array[li].g); 
				}
			}
		}
	}
	calc_coeff_deriv();
	free_2d_array<double>(temp);
	
	if (verbosity) {
		cout << "  data have been read from file " << file_name << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

h2_ph2_flower_data::h2_ph2_flower_data(const std::string &data_path, const energy_diagram *h2_di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, li, lf, f, nb, nb_lines;
	double **temp;

	string file_name;
	ifstream input;

	if (rounding(h2_di->mol.spin) == 0) {
		// 27 level data were completed by zero rows; coefficients for 28->i were calculated from i->28 and added;
		file_name = data_path + "coll_h2/coll_ph2_ph2.txt";
	}
	else if (rounding(h2_di->mol.spin) == 1) {
		// the zero row was added for the H2 transition 1 -> 1;
		file_name = data_path + "coll_h2/coll_oh2_ph2.txt";
	}
	else {
		file_name = data_path + "coll_h2/coll_h2_ph2.txt";
	}

	input.open(file_name.c_str(), std::ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
		exit(1);
	}
	
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> nb_lev;

	jmax = 10;
	imax = nb_lev*(nb_lev-1)/2;
	
	tgrid = new double [jmax];
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, jmax *imax *sizeof(double));
	
	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, jmax *imax *sizeof(double));

	temp = alloc_2d_array<double>(nb_lev*nb_lev, jmax);
	memset(*temp, 0, nb_lev*nb_lev *jmax *sizeof(double));

	nb_lines = nb_lev*nb_lev;
	tgrid[0] = 0.; // point at 0 K;

	for (j = 1; j < jmax; j++) {
		input >> tgrid[j];
	}
	for (i = 0; i < nb_lines; i++)
	{
		input >> li >> lf >> f >> f;
		for (j = 1; j < jmax; j++) {
			input >> temp[i][j];
		}
	}
	input.close();
		
	for (li = 1; li < nb_lev; li++) {
		for (lf = 0; lf < li; lf++)
		{
			i = li *nb_lev + lf;
			f = lf *nb_lev + li;
			nb = li*(li-1)/2 + lf;

			if (li < h2_di->nb_lev) {
				for (j = 1; j < jmax; j++)
				{
					coeff[nb][j] = 0.5*(temp[i][j] + temp[f][j] *exp((h2_di->lev_array[li].energy - h2_di->lev_array[lf].energy) *CM_INVERSE_TO_KELVINS/tgrid[j]) 
						*h2_di->lev_array[lf].g /h2_di->lev_array[li].g); 
				}
			}
		}
	}
	calc_coeff_deriv();
	free_2d_array<double>(temp);
	
	if (verbosity) {
		cout << "  data have been read from file " << file_name << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

// Rate coefficients involving level 24 of ortho-H2 (v=1,j) are set to zero, one can set them equal to those involving level 25 (v=3,j=7). 
// Rate coefficients involving level 27 of para-H2 (v=1,j) are set to zero, one can set them equal to those involving level 28 (v=3, j=8).
h2_he_flower_data::h2_he_flower_data(const std::string &data_path, const energy_diagram *h2_di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, li, lf, f, nb, nb_lines;
	double **temp;

	string file_name;
	ifstream input;

	if (rounding(h2_di->mol.spin) == 0) {
		file_name = data_path + "coll_h2/coll_ph2_he.txt";
	}
	else if (rounding(h2_di->mol.spin) == 1) {
		file_name = data_path + "coll_h2/coll_oh2_he.txt";
	}
	else {
		file_name = data_path + "coll_h2/coll_h2_he.txt";
	}

	input.open(file_name.c_str(), std::ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
		exit(1);
	}
	
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> nb_lev;
	
	jmax = 10;
	imax = nb_lev*(nb_lev-1)/2;
	
	tgrid = new double [jmax];
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, jmax *imax *sizeof(double));
	
	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, jmax *imax *sizeof(double));

	temp = alloc_2d_array<double>(nb_lev*nb_lev, jmax);
	memset(*temp, 0, nb_lev*nb_lev *jmax *sizeof(double));

	nb_lines = nb_lev*nb_lev;
	tgrid[0] = 0.; // additional point at 0 K;

	for (j = 1; j < jmax; j++) {
		input >> tgrid[j];
	}
	for (i = 0; i < nb_lines; i++)
	{
		input >> li >> lf >> f >> f;
		for (j = 1; j < jmax; j++) {
			input >> temp[i][j];
		}
	}
	input.close();
		
	for (li = 1; li < nb_lev; li++) {
		for (lf = 0; lf < li; lf++)
		{
			i = li *nb_lev + lf;
			f = lf *nb_lev + li;
			nb = li*(li-1)/2 + lf;

			if (li < h2_di->nb_lev) {
				for (j = 1; j < jmax; j++)
				{
					coeff[nb][j] = 0.5*(temp[i][j] + temp[f][j] *exp((h2_di->lev_array[li].energy - h2_di->lev_array[lf].energy) *CM_INVERSE_TO_KELVINS/tgrid[j]) 
						*h2_di->lev_array[lf].g /h2_di->lev_array[li].g); 
				}
			}
		}
	}
	calc_coeff_deriv();
	free_2d_array<double>(temp);
	
	if (verbosity) {
		cout << "  data have been read from file " << file_name << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

h2_h_wrathmall_data::h2_h_wrathmall_data(const std::string &data_path, const energy_diagram *h2_di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, k, li, lf, f, nb;
	double de, a, **temp;

	string file_name;
	ifstream input;

	if (rounding(h2_di->mol.spin) == 0) {
		file_name = data_path + "coll_h2/coll_ph2_h_w.txt";
	}
	else if (rounding(h2_di->mol.spin) == 1) {
		file_name = data_path + "coll_h2/coll_oh2_h_w.txt";
	}
	else {
		file_name = data_path + "coll_h2/coll_h2_h_w.txt";
	}
	
	input.open(file_name.c_str(), std::ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
		exit(1);
	}

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> nb_lev;

	jmax = 61; // point at 0 K;
	imax = nb_lev*(nb_lev-1)/2;
	
	tgrid = new double [jmax];

	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, jmax *imax *sizeof(double));
	
	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, jmax *imax *sizeof(double));

	temp = alloc_2d_array<double>(nb_lev*nb_lev, jmax);
	memset(*temp, 0, nb_lev*nb_lev *jmax *sizeof(double));

	tgrid[0] = 0.; // point at 0 K;		
	for (j = 1; j < jmax; j++) 
	{
		input >> tgrid[j];		
		for (lf = 0; lf < nb_lev; lf++) {
			for (li = 0; li < nb_lev; li++)
			{
				// the columns of the matrix correspond to the initial level of the transition and 
				// the rows correspond to the final level of the transition;
				i = li*nb_lev + lf;
				input >> temp[i][j];
				
				if (temp[i][j] < 0.) // some coefficient values in the data are < 0;
					temp[i][j] = 0.;
			}
		}
	}
	input.close();
		
	for (li = 1; li < nb_lev; li++) {
		for (lf = 0; lf < li; lf++)
		{
			i = li *nb_lev + lf;
			f = lf *nb_lev + li;
			nb = li*(li-1)/2 + lf;

			if (li < h2_di->nb_lev) {
				for (j = 1; j < jmax; j++)
				{
					coeff[nb][j] = 0.5*(temp[i][j] + temp[f][j] *exp((h2_di->lev_array[li].energy - h2_di->lev_array[lf].energy) *CM_INVERSE_TO_KELVINS/tgrid[j]) 
						*h2_di->lev_array[lf].g/h2_di->lev_array[li].g); 
				}
			}
		}
	}
	
	// the contribution of reactive channels to the collisions are added (Le Bourlot et al. 1999):
	for (li = 1; li < h2_di->nb_lev && li < nb_lev; li++) {
		for (lf = 0; lf < li; lf++)
		{
			nb = li*(li-1)/2 + lf;
			if (h2_di->lev_array[li].v == h2_di->lev_array[lf].v) {
				if (h2_di->lev_array[li].j - h2_di->lev_array[lf].j == 2) {
					for (j = 1; j < jmax; j++) {
						coeff[nb][j] += 8.e-11*exp(-3900./tgrid[j]); // Schofield, Planet. Space Sci. 15, p.643 (1967);
					}
				} 
				else if (h2_di->lev_array[li].j - h2_di->lev_array[lf].j == 1) {
					if (rounding(h2_di->lev_array[li].j)%2 == 0) {
						for (j = 1; j < jmax; j++) {
							coeff[nb][j] += 8.e-11*exp(-3900./tgrid[j]);
						}
					}
					else {
						for (j = 1; j < jmax; j++) {
							coeff[nb][j] += 2.67e-11*exp(-3900./tgrid[j]); // reduced by a factor of 3 when J_up is odd
						}
					}
				}					
			}
			else {
				for (j = 1; j < jmax; j++) {
					de = (3900. - (h2_di->lev_array[li].energy - h2_di->lev_array[lf].energy)*CM_INVERSE_TO_KELVINS)/tgrid[j];
					if (de < 0.) 
						de = 0.;
				
					if (abs(rounding(h2_di->lev_array[li].j - h2_di->lev_array[lf].j))%2 == 0) {		
						coeff[nb][j] += coeff[nb][j]*exp(-de);
					}
					else {
						f = h2_di->get_nb(h2_di->lev_array[lf].v, h2_di->lev_array[lf].j-1);
						k = h2_di->get_nb(h2_di->lev_array[lf].v, h2_di->lev_array[lf].j+1);
						
						a = 0.;
						if (f >= 0 && f < nb_lev) {
							f = li*(li-1)/2 + f;
							a = coeff[f][j];
						}
						if (k >= 0 && k < nb_lev) {	
							k = li*(li-1)/2 + k;
							a += coeff[k][j];
						}
						if (rounding(h2_di->lev_array[li].j)%2 == 0)
							coeff[nb][j] += 0.5*a*exp(-de);
						else
							coeff[nb][j] += 0.166666666666667*a*exp(-de);
					}
				}
			}
		}
	}
	calc_coeff_deriv();
	free_2d_array<double>(temp);
	
	if (verbosity) {
		cout << "  data have been read from file " << file_name << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

h2_h2_wan_data::h2_h2_wan_data(const std::string &path, const energy_diagram *h2_di, bool coll_partner_is_ortho, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, l, u, n, v, nb_lev_f, nb_lines, nb_coll_partners;
	int *lev_index_arr;
	double a;
	
	string file_name, str;
	stringstream ss;
	ifstream input;

	file_name = path + "coll_h2/coll_h2_h2_wan2018.txt";
	input.open(file_name.c_str(), std::ios_base::in);
	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
		exit(1);
	}
	// comment lines are read:
	do
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	while (text_line[0] == '!');

	// the line with molecule name is ignored
	
	do
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	while (text_line[0] == '!');
	
	// the line with molecule weight is ignored 

	do
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	while (text_line[0] == '!');

	ss.clear();
	ss.str(text_line);
	
	ss >> nb_lev_f; // nb of levels, given in the data file
	lev_index_arr = new int [nb_lev_f];
	
	nb_lev = (nb_lev_f < h2_di->nb_lev) ? nb_lev_f : h2_di->nb_lev;
	imax = nb_lev*(nb_lev-1)/2;

	do
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	while (text_line[0] == '!');
	
	for (i = 0; i < nb_lev_f; i++) {
		ss.clear();
		ss.str(text_line);
		
		ss >> l >> a >> l >> j >> v;
		lev_index_arr[i] = h2_di->get_nb(v, j);

		input.getline(text_line, MAX_TEXT_LINE_WIDTH); // note, the line is read in the cycle end
	}
	
	do
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	while (text_line[0] == '!');

	ss.clear();
	ss.str(text_line);
	ss >> nb_lines; // nb of spectroscopic lines
	
	do
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	while (text_line[0] == '!');

	for (i = 0; i < nb_lines; i++) {
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	}

	do
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	while (text_line[0] == '!');

	ss.clear();
	ss.str(text_line);
	ss >> nb_coll_partners; // nb of collisional partners

	for (n = 0; n < nb_coll_partners; n++) {
		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '!');

		// the line with data details is ignored

		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '!');

		ss.clear();
		ss.str(text_line);
		ss >> nb_lines; // nb of collisional transitions

		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '!');

		ss.clear();
		ss.str(text_line);
		
		ss >> jmax; // nb of temperature points
		jmax++; // one point for zero temperature
		
		if (n == 0) {	
			tgrid = new double [jmax];
			tgrid[0] = 0.;

			coeff = alloc_2d_array<double>(imax, jmax);
			memset(*coeff, 0, jmax *imax *sizeof(double));
	
			coeff_deriv = alloc_2d_array<double>(imax, jmax);
			memset(*coeff_deriv, 0, jmax *imax *sizeof(double));
		}

		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '!');

		ss.clear();
		ss.str(text_line);
		for (j = 1; j < jmax; j++) {	
			ss >> tgrid[j];
		}

		do
			input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		while (text_line[0] == '!');

		if ( ((n > 0) && coll_partner_is_ortho) || (n == 0 && !coll_partner_is_ortho) ) {
			for (i = 0; i < nb_lines; i++) {
				ss.clear();
				ss.str(text_line);
				
				ss >> j >> u >> l;
				u = lev_index_arr[u-1]; // the first level has index 0 in the code
				l = lev_index_arr[l-1];

				if (u >= 0 && l >= 0) {
					v = u*(u-1)/2 + l;
					for (j = 1; j < jmax; j++) {
						ss >> coeff[v][j];
					}
				}
				else {
					for (j = 1; j < jmax; j++) {
						ss >> a;
					}
				}
				input.getline(text_line, MAX_TEXT_LINE_WIDTH);
			}
		}
		else {
			for (i = 0; i < nb_lines; i++) {
				input.getline(text_line, MAX_TEXT_LINE_WIDTH);
			}
		}
	}
	input.close();

	calc_coeff_deriv();
	if (verbosity) {
		cout << "  data have been read from file " << file_name << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

h2_h_lique_data::h2_h_lique_data(const string &path, const energy_diagram *h2_di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, v, k, l;
	string file_name;
	ifstream input;
	
	file_name = path + "coll_h2/coll_h2_h_l.txt";
	input.open(file_name.c_str(), std::ios_base::in);

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
		exit(1);
	}

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	input >> nb_lev;
	jmax = 51; // +1 for zero point;
	imax = nb_lev*(nb_lev-1)/2;
	
	tgrid = new double [jmax];

	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, jmax *imax *sizeof(double));
	
	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, jmax *imax *sizeof(double));

	for (j = 0; j < jmax; j++) {
		tgrid[j] = 100.*j;
	}

	i = 0;
	while (i < imax)
	{
		input >> v >> j;
		k = h2_di->get_nb(v, j);

		input >> v >> j;
		l = h2_di->get_nb(v, j);

		l = k*(k-1)/2 + l; // it is assumed that k > l;
		for (j = 1; j < jmax; j++) {
			input >> coeff[l][j];
		}
		i++;
	}
	input.close();

	calc_coeff_deriv();
	if (verbosity) {
		cout << "  data have been read from file " << file_name << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

h2_h_bossion_data::h2_h_bossion_data(const string &path, const energy_diagram *h2_di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, k, l, j, n, vi, ji, vf, jf, nb_lines;
	double a;
	string file_name;
	ifstream input;
	
	file_name = path + "coll_h2/coll_h2_h_b.txt";
	input.open(file_name.c_str(), std::ios_base::in);

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
		exit(1);
	}
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	
	input >> nb_lines >> jmax; // check in the file the presence of these parameters;
	jmax++; // +1 for zero point;

	nb_lev = h2_di->nb_lev;
	imax = nb_lev*(nb_lev-1)/2;
	
	tgrid = new double [jmax];
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, jmax *imax *sizeof(double));
	
	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, jmax *imax *sizeof(double));

	tgrid[0] = 0.;
	for (j = 1; j < jmax; j++) {
		input >> tgrid[j];
	}
	
	n = 0;
	while (!input.eof() && n < nb_lines)
	{
		input >> vi >> ji >> vf >> jf;
		l = h2_di->get_nb(vi, ji); // initial level
		k = h2_di->get_nb(vf, jf); // final

		if (l > -1 && k > -1)
		{
			i = l*(l-1)/2 + k;
			for (j = 1; j < jmax; j++) 
			{
				input >> a;
				coeff[i][j] = a;
			}
		}
		else {
			for (j = 1; j < jmax; j++) {
				input >> a;
			}
		}
		n++;
	}
	input.close();
	
	calc_coeff_deriv();
	if (verbosity) {
		cout << "  data have been read from file " << file_name << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

// Maximal temperature for which coefficients are calculated is 10000 K;
h2_h_martin_data::h2_h_martin_data(const string &data_path, const energy_diagram *h2_di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, li, lf, nb_lines;
	double c1, c2, c3, c4, d, t, g;

	string file_name;
	ifstream input;

	nb_lev = h2_di->nb_lev;
	imax = nb_lev*(nb_lev-1)/2;
	jmax = 61;

	tgrid = new double [jmax];
	
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, imax*jmax*sizeof(double));
	
	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, jmax *imax *sizeof(double));

	tgrid[0] = 0.;
	for (j = 1; j < 21; j++) {
		tgrid[j] = tgrid[j-1] + 100.;
	}
	for (j = 21; j < jmax; j++) {
		tgrid[j] = tgrid[j-1] + 200.;
	}
	
	file_name = data_path + "coll_h2/coll_h2_h_m.txt";
	input.open(file_name.c_str(), std::ios_base::in);

	if (!input.is_open()) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
		exit(1);
	}
	
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> nb_lines;

	for (i = 0; i < nb_lines && i < imax; i++)
	{
		input >> li >> lf;
		input >> c1 >> c2 >> c3 >> c4;
		
		if (abs(rounding(h2_di->lev_array[li-1].j - h2_di->lev_array[lf-1].j))%2 == 0)
			d = 1.;
		else d = 0.3;

		for (j = 1; j < jmax; j++) 
		{	
			if (tgrid[j] > 50000.)
				t = d + 50.;
			else t = d + 0.001*tgrid[j];

			g = c1 + c2/t + c3/(t*t) + c4*t;
			coeff[i][j] = pow(10., g);
		}
	}
	input.close();

	calc_coeff_deriv();
	if (verbosity) {
		cout << "  data have been read from file " << file_name << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

h2_e_data::h2_e_data(const string &data_path, const energy_diagram *h2_di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, l, f, nb_lev_f;
	double a;

	string file_name;
	ifstream input;

	file_name = data_path + "coll_h2/coll_h2_e.txt";
	input.open(file_name.c_str(), std::ios_base::in);

	if (!input.is_open()) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
		exit(1);
	}
	
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	
	input >> nb_lev_f >> jmax;
	nb_lev = (nb_lev_f < h2_di->nb_lev) ? nb_lev_f : h2_di->nb_lev;

	jmax++; // temperature value at 0 K is added;
	imax = nb_lev*(nb_lev-1)/2;
	
	tgrid = new double [jmax];
	
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, imax*jmax*sizeof(double));
	
	coeff_deriv = alloc_2d_array<double>(imax, jmax);
	memset(*coeff_deriv, 0, jmax *imax *sizeof(double));

	tgrid[0] = 0.;
	for (j = 1; j < jmax; j++) {
		input >> tgrid[j];
	}

	l = 0;
	for (i = 0; i < nb_lev_f; i++) {
		for (f = 0; f < i; f++)
		{
			input >> j >> j >> j >> j;
			if (i < nb_lev && f < nb_lev)
			{
				for (j = 1; j < jmax; j++) {
					input >> coeff[l][j];
				}
				l++;
			}
			else {
				for (j = 1; j < jmax; j++) {
					input >> a;
				}
			}
		}
	}
	input.close();

	calc_coeff_deriv();
	if (verbosity) {
		cout << "  data have been read from file " << file_name << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

//
// The class that calculates collisional rates
//

h2_collisions::h2_collisions(const std::string &data_path, const energy_diagram *h2_di, int verbosity)
	: collisional_transitions()
{
	bool coll_partner_is_ortho;
	if (verbosity) {
		cout << "H2 collisional rate coefficients are being initializing..." << endl;
	}
	nb_lev = h2_di->nb_lev;

	coll_data.push_back( new h2_he_flower_data(data_path, h2_di, verbosity) );

	coll_data.push_back( new h2_h2_wan_data(data_path, h2_di, coll_partner_is_ortho = false, verbosity) );
	coll_data.push_back( new h2_ph2_flower_data(data_path, h2_di, verbosity) );
	
	coll_data.push_back( new h2_h2_wan_data(data_path, h2_di, coll_partner_is_ortho = true, verbosity) );
	coll_data.push_back( new h2_oh2_flower_data(data_path, h2_di, verbosity) );

	coll_data.push_back( new h2_h_lique_data(data_path, h2_di, verbosity) );
	coll_data.push_back( new h2_h_bossion_data(data_path, h2_di, verbosity) );	
	coll_data.push_back( new h2_h_wrathmall_data(data_path, h2_di, verbosity) );
	coll_data.push_back( new h2_h_martin_data(data_path, h2_di, verbosity) );

	nb1 = (int) coll_data.size();

	coll_data.push_back( new h2_e_data(data_path, h2_di, verbosity) );
	nb2 = (int) coll_data.size();

	// data on H+ collisions must be here;
	nb3 = (int) coll_data.size();

	max_temp = new double [nb3];
	for (int i = 0; i < nb3; i++) {
		max_temp[i] = coll_data[i]->get_max_temp();
	}
}

void h2_collisions::check_spline(int ilev, int flev, const std::string & fname) const
{
#if H2_COLL_CUBIC_SPLINE
	dynamic_cast<collision_data_cub_spline*>(coll_data[0])->check_spline(ilev, flev, fname);
#endif
}

void h2_collisions::set_gas_param(double temp_neutrals, double temp_el, double he_conc, double ph2_conc, double oh2_conc, 
	double h_conc, double el_conc, double *&concentration, int *&indices) const
{
	// must be called first:
	collisional_transitions::set_gas_param(temp_neutrals, temp_el, he_conc, ph2_conc, oh2_conc, h_conc, el_conc, concentration, 
		indices);

	concentration[0] = he_conc;
	concentration[1] = concentration[2] = ph2_conc;
	concentration[3] = concentration[4] = oh2_conc;
	concentration[5] = concentration[6] = concentration[7] = concentration[8] = h_conc;
}

// The energy of the first level is higher, up_lev.nb > low_lev.nb;
void h2_collisions::get_rate_neutrals(const energy_level &up_lev, const energy_level &low_lev, double &down_rate, 
	double &up_rate, double temp_neutrals, const double *concentration, const int *indices) const
{
	down_rate = 0.;
	if (rounding(low_lev.spin - up_lev.spin) == 0) {
		// collisions with He
		if (up_lev.nb < coll_data[0]->nb_lev) {
			down_rate += coll_data[0]->get_rate(up_lev.nb, low_lev.nb, indices[0], (temp_neutrals < max_temp[0]) ? temp_neutrals : max_temp[0]) 
				*concentration[0];	
		}
		// collisions with p-H2
		if (low_lev.v == 0 && up_lev.v == 0 && up_lev.nb < coll_data[1]->nb_lev) {
			down_rate += coll_data[1]->get_rate(up_lev.nb, low_lev.nb, indices[1], (temp_neutrals < max_temp[1]) ? temp_neutrals : max_temp[1]) 
				*concentration[1];	
		}
		else if (up_lev.nb < coll_data[2]->nb_lev) {
			down_rate += coll_data[2]->get_rate(up_lev.nb, low_lev.nb, indices[2], (temp_neutrals < max_temp[2]) ? temp_neutrals : max_temp[2]) 
				*concentration[2];	
		}
		// collisions with o-H2
		if (low_lev.v == 0 && up_lev.v == 0 && up_lev.nb < coll_data[3]->nb_lev) {
			down_rate += coll_data[3]->get_rate(up_lev.nb, low_lev.nb, indices[3], (temp_neutrals < max_temp[3]) ? temp_neutrals : max_temp[3]) 
				*concentration[3];	
		}
		else if (up_lev.nb < coll_data[4]->nb_lev) {
			down_rate += coll_data[4]->get_rate(up_lev.nb, low_lev.nb, indices[4], (temp_neutrals < max_temp[4]) ? temp_neutrals : max_temp[4]) 
				*concentration[4];	
		}
	}
#if (H2_H_COLL_DATA == 0)
	// H2-H data by Wrathmall et al. (2007):
	if (up_lev.nb < coll_data[7]->nb_lev) {
		down_rate += coll_data[7]->get_rate(up_lev.nb, low_lev.nb, indices[7], (temp_neutrals < max_temp[7]) ? temp_neutrals : max_temp[7]) 
			*concentration[7];
	}
	// H2-H data by Mandy & Martin (1995):
	else if (up_lev.nb < coll_data[8]->nb_lev) {
		down_rate += coll_data[8]->get_rate(up_lev.nb, low_lev.nb, indices[8], (temp_neutrals < max_temp[8]) ? temp_neutrals : max_temp[8]) 
			*concentration[8];
	}
#elif (H2_H_COLL_DATA == 1)
	// H2-H data by Lique (2015) is used:
	if (up_lev.nb < coll_data[5]->nb_lev) {
		down_rate += coll_data[5]->get_rate(up_lev.nb, low_lev.nb, indices[5], (temp_neutrals < max_temp[5]) ? temp_neutrals : max_temp[5]) 
			*concentration[5];
	}
	// H2-H data by Wrathmall et al. (2007):
	else if (up_lev.nb < coll_data[7]->nb_lev) {
		down_rate += coll_data[7]->get_rate(up_lev.nb, low_lev.nb, indices[7], (temp_neutrals < max_temp[7]) ? temp_neutrals : max_temp[7]) 
			*concentration[7];
	}
	// H2-H data by Mandy & Martin (1995):
	else if (up_lev.nb < coll_data[8]->nb_lev) {
		down_rate += coll_data[8]->get_rate(up_lev.nb, low_lev.nb, indices[8], (temp_neutrals < max_temp[8]) ? temp_neutrals : max_temp[8]) 
			*concentration[8];
	}
#elif (H2_H_COLL_DATA >= 2)
	// H2-H data by Lique (2015) is used:
	if (up_lev.nb < coll_data[5]->nb_lev) {
		down_rate += coll_data[5]->get_rate(up_lev.nb, low_lev.nb, indices[5], (temp_neutrals < max_temp[5]) ? temp_neutrals : max_temp[5]) 
			*concentration[5];
	}
	// H2-H data by Bossion et al. MNRAS 480, p.3718, 2018;
	else if (up_lev.nb < coll_data[6]->nb_lev) {
		down_rate += coll_data[6]->get_rate(up_lev.nb, low_lev.nb, indices[6], (temp_neutrals < max_temp[6]) ? temp_neutrals : max_temp[6]) 
			*concentration[6];
	}
#endif

	if (down_rate > MIN_COLLISION_RATE) {
		up_rate = down_rate *exp((low_lev.energy - up_lev.energy)*CM_INVERSE_TO_KELVINS/temp_neutrals) *up_lev.g /low_lev.g;
	}
	else down_rate = up_rate = 0.;
}

//
// H2 dissociation
//

h2_h_dissociation_data::h2_h_dissociation_data(const std::string & path, const energy_diagram *h2_di, int verbosity)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, n, vi, ji, nb_lines;
	double a;
	
	string file_name;
	ifstream input;
	
	file_name = path + "coll_h2/diss_h2_h_b.txt";
	input.open(file_name.c_str(), std::ios_base::in);

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
		exit(1);
	}
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	
	input >> nb_lines >> jmax; // check in the file the presence of these parameters;
	jmax++; // +1 for zero point;
	nb_lev = h2_di->nb_lev;
	
	tgrid = new double [jmax];
	coeff = alloc_2d_array<double>(nb_lev, jmax);
	memset(*coeff, 0, nb_lev*jmax*sizeof(double));
	
	coeff_deriv = alloc_2d_array<double>(nb_lev, jmax);
	memset(*coeff_deriv, 0, nb_lev*jmax*sizeof(double));

	tgrid[0] = 0.;
	for (j = 1; j < jmax; j++) {
		input >> tgrid[j];
	}
	
	n = 0;
	while (!input.eof() && n < nb_lines)
	{
		input >> vi >> ji;
		i = h2_di->get_nb(vi, ji);
	
		if (i > -1) {
			for (j = 1; j < jmax; j++) {
				input >> coeff[i][j];
			}	
		}
		else {
			for (j = 1; j < jmax; j++) {
				input >> a;
			}
		}
		n++;
	}
	input.close();
	
	calc_coeff_deriv();
	if (verbosity) {
		cout << "  data have been read from file " << file_name << endl
			<< "  temperature range " << (int) tgrid[1] << " - " << (int) tgrid[jmax-1] << endl;
	}
}

//
// H2 excitation processes
//

h2_grain_formation::h2_grain_formation(const energy_diagram *h2_di)
{
	int i;
	double sum = 0.;

	nb_lev = h2_di->nb_lev;
	popul = new double [nb_lev];
	average_energy = 4.4781 *EV_TO_CM_INVERSE/3.; // dissociation energy of H2 - 36118.11 cm-1 = 4.4781 eV  

	for (i = 0; i < nb_lev; i++) 
	{
		popul[i] = h2_di->lev_array[i].g *exp(-h2_di->lev_array[i].energy/average_energy);
		sum += popul[i];
	}
	for (i = 0; i < nb_lev; i++) {
		popul[i] /= sum;
	}
}

h2_grain_formation::~h2_grain_formation() {
	delete [] popul;
}

h2_gasphase_formation::h2_gasphase_formation(const energy_diagram *h2_di)
{
	int i, j;
	double gtemp, sum = 0.;

	max_temp = 10000.; // in K
	step = 100.; // in K
	step_inv = 1./step;

	nb_lev = h2_di->nb_lev;
	nb_temp = (int) (max_temp/step) + 1;
	
	popul = alloc_2d_array<double>(nb_lev, nb_temp);
	memset(*popul, 0, nb_lev*nb_temp*sizeof(double));
	
	popul[0][0] = 1.; // zero temperature distribution;
	for (j = 1; j < nb_temp; j++) {
		gtemp = step *j/CM_INVERSE_TO_KELVINS; // in cm-1
		for (i = 0; i < nb_lev; i++) 
		{
			popul[i][j] = h2_di->lev_array[i].g *exp(-h2_di->lev_array[i].energy/gtemp);
			sum += popul[i][j];
		}
		for (i = 0; i < nb_lev; i++) {
			popul[i][j] /= sum;
		}
	}
}

double h2_gasphase_formation::get_efficiency(int l, double gtemp) const
{
	int j;
	double x;

	if (gtemp > max_temp) 
		gtemp = max_temp;
	
	x = step_inv*gtemp;
	j = (int) x;
	x = x - j;
	return (1 - x)*popul[l][j] + x*popul[l][j+1];
}

h2_gasphase_formation::~h2_gasphase_formation() {
	free_2d_array(popul);
}

h2_excit_cosmic_rays::h2_excit_cosmic_rays(const std::string &data_path, const energy_diagram *h2_di, int verbosity) 
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, l, v, li, lf, dj, imax, jmax;
	double a;

	string file_name = data_path + "xray_efficiencies.txt";
	ifstream input(file_name.c_str());
	
	nb_ion = 4;
	ion_grid = new double [nb_ion];
	
	nb_init_lev = h2_di->get_nb(0, 9.) + 1; // be carefull here, the level (v,j) = (1,0) is below (0,9);
	nb_fin_lev = h2_di->nb_lev;
	
	eff = alloc_2d_array<double>(nb_ion, nb_fin_lev*nb_init_lev);
	memset(*eff, 0, nb_ion*nb_fin_lev*nb_init_lev*sizeof(double));

	exit_eff = alloc_2d_array<double>(nb_ion, nb_init_lev);
	memset(*exit_eff, 0, nb_ion*nb_init_lev*sizeof(double));
	
	if (input) 
	{
		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input >> imax >> jmax;

		for (l = 0; l < nb_ion; l++)
		{
			input >> ion_grid[l];
			for (i = 0; i < imax; i++)
			{
				input >> v >> dj;
				for (j = 0; j < jmax; j++) 
				{
					input >> a;		
					li = h2_di->get_nb(0, (double) j);
					lf = h2_di->get_nb(v, (double) (j+dj));

					if (li != -1 && lf != -1)
					{
						eff[l][li*nb_fin_lev + lf] = a;
						exit_eff[l][li] += a;
					}
				}
			}
		}
	}
	else {
		cout << "Error in " << SOURCE_NAME << ": can't open file " << file_name << endl;
		exit(1);
	}
	input.close();

	if (verbosity) 
		cout << "  data have been read from file " << file_name << endl;
}

h2_excit_cosmic_rays::~h2_excit_cosmic_rays()
{
	delete [] ion_grid;
	free_2d_array<double>(eff);
	free_2d_array<double>(exit_eff);
}

// The maximal value of the index is nb_ion-2; the minimal value is -1;
// normal-log approximation is used;
void h2_excit_cosmic_rays::set_ionization(double ionization, int &index, double &param) const
{ 
	int i;
	for (i = 0; i < nb_ion-1; i++) {
		if (ionization < ion_grid[i])
			break;
	}
	i--;

	if (i >= 0)
		param = log(ionization/ion_grid[i]) /log(ion_grid[i+1]/ion_grid[i]);
	else param = 0.;
	index = i;
}

double h2_excit_cosmic_rays::get_efficiency(int il, int fl, int index, double param) const
{
	if (il > nb_init_lev-1) return 0.;
	double answ;

	if (index < 0) {
		answ = eff[0][il*nb_fin_lev+fl];
	}
	else if (index < nb_ion-1) {
		answ = eff[index][il*nb_fin_lev + fl] + param *(eff[index+1][il*nb_fin_lev + fl] - eff[index][il*nb_fin_lev + fl]);
	}
	else // the case of large ionization fractions, > 0.01, must be considered accurately;
		answ = 0.;
	
	return answ;
}

double h2_excit_cosmic_rays::get_efficiency(int l, int index, double param) const
{
	if (l > nb_init_lev-1) return 0.;
	
	double answ;
	if (index < 0) {
		answ = exit_eff[0][l];
	}
	else if (index < nb_ion-1) {
		answ = exit_eff[index][l] + param*(exit_eff[index+1][l] - exit_eff[index][l]);
	}
	else // the case of large ionization fractions, > 0.01, must be considered accurately;
		answ = 0.;

	return answ;
}

//
// Functions
//

void h2_coll_data_process(const std::string &data_path)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, f, k, l, m, jmax, nb_lev_h2, nb_lev_ph2, nb_lev_oh2;
	double *tgrid, **temp;

	string fnp, fno, fn;
	string n1_arr[3] = {"coll_ph2_oh2.txt", "coll_ph2_ph2.txt", "coll_ph2_he.txt"};
	string n2_arr[3] = {"coll_oh2_oh2.txt", "coll_oh2_ph2.txt", "coll_oh2_he.txt"};
	string n3_arr[3] = {"coll_h2_oh2.txt", "coll_h2_ph2.txt", "coll_h2_he.txt"};
	
	ifstream input;
	ofstream out;
	nb_lev_h2 = nb_lev_ph2 = nb_lev_oh2 = 320;

	molecule ph2_mol("ph2", 1, 2.*ATOMIC_MASS_UNIT, 0);
	molecule oh2_mol("oh2", 1, 2.*ATOMIC_MASS_UNIT, 1);
	molecule h2_mol("h2", 1, 2.*ATOMIC_MASS_UNIT);
	 
	h2_diagram *ph2_di
		= new h2_diagram(data_path, ph2_mol, nb_lev_ph2);

	h2_diagram *oh2_di 
		= new h2_diagram(data_path, oh2_mol, nb_lev_oh2);
	
	h2_diagram *h2_di 
		= new h2_diagram(data_path, h2_mol, nb_lev_h2);

	// Nb of temperatures for the data in question:
	jmax = 9;
	tgrid = new double [jmax];
	temp = alloc_2d_array<double>(nb_lev_h2*nb_lev_h2, jmax);
	
	for (m = 0; m < 3; m++) 
	{
		memset(*temp, 0, nb_lev_h2*nb_lev_h2*jmax *sizeof(double));

		fnp = data_path + "coll_h2/";
		fnp += n1_arr[m];
	
		fno = data_path + "coll_h2/";
		fno += n2_arr[m];

		input.open(fnp.c_str(), std::ios_base::in);
		input.seekg (0, ios_base::beg);

		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input >> nb_lev_ph2;

		for (j = 0; j < jmax; j++) {
			input >> tgrid[j];
		}
		for (l = 0; l < nb_lev_ph2*nb_lev_ph2; l++)
		{
			input >> i >> f >> j >> j;
			// Calculation of the level numbers in the total H2 level list:
			i = h2_di->get_nb(ph2_di->lev_array[i-1].v, ph2_di->lev_array[i-1].j);
			f = h2_di->get_nb(ph2_di->lev_array[f-1].v, ph2_di->lev_array[f-1].j);

			for (j = 0; j < jmax; j++) {	
				input >> temp[i*nb_lev_h2+f][j];
			}
		}
		input.close();

		input.open(fno.c_str(), std::ios_base::in);
		input.seekg (0, ios_base::beg);

		input.getline(text_line, MAX_TEXT_LINE_WIDTH);
		input >> nb_lev_oh2;

		for (j = 0; j < jmax; j++) {
			input >> tgrid[j];
		}
		for (l = 0; l < nb_lev_oh2*nb_lev_oh2; l++)
		{
			input >> i >> f >> j >> j;
			i = h2_di->get_nb(oh2_di->lev_array[i-1].v, oh2_di->lev_array[i-1].j);
			f = h2_di->get_nb(oh2_di->lev_array[f-1].v, oh2_di->lev_array[f-1].j);

			for (j = 0; j < jmax; j++) {	
				input >> temp[i*nb_lev_h2+f][j];
			}
		}
		input.close();

		// Calculation of the number of the highest level, for which rate coefficients are available:
		i = h2_di->get_nb(ph2_di->lev_array[nb_lev_ph2-1].v, ph2_di->lev_array[nb_lev_ph2-1].j);
		j = h2_di->get_nb(oh2_di->lev_array[nb_lev_oh2-1].v, oh2_di->lev_array[nb_lev_oh2-1].j);
		k = (i > j ) ? i+1 : j+1;
	
		fn = data_path + "coll_h2/";
		fn += n3_arr[m];

		out.open(fn.c_str());
		out << "# joined data for o-H2 and p-H2;" << endl << k << endl;

		for (j = 0; j < jmax; j++) {
			out << left << setw(12) << tgrid[j];
		}
		out << endl;
	
		for (i = 0; i < k; i++) {
			for (f = 0; f < k; f++) 
			{
				out << left << setw(5) << i+1 << setw(5) << f+1 << setw(5) << "1" << setw(5) << "1";
				for (j = 0; j < jmax; j++) {
					 out << left << setw(12) << temp[i*nb_lev_h2+f][j];
				}
				if (i < k-1 || f < k-1) 
					out << endl;
			}
		}
		out.close();
	}

	delete [] tgrid;
	free_2d_array(temp);

	jmax = 60;
	tgrid = new double [jmax];

	temp = alloc_2d_array<double>(nb_lev_h2*nb_lev_h2, jmax);
	memset(*temp, 0, nb_lev_h2*nb_lev_h2*jmax *sizeof(double));

	fnp = data_path + "coll_h2/coll_ph2_h_w.txt";
	fno = data_path + "coll_h2/coll_oh2_h_w.txt";
	
	input.open(fnp.c_str(), std::ios_base::in);
	input.seekg(0, ios_base::beg);

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> nb_lev_ph2;

	for (j = 0; j < jmax; j++) 
	{
		input >> tgrid[j];
		for (f = 0; f < nb_lev_ph2; f++) {
			for (i = 0; i < nb_lev_ph2; i++)
			{
				l = h2_di->get_nb(ph2_di->lev_array[i].v, ph2_di->lev_array[i].j); // initial level
				m = h2_di->get_nb(ph2_di->lev_array[f].v, ph2_di->lev_array[f].j);
				input >> temp[l*nb_lev_h2+m][j];
			}
		}
	}
	input.close();

	input.open(fno.c_str(), std::ios_base::in);
	input.seekg (0, ios_base::beg);

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> nb_lev_oh2;

	for (j = 0; j < jmax; j++) 
	{
		input >> tgrid[j];
		for (f = 0; f < nb_lev_oh2; f++) {
			for (i = 0; i < nb_lev_oh2; i++)
			{
				l = h2_di->get_nb(oh2_di->lev_array[i].v, oh2_di->lev_array[i].j);
				m = h2_di->get_nb(oh2_di->lev_array[f].v, oh2_di->lev_array[f].j);
				input >> temp[l*nb_lev_h2+m][j];
			}
		}
	}
	input.close();

	i = h2_di->get_nb(ph2_di->lev_array[nb_lev_ph2-1].v, ph2_di->lev_array[nb_lev_ph2-1].j);
	j = h2_di->get_nb(oh2_di->lev_array[nb_lev_oh2-1].v, oh2_di->lev_array[nb_lev_oh2-1].j);
	k = (i > j ) ? i+1 : j+1; // 
	
	fn = data_path + "coll_h2/coll_h2_h_w.txt";	
	out.open(fn.c_str());
	out.precision(4);

	out << "# S.A. Wrathmall et al., Mon. Not. R. Astron. Soc. 382, 133 (2007); http://ccp7.dur.ac.uk/pubs.html" << endl
		<< "# columns correspond to the initial level, rows - to the final level; joined data for o-H2 and p-H2;" << endl << k << endl;
	for (j = 0; j < jmax; j++) 
	{
		out << tgrid[j] << endl;
		for (f = 0; f < k; f++) {
			for (i = 0; i < k; i++) {
				out << left << setw(12) << temp[i*nb_lev_h2+f][j];
			}
			if (f < k-1 || j < jmax-1) 
				out << endl;
		}
	}
	out.close();

	delete ph2_di;
	delete oh2_di;
	delete h2_di;

	delete [] tgrid;
	free_2d_array(temp);
}

void h2_coll_data_process_bossion(const std::string &data_path)
{
	int k, vi, vf, ji, jf, nbt, nb_lev;
	double r;
	string fname;
	ifstream input;
	ofstream output;

	nb_lev = 318; //
	molecule h2_mol("H2", 1, 2.*ATOMIC_MASS_UNIT);
	
	h2_diagram *h2_di 
		= new h2_diagram(data_path, h2_mol, nb_lev);
 
	// please, check the data format, presence of the comments, nb of temperature values;
	nbt = 100;
	fname = data_path + "coll_h2/bossion_2018_h_h2_sts_rates_100-10000.txt";
	input.open(fname.c_str());

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << fname << endl;
		exit(1);
	}

	fname = data_path + "coll_h2/coll_h2_h_b_.txt";
	output.open(fname.c_str());
	
	output << scientific;
	output.precision(2);
	output << "# The table contains the H2-H collisional rate coefficients. Bossion et al. MNRAS 480, p.3718, 2018;" << endl
		<< "# vi ji -> vf jf, k(cm3 s-1) (T); nb of data lines, nb of temperatures (T in K):" << endl;

	output << nbt << endl << left << setw(16) << "";
	for (k = 1; k <= nbt; k++) {
		output << left << setw(10) << k*100;
	}
	
	while (!input.eof())
	{
		input >> vi >> ji >> vf >> jf;
		if (input.eof())
			break;

		if (h2_di->get_nb(vi, ji) != -1 && h2_di->get_nb(vf, jf) != -1 && h2_di->get_nb(vi, ji) > h2_di->get_nb(vf, jf)) {
			output << endl << left << setw(4) << vi << setw(4) << ji << setw(4) << vf << setw(4) << jf;

			for (k = 0; k < nbt; k++) {
				input >> r;
				r = (r > 0.) ? r : 0.;
				output << left << setw(10) << r;
			}
		}
		else {
			for (k = 0; k < nbt; k++) {
				input >> r;
			}
		}
	}
	input.close();
	output.close();

	// please, check the data format, presence of the comments, nb of temperature values;
	nbt = 100;
	fname = data_path + "coll_h2/bossion_2018_h_h2_diss_100-10000.txt";
	input.open(fname.c_str());

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with dissociation data " << fname << endl;
		exit(1);
	}

	fname = data_path + "coll_h2/diss_h2_h_b.txt";
	output.open(fname.c_str());
	
	output << scientific;
	output.precision(3);
	output << "# The table contains the H2+H->H+H+H dissociation rate coefficients. Bossion et al. MNRAS 480, p.3718, 2018;" << endl
		<< "# vi ji, k(cm3 s-1) (T), temperature in K; " << endl;

	output << nbt << endl << left << setw(8) << "";
	for (k = 1; k <= nbt; k++) {
		output << left << setw(11) << k*100;
	}
	
	while (!input.eof())
	{
		input >> vi >> ji;
		if (input.eof())
			break;

		output << endl << left << setw(4) << vi << setw(4) << ji;

		for (k = 0; k < nbt; k++) {
			input >> r;
			r = (r > 0.) ? r : 0.;
			output << left << setw(11) << r;
		}
	}
	input.close();
	output.close(); 
}

void h2_coll_data_check(const std::string &data_path)
{
	char text_line[MAX_TEXT_LINE_WIDTH];
	int i, j, l, v, k, imax, jmax, nb_lev;
	int *indices(0);
	double o, p, tot_h_conc;
	double **coeff, *tgrid, *rates, *pop, *concentration(0);
	
	string file_name;
	ifstream input;
	
	nb_lev = 54; //
	molecule h2_mol("H2", 1, 2.*ATOMIC_MASS_UNIT);
	
	h2_diagram *h2_di 
		= new h2_diagram(data_path, h2_mol, nb_lev);
 
	h2_einstein_coeff *h2_einst 
		= new h2_einstein_coeff(data_path, h2_di);
	
	h2_collisions *h2_coll 
		= new h2_collisions(data_path, h2_di);

	
	file_name = data_path + "coll_h2/coll_h2_h_l.txt";
	input.open(file_name.c_str(), std::ios_base::in);

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
		exit(1);
	}

	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);

	input >> k;
	jmax = 51; // +1 for zero point;
	imax = nb_lev*(nb_lev-1)/2;
	
	tgrid = new double [jmax];
	coeff = alloc_2d_array<double>(imax, jmax);
	memset(*coeff, 0, jmax *imax *sizeof(double));
	
	for (j = 0; j < jmax; j++) {
		tgrid[j] = 100.*j;
	}
	
	for (l = 0; l < imax; l++)
	{
		input >> v >> j;
		i = h2_di->get_nb(v, j);

		input >> v >> j;
		k = h2_di->get_nb(v, j);

		if (i <= k)
			cout << "!";

		i = i*(i-1)/2 + k;
		for (j = 1; j < jmax; j++) {
			input >> coeff[i][j];
		}
	}
	input.close();

	rates = new double [nb_lev];
	memset(rates, 0, nb_lev*sizeof(double));

	pop = new double [nb_lev];
	memset(pop, 0, nb_lev*sizeof(double));

	j = 44;
	tot_h_conc = 2.45e+4;
//	boltzmann_populations(tgrid[j], pop, h2_di);
	
	file_name = "h2_pop_test.txt";
	input.open(file_name.c_str(), std::ios_base::in);

	if (!input) {
		cout << "Error in " << SOURCE_NAME << ": can't open file with collisional data " << file_name << endl;
		exit(1);
	}
	input.getline(text_line, MAX_TEXT_LINE_WIDTH);
	input >> k;
	
	for (i = 0; i < k && i < nb_lev; i++) {
		input >> v >> pop[i];
	}

	h2_coll->set_gas_param(tgrid[j], tgrid[j], 0.1*tot_h_conc, 0.125*tot_h_conc, 0.375*tot_h_conc, 0.002*tot_h_conc, 1.e-7*tot_h_conc,
		concentration, indices);

//	opt_thin_pop(pop, h2_di, h2_einst, h2_coll, tgrid[j], tgrid[j], concentration, indices);
/*	double down1, down2, up1, up2;
	double **matrix = alloc_2d_array<double>(nb_lev, nb_lev);
	memset(*matrix, 0, nb_lev*nb_lev*sizeof(double));
	
	for (v = 1; v < nb_lev; v++) {
		for (k = 0; k < v; k++)
		{
			h2_coll->get_rate_neutrals(h2_di->lev_array[v], h2_di->lev_array[k], down1, up1,
				tgrid[j], concentration, indices);

			h2_coll->get_rate_electrons(h2_di->lev_array[v], h2_di->lev_array[k], down2, up2,
				tgrid[j], concentration, indices);

			matrix[k][v] = h2_einst->arr[v][k] + down1 + down2;	// v->k
			matrix[v][v] -= h2_einst->arr[v][k] + down1 + down2;
		
			matrix[v][k] = up1 + up2;	// k->v
			matrix[k][k] -= up1 + up2;
		}
	}

	for (v = 0; v < nb_lev; v++) {
		matrix[0][v] = 1.;
	}
	pop[0] = 1.;
	
	lu_matrix_solve(matrix, pop, nb_lev);
	free_2d_array<double>(matrix);

*/

	j = 44;
	for (v = 1; v < nb_lev; v++) {
		for (k = 0; k < v; k++) 
		{
			if (abs(rounding(h2_di->lev_array[v].spin - h2_di->lev_array[k].spin)) > DBL_EPSILON) 
			{
				i = v*(v-1)/2 + k;
				p = (coeff[i][j] + h2_einst->arr[v][k])*pop[v]; // radiative transitions does not matter here
				rates[v] -= p;
				rates[k] += p;

				p = coeff[i][j]*pop[k] 
					*exp((h2_di->lev_array[k].energy - h2_di->lev_array[v].energy)*CM_INVERSE_TO_KELVINS/tgrid[j]) 
					*h2_di->lev_array[v].g /h2_di->lev_array[k].g;
					
				rates[v] += p;
				rates[k] -= p;
			}
		}
	}
	o = p = 0.;
	for (v = 0; v < nb_lev; v++) {
		if (rounding(h2_di->lev_array[v].spin) == 0)
			o += rates[v];
		else p += rates[v];
	}
	cout << "H2 molecule excitation by H - checking ortho-para conversion, gas temperature (K) - " << tgrid[j] << endl
		<< "Ortho to para conversion, rate [cm3/s] - " << o << endl
		<< "Para to ortho conversion, rate [cm3/s] - " << p << endl;

	for (v = 0; v < nb_lev; v++) {
		cout << left << setw(5) << h2_di->lev_array[v].v << setw(5) << h2_di->lev_array[v].j << setw(12) << rates[v] << endl;
	}

	delete [] tgrid;
	delete [] rates;
	free_2d_array(coeff);
}
